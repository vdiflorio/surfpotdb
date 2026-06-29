#!/usr/bin/env python3
"""
ngpb_potential - reconstruct the NextGenPB electrostatic potential (and field) at
ANY point inside a structure, starting from the surface VTP enriched with the
per-node polarization charge `q_pol`.

It reproduces the `pot_field_fast` routine of the NextGenPB solver by summing
contributions defined ON THE NODES of the triangulated surface:

    phi(r) = phi_c(r) + phi_p(r) + phi_i(r)

    phi_c(r) = (1/eps_in) * sum_a  q_a / |r - r_a|              (Coulomb, atoms)
    phi_p(r) = C_pol      * sum_V  q_pol(V) / |r - V|           (polarization)
    phi_i(r) = (1/4pi)    * sum_V  phi(V) * (V-r).N(V)/|V-r|^3 * W_V
                         - C_react * sum_V q_pol(V) / |r - V|   (ionic / reaction)

where:
    eps_in/eps_out = 4pi e_0 eps_r kb T Angs / e^2
    C_pol   = (1/eps_out - 1/eps_in)/(4pi)
    C_react = (1/eps_out)/(4pi)
    W_V     = 1/3 * sum of the areas of the triangles incident to node V

The electric field E is obtained from the same sums differentiated w.r.t. r.
V, phi(V), N(V), q_pol(V) and the triangles are read from the enriched VTP;
q_a, r_a from the PQR file. The potential is in reduced units kT/e (the same as
the VTP `phi` array): multiply by kt_over_e_mV(T) to get millivolts.

Dependencies: numpy, vtk.

CLI usage:
    python ngpb_potential.py structure.vtp mostocc.pqr            # -> on atoms
    python ngpb_potential.py structure.vtp mostocc.pqr --points pts.xyz
    python ngpb_potential.py structure.vtp mostocc.pqr -o out.dat --field

Library usage:
    from ngpb_potential import SurfacePotential
    sp = SurfacePotential.from_vtp("structure.vtp", pqr="mostocc.pqr")
    phi, E = sp.on_atoms()              # default: atom centers
    phi, E = sp.evaluate(points)        # arbitrary points (N,3)
"""

import argparse
import numpy as np

# -- Physical constants (runner/NextGenPB/include/pb_class.h) ------------------
E_0   = 8.85418781762e-12   # vacuum permittivity [F/m]
KB    = 1.380649e-23        # Boltzmann constant [J/K]
ECH   = 1.602176634e-19     # elementary charge [C]
N_AV  = 6.022e23            # Avogadro number [1/mol]
ANGS  = 1e-10               # Angstrom [m]
PI    = 3.14159265358979323846
INV_4PI = 1.0 / (4.0 * PI)


def kt_over_e_mV(T: float) -> float:
    """kT/e in millivolts (to convert reduced phi -> mV)."""
    return KB * T / ECH * 1.0e3


class SurfacePotential:
    """Reconstructs phi/E from an enriched surface VTP + PQR charges."""

    def __init__(self, nodes, phi_surf, normals, q_pol, w_area,
                 atom_pos, atom_q, *, eps_in_r=2.0, eps_out_r=80.0,
                 T=298.15, ionic_strength=0.145):
        self.V = np.ascontiguousarray(nodes, dtype=np.float64)        # (Nv,3)
        self.phi_s = np.ascontiguousarray(phi_surf, dtype=np.float64) # (Nv,)
        self.N = np.ascontiguousarray(normals, dtype=np.float64)      # (Nv,3)
        self.q_pol = np.ascontiguousarray(q_pol, dtype=np.float64)    # (Nv,)
        self.W = np.ascontiguousarray(w_area, dtype=np.float64)       # (Nv,)
        self.atom_pos = np.ascontiguousarray(atom_pos, dtype=np.float64)
        self.atom_q = np.ascontiguousarray(atom_q, dtype=np.float64)
        self.T = T

        eps_in  = 4.0 * PI * E_0 * eps_in_r  * KB * T * ANGS / (ECH * ECH)
        eps_out = 4.0 * PI * E_0 * eps_out_r * KB * T * ANGS / (ECH * ECH)
        self.den_in  = 1.0 / eps_in
        self.C_pol   = (1.0 / eps_out - 1.0 / eps_in) * INV_4PI
        self.C_react = (1.0 / eps_out) * INV_4PI

        C0 = 1.0e3 * N_AV * ionic_strength
        k2 = 2.0 * C0 * ANGS * ANGS * ECH * ECH / (E_0 * eps_out_r * KB * T)
        self.has_salt = np.sqrt(k2) > 1.0e-5

    # -- construction from the enriched VTP ------------------------------------
    @classmethod
    def from_vtp(cls, vtp_path, pqr, **kw):
        import vtk
        try:
            from vtkmodules.util.numpy_support import vtk_to_numpy
        except ImportError:
            from vtk.util.numpy_support import vtk_to_numpy

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(str(vtp_path))
        reader.Update()
        pd = reader.GetOutput()

        nodes = vtk_to_numpy(pd.GetPoints().GetData()).astype(np.float64)
        pdata = pd.GetPointData()
        phi_s = vtk_to_numpy(pdata.GetArray("phi")).astype(np.float64)
        qpol_arr = pdata.GetArray("q_pol")
        if qpol_arr is None:
            raise ValueError(
                "The VTP does not contain the 'q_pol' array. Use an enriched "
                "VTP (scripts/enrich_vtp.py).")
        q_pol = vtk_to_numpy(qpol_arr).astype(np.float64)

        normals_arr = pdata.GetArray("Normals") or pdata.GetNormals()
        if normals_arr is None:
            raise ValueError("The VTP does not contain normals ('Normals').")
        normals = vtk_to_numpy(normals_arr).astype(np.float64)

        # triangle connectivity -> per-node effective area W_V
        polys = vtk_to_numpy(pd.GetPolys().GetData())
        tris = polys.reshape(-1, 4)[:, 1:]   # [3,i,j,k] -> i,j,k (triangle mesh)
        w_area = cls._vertex_areas(nodes, tris)

        atom_pos, atom_q = cls._read_pqr(pqr)
        return cls(nodes, phi_s, normals, q_pol, w_area, atom_pos, atom_q, **kw)

    @staticmethod
    def _vertex_areas(nodes, tris):
        a = nodes[tris[:, 0]]
        b = nodes[tris[:, 1]]
        c = nodes[tris[:, 2]]
        area = 0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1)  # (Ntri,)
        w = np.zeros(nodes.shape[0], dtype=np.float64)
        third = area / 3.0
        np.add.at(w, tris[:, 0], third)
        np.add.at(w, tris[:, 1], third)
        np.add.at(w, tris[:, 2], third)
        return w

    @staticmethod
    def _read_pqr(path):
        pos, q = [], []
        with open(path) as fh:
            for line in fh:
                if line.startswith(("ATOM", "HETATM")):
                    t = line.split()
                    # ... x y z charge radius  -> last 5 fields
                    x, y, z, charge = (float(t[-5]), float(t[-4]),
                                       float(t[-3]), float(t[-2]))
                    pos.append((x, y, z))
                    q.append(charge)
        return np.array(pos, dtype=np.float64), np.array(q, dtype=np.float64)

    # -- evaluation ------------------------------------------------------------
    def evaluate(self, points, field=True, block=64):
        """phi (and optionally E) at `points` (N,3).

        Points that coincide EXACTLY with an atom exclude the Coulomb self-term
        (as in pot_field_fast): convenient for evaluating on atoms.
        """
        pts = np.atleast_2d(np.asarray(points, dtype=np.float64))
        n = pts.shape[0]
        phi = np.zeros(n)
        E = np.zeros((n, 3)) if field else None

        for s in range(0, n, block):
            r = pts[s:s + block]
            self._accumulate_coulomb(r, phi, E, s)
            self._accumulate_surface(r, phi, E, s)
        return (phi, E) if field else phi

    def on_atoms(self, field=True, block=64):
        """Default: phi (and E) at the PQR atom centers."""
        return self.evaluate(self.atom_pos, field=field, block=block)

    def components(self, points=None, block=64):
        """Return the three separate contributions (phi_c, phi_p, phi_i) - handy
        for comparison/validation against pot_field.dat. Default: on atoms."""
        pts = self.atom_pos if points is None else \
            np.atleast_2d(np.asarray(points, dtype=np.float64))
        n = pts.shape[0]
        phi_c = np.zeros(n); phi_p = np.zeros(n); phi_i = np.zeros(n)
        for s in range(0, n, block):
            r = pts[s:s + block]; b = r.shape[0]
            # Coulomb (self-term excluded where r=0)
            d = r[:, None, :] - self.atom_pos[None, :, :]
            r2 = np.einsum("ijk,ijk->ij", d, d)
            inv_r = np.zeros_like(r2); nz = r2 > 0.0
            inv_r[nz] = 1.0 / np.sqrt(r2[nz])
            phi_c[s:s + b] = self.den_in * (inv_r * self.atom_q[None, :]).sum(1)
            # surface
            d_rV = r[:, None, :] - self.V[None, :, :]
            inv_rv = 1.0 / np.sqrt(np.einsum("ijk,ijk->ij", d_rV, d_rV))
            sum_qr = (self.q_pol[None, :] * inv_rv).sum(1)
            phi_p[s:s + b] = self.C_pol * sum_qr
            phi_i[s:s + b] = -self.C_react * sum_qr
            if self.has_salt:
                inv_rv3 = inv_rv ** 3
                dot = np.einsum("ijk,jk->ij", -d_rV, self.N)   # (V-r).N
                geo = self.phi_s[None, :] * dot * inv_rv3 * self.W[None, :]
                phi_i[s:s + b] += INV_4PI * geo.sum(1)
        return phi_c, phi_p, phi_i

    def _accumulate_coulomb(self, r, phi, E, off):
        d = r[:, None, :] - self.atom_pos[None, :, :]     # (b,Na,3) = r - r_a
        r2 = np.einsum("ijk,ijk->ij", d, d)
        inv_r = np.zeros_like(r2)
        nz = r2 > 0.0                                     # excludes self (r=0)
        inv_r[nz] = 1.0 / np.sqrt(r2[nz])
        phi[off:off + r.shape[0]] += self.den_in * (inv_r * self.atom_q[None, :]).sum(1)
        if E is not None:
            inv_r3 = inv_r ** 3
            coef = (self.atom_q[None, :] * inv_r3)[:, :, None]
            E[off:off + r.shape[0]] += self.den_in * (coef * d).sum(1)

    def _accumulate_surface(self, r, phi, E, off):
        b = r.shape[0]
        d_rV = r[:, None, :] - self.V[None, :, :]         # (b,Nv,3) = r - V
        dV = -d_rV                                        # V - r
        r2 = np.einsum("ijk,ijk->ij", d_rV, d_rV)
        inv_r = 1.0 / np.sqrt(r2)
        inv_r3 = inv_r ** 3

        sum_qr = (self.q_pol[None, :] * inv_r).sum(1)
        phi_p = self.C_pol * sum_qr
        phi_i_corr = -self.C_react * sum_qr

        if self.has_salt:
            dot = np.einsum("ijk,jk->ij", dV, self.N)     # (V-r).N
            geo = self.phi_s[None, :] * dot * inv_r3 * self.W[None, :]
            phi_i_geo = INV_4PI * geo.sum(1)
        else:
            phi_i_geo = 0.0

        phi[off:off + b] += phi_p + phi_i_geo + phi_i_corr

        if E is not None:
            qc = (self.q_pol[None, :] * inv_r3)[:, :, None]
            # E_p = C_pol sum q_pol (r-V)/r^3 ; ionic correction -C_react * same
            E[off:off + b] += (self.C_pol - self.C_react) * (qc * d_rV).sum(1)
            if self.has_salt:
                inv_r5 = inv_r3 * (inv_r ** 2)
                factor = (self.phi_s * self.W)[None, :]            # (1,Nv)
                dot = np.einsum("ijk,jk->ij", dV, self.N)          # (b,Nv)
                term = (-3.0 * dV * (dot * inv_r5)[:, :, None]
                        + self.N[None, :, :] * inv_r3[:, :, None])
                E[off:off + b] += INV_4PI * (factor[:, :, None] * term).sum(1)


def _read_points(path):
    pts = []
    with open(path) as fh:
        for line in fh:
            t = line.split()
            if len(t) >= 3:
                try:
                    pts.append((float(t[-3]), float(t[-2]), float(t[-1])))
                except ValueError:
                    continue
    return np.array(pts, dtype=np.float64)


def main():
    ap = argparse.ArgumentParser(
        description="NextGenPB electrostatic potential/field at interior points.")
    ap.add_argument("vtp", help="enriched VTP (with q_pol)")
    ap.add_argument("pqr", help="PQR file (atom charges/positions)")
    ap.add_argument("--points", help="file with x y z points (default: atoms)")
    ap.add_argument("-o", "--out", help="output file (default: stdout)")
    ap.add_argument("--field", action="store_true", help="also output field E")
    ap.add_argument("--eps-in", type=float, default=2.0)
    ap.add_argument("--eps-out", type=float, default=80.0)
    ap.add_argument("--T", type=float, default=298.15)
    ap.add_argument("--ionic-strength", type=float, default=0.145)
    args = ap.parse_args()

    sp = SurfacePotential.from_vtp(
        args.vtp, pqr=args.pqr, eps_in_r=args.eps_in, eps_out_r=args.eps_out,
        T=args.T, ionic_strength=args.ionic_strength)

    pts = _read_points(args.points) if args.points else sp.atom_pos
    out = sp.evaluate(pts, field=args.field)
    phi, E = (out if args.field else (out, None))

    import sys
    fh = open(args.out, "w") if args.out else sys.stdout
    hdr = "# index  x  y  z  phi[kT/e]" + ("  Ex  Ey  Ez" if args.field else "")
    fh.write(hdr + "\n")
    for i in range(pts.shape[0]):
        row = f"{i+1:6d}  {pts[i,0]:.4f}  {pts[i,1]:.4f}  {pts[i,2]:.4f}  {phi[i]:.8e}"
        if args.field:
            row += f"  {E[i,0]:.8e}  {E[i,1]:.8e}  {E[i,2]:.8e}"
        fh.write(row + "\n")
    if args.out:
        fh.close()
        print(f"Wrote {pts.shape[0]} points to {args.out}")


if __name__ == "__main__":
    main()
