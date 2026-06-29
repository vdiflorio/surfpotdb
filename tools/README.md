# ngpb_potential — electrostatic potential at interior points

Reconstructs the electrostatic potential (and electric field) computed by
[NextGenPB](https://github.com/vdiflorio/NextGenPB) at **any point inside** a
structure, starting from the surface VTP file **enriched** with the per-node
polarization charge `q_pol`.

The published VTP files contain, for each node of the triangulated surface:
`phi` (surface potential), `Normals`, and `q_pol` (polarization charge,
= flux of **D** through the surface). With this data plus the atomic charges
(PQR file) the `pot_field_fast` routine of the solver is reproduced exactly.

## Physics

The potential at a point `r` is the sum of three contributions, all defined
**on the nodes** of the surface:

```
phi(r) = phi_c(r) + phi_p(r) + phi_i(r)

phi_c(r) = (1/eps_in) · Σ_a q_a / |r − r_a|                       Coulomb (atoms)
phi_p(r) = C_pol      · Σ_V q_pol(V) / |r − V|                    polarization
phi_i(r) = (1/4π)     · Σ_V phi(V)·(V−r)·N(V)/|V−r|³ · W_V
                      − C_react · Σ_V q_pol(V) / |r − V|          ionic / reaction
```

with `eps = 4π·e₀·εr·kB·T·Å/e²`, `C_pol = (1/eps_out − 1/eps_in)/(4π)`,
`C_react = (1/eps_out)/(4π)`, and `W_V` = 1/3 of the sum of the areas of the
triangles incident to node `V` (computed from the mesh connectivity).

The ionic term is included only when the ionic strength is > 0.

**Units**: `phi` is in reduced units `kT/e` (same as the VTP `phi` array).
To convert to millivolts multiply by `kT/e ≈ 25.7 mV` at 298.15 K
(`ngpb_potential.kt_over_e_mV(T)`).

## Dependencies

```
pip install numpy vtk
```

## CLI usage

```bash
# Potential at atom centers (default), output to stdout
python ngpb_potential.py structure.vtp mostocc.pqr

# With electric field, to a file
python ngpb_potential.py structure.vtp mostocc.pqr --field -o out.dat

# At arbitrary points (file with "x y z" lines)
python ngpb_potential.py structure.vtp mostocc.pqr --points points.xyz
```

Model parameters (defaults = the dataset values): `--eps-in 2.0 --eps-out 80.0
--T 298.15 --ionic-strength 0.145`.

## Library usage

```python
from ngpb_potential import SurfacePotential

sp = SurfacePotential.from_vtp("structure.vtp", pqr="mostocc.pqr")

phi, E = sp.on_atoms()            # default: atom centers -> phi (Na,), E (Na,3)
phi, E = sp.evaluate(points)      # arbitrary points (N,3)
phi_c, phi_p, phi_i = sp.components(points)   # separate contributions
```

Points that coincide exactly with an atom automatically exclude the Coulomb
self-term (as in `pot_field_fast`).

## Input files

- **Enriched VTP**: must contain the `q_pol` array. Plain VTP files (only `phi`)
  are not enough. The enrichment is produced by `scripts/enrich_vtp.py` from
  `vertexdata.csv`.
- **PQR**: atom positions and charges (e.g. `mostocc.pqr`, the MCCE protonated
  structure).
