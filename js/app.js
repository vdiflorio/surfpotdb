/* NextGenPB surface electrostatics DB — static site logic.
 *
 * Reads two static JSON files (no server query):
 *   data/results.json  -> {generated, count, structures:[...]}
 *   data/resolve.json  -> {pdb_code: representative_pdb_code}  (members + reps)
 *
 * Cluster resolution (see WEB_DATABASE_PLAN §3):
 *   1. code is a computed representative          -> show it
 *   2. code is a cluster member (in resolve.json) -> show its representative + note
 *   3. otherwise                                  -> not found
 */

"use strict";

const DATA_BASE = "./data/"; // relative -> safe under GitHub Pages /repo/ subpath

// Populated on load:
let RESULTS = [];            // array of structure objects
let RESULTS_MAP = new Map(); // pdb_code -> structure
let RESOLVE = {};            // pdb_code -> representative pdb_code

// ---------- formatting helpers ----------
function fmt(value, maxFrac) {
  if (value === null || value === undefined || value === "") return "—";
  const n = Number(value);
  if (!isFinite(n)) return "—";
  return n.toLocaleString("en-US", { maximumFractionDigits: maxFrac });
}
const fmtE = (v) => fmt(v, 1);   // energies [kT]
const fmtPhi = (v) => fmt(v, 2); // potentials [kT/e]
const fmtRes = (v) => (v === null || v === undefined ? "—" : fmt(v, 2));

function escapeHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, (c) => (
    { "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]
  ));
}

// UniProt accession(s): show the first (linked to uniprot.org); any others go in
// the tooltip with a "+N" marker. Value is pipe-separated in the data.
function uniprotHtml(value) {
  const accs = String(value ?? "").split("|").filter(Boolean);
  if (!accs.length) return "";
  const first = accs[0];
  const title = accs.length > 1 ? ` title="${escapeHtml(accs.join(", "))}"` : "";
  const more = accs.length > 1 ? `<span class="uni-more">+${accs.length - 1}</span>` : "";
  return `<a class="uniprot" href="https://www.uniprot.org/uniprotkb/${encodeURIComponent(first)}"`
    + ` target="_blank" rel="noopener"${title}>${escapeHtml(first)}</a>${more}`;
}

// download link or "in arrivo" badge
function dlButton(url, label) {
  if (url) {
    return `<a class="btn-dl" href="${escapeHtml(url)}" target="_blank" rel="noopener">${label}</a>`;
  }
  return `<span class="badge-pending" title="Available after Zenodo upload">${label} — coming soon</span>`;
}

// ---------- resolver (lookup box) ----------
function resolve(rawCode) {
  const code = String(rawCode || "").trim().toLowerCase();
  const out = document.getElementById("resolve-result");
  out.hidden = false;

  if (!code) {
    out.innerHTML = `<div class="error">Please enter a PDB code.</div>`;
    return;
  }

  // 1. directly computed representative
  if (RESULTS_MAP.has(code)) {
    out.innerHTML = renderCard(RESULTS_MAP.get(code), null);
    return;
  }

  // 2. cluster member -> representative
  if (Object.prototype.hasOwnProperty.call(RESOLVE, code)) {
    const rep = RESOLVE[code];
    if (RESULTS_MAP.has(rep)) {
      const s = RESULTS_MAP.get(rep);
      const note = `<div class="note"><strong>${escapeHtml(code.toUpperCase())}</strong> `
        + `belongs to cluster <strong>${escapeHtml(s.cluster_id)}</strong>, `
        + `represented by <strong>${escapeHtml(rep.toUpperCase())}</strong>. `
        + `Showing the representative's results.</div>`;
      out.innerHTML = renderCard(s, note);
      return;
    }
    // representative known but not yet computed
    out.innerHTML = `<div class="note"><strong>${escapeHtml(code.toUpperCase())}</strong> `
      + `belongs to a cluster represented by <strong>${escapeHtml(rep.toUpperCase())}</strong>, `
      + `which has not been computed yet. Results not available.</div>`;
    return;
  }

  // 3. not found
  out.innerHTML = `<div class="error">PDB code <strong>${escapeHtml(code.toUpperCase())}</strong> `
    + `is not in the dataset (neither a computed structure nor a known cluster member).</div>`;
}

function metric(label, value, unit) {
  return `<div class="metric"><span class="label">${label}</span>`
    + `<span class="value">${value}${unit ? ` <span class="unit">${unit}</span>` : ""}</span></div>`;
}

// Human-readable protein class for a structure: GO functional buckets, plus the
// EC enzyme subclass(es) when present. Empty string if nothing is annotated.
function classText(s) {
  const fc = String(s.functional_class ?? "").split("|").filter(Boolean);
  const ec = String(s.ec_class ?? "").split("|").filter(Boolean)
    .map((d) => EC_LABELS[d] || `EC ${d}`);
  const parts = fc.slice();
  if (ec.length) parts.push(ec.join(", "));
  return parts.join(" · ");
}

function renderCard(s, noteHtml) {
  const proteinLine = [s.protein_name, s.release_date ? `released ${s.release_date}` : null]
    .filter(Boolean).map(escapeHtml).join(" · ");

  const uni = uniprotHtml(s.uniprot);
  const klass = classText(s);

  return (noteHtml || "")
    + `<h3 class="card-title"><span class="code">${escapeHtml(s.pdb_code)}</span>`
    +   (uni ? ` <span class="card-uniprot">${uni}</span>` : "")
    + `</h3>`
    + `<p class="card-sub">${proteinLine || "—"}</p>`
    + (klass ? `<p class="card-class"><span class="card-class-label">Class</span> ${escapeHtml(klass)}</p>` : "")
    + `<div class="metrics">`
    +   metric("Polarization", fmtE(s.polarization_energy), "kT")
    +   metric("Ionic", fmtE(s.ionic_energy), "kT")
    +   metric("Coulombic", fmtE(s.coulombic_energy), "kT")
    +   metric("Total", fmtE(s.total_energy), "kT")
    +   metric("Potential min", fmtPhi(s.surface_potential_min), "kT/e")
    +   metric("Potential max", fmtPhi(s.surface_potential_max), "kT/e")
    +   metric("Resolution", fmtRes(s.resolution), "Å")
    +   metric("Cluster", `${escapeHtml(s.cluster_id)}`, `size ${escapeHtml(s.cluster_size)}`)
    +   metric("Atoms", fmt(s.n_atoms, 0), "")
    +   metric("Residues", fmt(s.n_residues, 0), "")
    + `</div>`
    + `<div class="dl-row">`
    +   dlButton(s.vtp_url, "Download .vtp (ParaView)")
    +   dlButton(s.pqr_url, "Download .pqr (PyMOL/VMD)")
    + `</div>`;
}

// ---------- DataTables ----------
function numCol(render) {
  return { className: "num", render: render };
}

function miniDl(url, ext) {
  if (url) return `<a class="dl-mini" href="${escapeHtml(url)}" target="_blank" rel="noopener">${ext}</a>`;
  return `<span class="dl-mini disabled" title="coming soon">${ext}</span>`;
}

// ---------- numeric range filters ----------
const NUMERIC_FILTERS = [
  { key: "surface_potential_min", label: "Surface potential (min)", unit: "kT/e" },
  { key: "surface_potential_max", label: "Surface potential (max)", unit: "kT/e" },
  { key: "n_residues", label: "Residues", unit: "" },
  { key: "cluster_size", label: "Cluster size", unit: "" },
];
let FILTER_INPUTS = []; // [{key, min, max}]

// ---------- categorical class filter (functional_class / ec_class) ----------
// Each option encodes which row field it matches and the value to look for, e.g.
// {field:"functional_class", value:"Receptor"} or {field:"ec_class", value:"3"}.
// "__none__" matches rows whose functional_class is empty (no MF annotation).
const FUNCTIONAL_BUCKETS = [
  "Enzyme", "Receptor", "Transporter", "Transcription regulator",
  "Regulator/inhibitor", "Toxin", "Structural", "Electron transfer",
  "Antioxidant", "Translation regulator", "Molecular adaptor", "Binding",
  "Other function",
];
const EC_LABELS = {
  "1": "Oxidoreductases (EC 1)", "2": "Transferases (EC 2)",
  "3": "Hydrolases (EC 3)", "4": "Lyases (EC 4)", "5": "Isomerases (EC 5)",
  "6": "Ligases (EC 6)", "7": "Translocases (EC 7)",
};
let CLASS_FILTER = null; // {field, value} | null

// Active-filter badge (numeric inputs + categorical class). Shared by both
// filter builders so the count stays correct whichever control changes.
function updateFilterCount() {
  const countEl = document.getElementById("filter-count");
  if (!countEl) return;
  let n = FILTER_INPUTS.filter((f) => f.min.value !== "" || f.max.value !== "").length;
  if (CLASS_FILTER) n += 1;
  countEl.textContent = `${n} active`;
  countEl.hidden = n === 0;
}

function rowHasClass(rowData, field, value) {
  const raw = String(rowData[field] ?? "");
  if (field === "functional_class" && value === "__none__") return raw === "";
  if (raw === "") return false;
  return raw.split("|").includes(value);
}

function buildClassFilter(table) {
  const sel = document.getElementById("class-filter");
  if (!sel) return;

  // counts per option, computed from the loaded data
  const fcCount = new Map(), ecCount = new Map();
  let noneCount = 0;
  for (const r of RESULTS) {
    const fc = String(r.functional_class ?? "");
    if (fc === "") noneCount++;
    else fc.split("|").forEach((b) => fcCount.set(b, (fcCount.get(b) || 0) + 1));
    const ec = String(r.ec_class ?? "");
    if (ec) ec.split("|").forEach((d) => ecCount.set(d, (ecCount.get(d) || 0) + 1));
  }

  const mkOpt = (field, value, label, n) => {
    const o = document.createElement("option");
    o.value = `${field}:${value}`;
    o.textContent = n != null ? `${label} (${n.toLocaleString("en-US")})` : label;
    return o;
  };

  const gFunc = document.createElement("optgroup");
  gFunc.label = "By molecular function";
  for (const b of FUNCTIONAL_BUCKETS) {
    if (fcCount.get(b)) gFunc.append(mkOpt("functional_class", b, b, fcCount.get(b)));
  }
  sel.append(gFunc);

  const gEc = document.createElement("optgroup");
  gEc.label = "By enzyme class (EC)";
  for (const d of ["1", "2", "3", "4", "5", "6", "7"]) {
    if (ecCount.get(d)) gEc.append(mkOpt("ec_class", d, EC_LABELS[d], ecCount.get(d)));
  }
  sel.append(gEc);

  const gOther = document.createElement("optgroup");
  gOther.label = "Other";
  gOther.append(mkOpt("functional_class", "__none__", "No functional annotation", noneCount));
  sel.append(gOther);

  sel.addEventListener("change", () => {
    const v = sel.value;
    if (!v) { CLASS_FILTER = null; }
    else {
      const i = v.indexOf(":");
      CLASS_FILTER = { field: v.slice(0, i), value: v.slice(i + 1) };
    }
    table.draw();
    updateFilterCount();
  });
}

function buildFilters(table) {
  const container = document.getElementById("numeric-filters");
  const updateCount = updateFilterCount;

  FILTER_INPUTS = NUMERIC_FILTERS.map((f) => {
    const wrap = document.createElement("div");
    wrap.className = "filter-item";
    wrap.innerHTML = `<span class="filter-label">${f.label} <span class="unit">${f.unit}</span></span>`;

    const row = document.createElement("div");
    row.className = "filter-inputs";
    const min = document.createElement("input");
    const max = document.createElement("input");
    for (const [el, ph] of [[min, "min"], [max, "max"]]) {
      el.type = "number";
      el.step = "any";
      el.placeholder = ph;
      el.className = "filter-num";
      el.addEventListener("input", () => { table.draw(); updateCount(); });
    }
    row.append(min, max);
    wrap.append(row);
    container.append(wrap);
    return { key: f.key, min, max };
  });

  // custom range search: reads numeric values from the row's data object
  DataTable.ext.search.push((settings, searchData, dataIndex, rowData) => {
    if (CLASS_FILTER && !rowHasClass(rowData, CLASS_FILTER.field, CLASS_FILTER.value)) {
      return false;
    }
    for (const f of FILTER_INPUTS) {
      const v = Number(rowData[f.key]);
      const mn = f.min.value === "" ? NaN : parseFloat(f.min.value);
      const mx = f.max.value === "" ? NaN : parseFloat(f.max.value);
      if (!isNaN(mn) && (!isFinite(v) || v < mn)) return false;
      if (!isNaN(mx) && (!isFinite(v) || v > mx)) return false;
    }
    return true;
  });

  document.getElementById("filter-reset").addEventListener("click", () => {
    FILTER_INPUTS.forEach((f) => { f.min.value = ""; f.max.value = ""; });
    const sel = document.getElementById("class-filter");
    if (sel) sel.value = "";
    CLASS_FILTER = null;
    table.draw();
    updateCount();
  });
}

function initTable() {
  const table = new DataTable("#results-table", {
    data: RESULTS,
    deferRender: true,
    pageLength: 25,
    order: [[0, "asc"]],
    columns: [
      { // PDB code — clickable -> loads into resolver
        data: "pdb_code",
        render: (d, type) => {
          if (type === "display") {
            return `<a class="code-link" data-code="${escapeHtml(d)}">${escapeHtml(d)}</a>`;
          }
          return d;
        },
      },
      { // UniProt — first accession linked, rest in tooltip
        data: "uniprot",
        render: (d, type) => (type === "display" ? (uniprotHtml(d) || "—") : (d || "")),
      },
      { data: "protein_name", defaultContent: "" },
      Object.assign({ data: "polarization_energy" }, numCol((d, t) => (t === "display" ? fmtE(d) : d ?? 0))),
      Object.assign({ data: "ionic_energy" }, numCol((d, t) => (t === "display" ? fmtE(d) : d ?? 0))),
      Object.assign({ data: "coulombic_energy" }, numCol((d, t) => (t === "display" ? fmtE(d) : d ?? 0))),
      Object.assign({ data: "total_energy" }, numCol((d, t) => (t === "display" ? fmtE(d) : d ?? 0))),
      Object.assign({ data: "surface_potential_min" }, numCol((d, t) => (t === "display" ? fmtPhi(d) : d ?? 0))),
      Object.assign({ data: "surface_potential_max" }, numCol((d, t) => (t === "display" ? fmtPhi(d) : d ?? 0))),
      Object.assign({ data: "n_residues" }, numCol((d, t) => (t === "display" ? fmt(d, 0) : d ?? 0))),
      { // downloads
        data: null,
        orderable: false,
        searchable: false,
        render: (row) => `${miniDl(row.vtp_url, "vtp")} · ${miniDl(row.pqr_url, "pqr")}`,
      },
    ],
  });

  // click on a code in the table -> populate and run the resolver
  document.querySelector("#results-table").addEventListener("click", (ev) => {
    const a = ev.target.closest(".code-link");
    if (!a) return;
    const code = a.getAttribute("data-code");
    document.getElementById("pdb-input").value = code;
    resolve(code);
    document.getElementById("resolver-panel").scrollIntoView({ behavior: "smooth", block: "start" });
  });

  return table;
}

// ---------- bootstrap ----------
async function main() {
  const form = document.getElementById("resolve-form");
  form.addEventListener("submit", (ev) => {
    ev.preventDefault();
    resolve(document.getElementById("pdb-input").value);
  });

  try {
    const [results, resolveMap] = await Promise.all([
      fetch(DATA_BASE + "results.json").then((r) => {
        if (!r.ok) throw new Error("results.json: HTTP " + r.status);
        return r.json();
      }),
      fetch(DATA_BASE + "resolve.json").then((r) => {
        if (!r.ok) throw new Error("resolve.json: HTTP " + r.status);
        return r.json();
      }),
    ]);

    RESULTS = results.structures || [];
    RESOLVE = resolveMap || {};
    RESULTS_MAP = new Map(RESULTS.map((s) => [s.pdb_code, s]));

    document.getElementById("count").textContent = (results.count ?? RESULTS.length).toLocaleString("en-US");
    document.getElementById("generated").textContent = results.generated || "—";

    const table = initTable();
    buildClassFilter(table);
    buildFilters(table);
  } catch (err) {
    const out = document.getElementById("resolve-result");
    out.hidden = false;
    out.innerHTML = `<div class="error">Failed to load data: ${escapeHtml(err.message)}.<br>`
      + `If you opened the file directly (file://), start a local HTTP server instead — see README.</div>`;
    console.error(err);
  }
}

document.addEventListener("DOMContentLoaded", main);
