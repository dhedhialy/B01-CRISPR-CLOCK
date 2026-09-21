# Source Lock — Public Lineage-Recording Datasets

**Status:** LOCKED  
**Package mirrors:** `b01_clock.data.adapters.convexml`, `b01_clock.data.provenance`, GEO/manual loaders  
**Rule:** Headline external-and-wet-lab certificates cite only sources below. New sources require an amendment + audit checklist pass (`b01_clock audit`).

---

## 1. Primary benchmark — ConvexML

| Field | Lock |
|-------|------|
| **Name** | ConvexML lineage-recording benchmark / processed trees |
| **Accession** | Dryad DOI **`10.5061/dryad.qrfj6q5nz`** |
| **Landing** | https://doi.org/10.5061/dryad.qrfj6q5nz |
| **API** | `https://datadryad.org/api/v2/datasets/doi:10.5061%2Fdryad.qrfj6q5nz` |
| **Local path** | `data_store/convexml/` (auto) or `data_store/convexml/manual/` (hand-placed) |
| **Role** | **Primary external benchmark** for timing certificates; first external transfer after synthetic scout. **Simulated** (Cassiopeia-generated per deposit README) — carries known ground-truth times for coverage checks; must never be labeled wet-lab |
| **What we use** | Processed character matrices + associated trees (or tree-compatible leaf labels); edit/unedited/missing encodings mapped to `{0,1,-1}` |
| **What we do not use** | ConvexML’s own chronogram / convex-optimization **point timing** as ground truth for absolute ages; their outputs may appear only as **comparators**, not as oracle \(T\)-scaled truth |
| **Provenance** | SHA-256 of downloaded zip / manual files via `ProvenanceRecord`; store `DOWNLOAD_STATUS.json` |

**Scientific use note.** ConvexML targets scalable reconstruction under recorder models. B01 reuses the **same simulated public benchmark** to ask a different question: which timings remain in the profile-feasible set once rate bounds and dated \(T\) are locked — and, because the deposit is simulated, whether those intervals **cover the known ground truth** (they do not: 1/9; see `07_REAL_DATA/convexml_ground_truth_coverage.json`).

---

## 2. Sequential-recorder transfer — SciPhy

| Field | Lock |
|-------|------|
| **Name** | SciPhy (scientific phylogenetics / sequential lineage recorder cohort) |
| **Accession** | GEO **`GSE315827`** |
| **Landing** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE315827 |
| **Local path** | `data_store/sciphy/` |
| **Role** | **Transfer** test: sequential / multi-state or ordered editing architectures beyond binary irreversible scout defaults |
| **What we use** | Public processed count matrices / allele tables that can be collapsed or mapped to irreversible exposure summaries compatible with \(p=1-e^{-rt}\) (or documented multi-state extension in sensitivity) |
| **What we do not use** | Unpublished supplemental chronologies; any timing labels not in the GEO release |
| **Analysis posture** | Secondary: report \(\mathcal{I}_v\) and identification status; do not pool with ConvexML for a single “mega-interval” |

---

## 3. External architecture robustness — PALINCODE

| Field | Lock |
|-------|------|
| **Name** | PALINCODE lineage recording |
| **Accession** | GEO **`GSE327634`** |
| **Landing** | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE327634 |
| **Local path** | `data_store/palincode/` |
| **Role** | **Architecture robustness**: recorder design / barcode chemistry distinct from ConvexML and SciPhy |
| **What we use** | Public barcode/allele observation tables mappable to leaf×site states |
| **What we do not use** | Vendor-internal pipelines as silent preprocessing without a frozen script hash |
| **Analysis posture** | Pre-registered external stress test; failure here is **sensitivity**, not automatic kill of ConvexML results (document in referee table) |

---

## 4. Wet-lab certificate stress — intMEMOIR

| Field | Lock |
|-------|------|
| **Name** | intMEMOIR (integrase / MEMOIR-family lineage recorder) |
| **Accession** | Public processed trees + character matrices as released by the originating study (cite versioned archive in provenance JSON at download time) |
| **Local path** | `data_store/intmemoir/` |
| **Role** | Wet-lab certificate stress on an established imaging / integrase recorder lineage with public trees (only certified wet-lab tier once integrated) |
| **What we use** | Leaf barcodes, site edit states, and published tree topologies used as **audited topology inputs** (not re-estimated in headline unless topology ensemble is explicitly on) |
| **What we do not use** | Informal “known developmental ages” as soft constraints unless entered as an explicit dated-anchor amendment |
| **Analysis posture** | Demonstrate certificate emission on non-synthetic trees; compare identification regime to scout |

---

## 5. Source hierarchy (locked)

```
Synthetic exact scout  →  ConvexML (primary)  →  intMEMOIR (real trees)
                                              →  SciPhy (sequential transfer)
                                              →  PALINCODE (architecture robustness)
```

Headline wet-lab figure defaults to ConvexML **only if** certified; until then the external tier is the simulated ConvexML benchmark alone (coverage 1/9, reported in manuscript §5.4). SciPhy / PALINCODE / intMEMOIR appear as transfer / robustness panels unless a source fails audit.

---

## 6. Shared ingestion rules

1. **States:** `UNEDITED=0`, `EDITED=1`, `MISSING=-1` (`b01_clock.model.observations`).
2. **Irreversibility:** ancestral preparation forbids \(1\to 0\); such transitions yield \(-\infty\) log-likelihood.
3. **Dated \(T\):** must come from experimental metadata (hours / generations) or be set to \(T=1\) with all reports in units of \(T\); never invent absolute hours.
4. **Trees:** Newick or edge lists; leaf names must match matrix keys after a frozen alias map.
5. **Provenance minimum:** DOI/GEO accession, download UTC, file SHA-256, adapter version, alias map hash.
6. **Manual fallback:** If Dryad/GEO API fails, place files under `*/manual/` and record `DOWNLOAD_STATUS.json` with `ok: false` and operator note — still valid if hashes match published files.

---

## 7. Explicitly out of lock

- Private lab datasets without accession.
- Simulated data labeled as “real.”
- Reconstructed timings from third-party clock software used as **truth**.
- Any source requiring paywalled raw reads without a public processed release path.

---

## 8. Amendment protocol

To add a source: update this file, add adapter under `b01_clock/data/adapters/`, extend audit checklist, re-run holdout lock on that source’s site inventory, and add a row to `11_MANUSCRIPT/REFEREE_TABLE.md` for dataset-shift attacks.
