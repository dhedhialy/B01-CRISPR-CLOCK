# Timing Certificate JSON Schema

**Status:** LOCKED emission contract  
**Package:** `b01_clock.solver.certificates.CertificateReport` (extended fields below)  
**Version:** `b01-certificate-1.0`

---

## 1. Purpose

A **timing certificate** is a machine-readable record of a sharp profile interval \(\mathcal{I}_v(D)\) under locked constraints \(\mathcal{C}\), including identification status and provenance needed to audit leakage, invariance, and topology assumptions.

---

## 2. Top-level object

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://b01-clock.local/schemas/timing_certificate_v1.json",
  "title": "B01TimingCertificate",
  "type": "object",
  "required": [
    "schema_version",
    "method",
    "total_time",
    "confidence",
    "rate_bounds",
    "constraints",
    "split",
    "targets",
    "diagnostics"
  ],
  "properties": {
    "schema_version": { "const": "b01-certificate-1.0" },
    "method": {
      "type": "string",
      "enum": ["exact_profile_grid", "scalable_outer", "topology_ensemble_union"]
    },
    "created_utc": { "type": "string", "format": "date-time" },
    "seed": { "type": ["integer", "null"] },
    "total_time": { "type": "number", "exclusiveMinimum": 0 },
    "confidence": { "type": "number", "minimum": 0, "maximum": 1 },
    "topology_id": { "type": "string" },
    "tree_newick": { "type": "string" },
    "dataset": {
      "type": "object",
      "required": ["source_id"],
      "properties": {
        "source_id": {
          "type": "string",
          "enum": ["synthetic", "convexml", "sciphy", "palincode", "intmemoir", "other"]
        },
        "accession": { "type": "string" },
        "provenance_sha256": { "type": "string" },
        "adapter": { "type": "string" }
      }
    },
    "rate_bounds": {
      "type": "object",
      "required": ["mode", "global"],
      "properties": {
        "mode": { "type": "string", "enum": ["global", "per_site"] },
        "global": {
          "type": "object",
          "required": ["min", "max"],
          "properties": {
            "min": { "type": "number", "exclusiveMinimum": 0 },
            "max": { "type": "number", "exclusiveMinimum": 0 }
          }
        },
        "per_site": {
          "type": "array",
          "items": {
            "type": "object",
            "required": ["site", "min", "max"],
            "properties": {
              "site": { "type": "integer", "minimum": 0 },
              "min": { "type": "number" },
              "max": { "type": "number" }
            }
          }
        },
        "source": {
          "type": "string",
          "enum": ["calibration_holdout", "config_lock", "oracle_sim", "sensitivity"]
        }
      }
    },
    "constraints": {
      "type": "object",
      "required": ["ultrametric", "dated_anchor", "irreversible"],
      "properties": {
        "ultrametric": { "type": "boolean" },
        "dated_anchor": { "type": "boolean" },
        "irreversible": { "type": "boolean" },
        "rate_time_guardrail": { "type": "boolean" },
        "notes": { "type": "string" }
      }
    },
    "split": {
      "type": "object",
      "required": ["locked_before_headline", "calib_sites", "analysis_sites", "seed"],
      "properties": {
        "locked_before_headline": { "const": true },
        "calib_frac": { "type": "number" },
        "calib_sites": { "type": "array", "items": { "type": "integer" } },
        "analysis_sites": { "type": "array", "items": { "type": "integer" } },
        "seed": { "type": "integer" },
        "degenerate_split": { "type": "boolean" },
        "note": { "type": "string" }
      }
    },
    "targets": {
      "type": "array",
      "minItems": 1,
      "items": { "$ref": "#/$defs/target_certificate" }
    },
    "diagnostics": { "$ref": "#/$defs/diagnostics" },
    "notes": {
      "type": "array",
      "items": { "type": "string" }
    }
  },
  "$defs": {
    "identification_status": {
      "type": "string",
      "enum": [
        "unidentified",
        "partially_identified",
        "point_identified",
        "infeasible",
        "abstain"
      ]
    },
    "target_certificate": {
      "type": "object",
      "required": [
        "target",
        "lower",
        "upper",
        "mle",
        "status",
        "width_frac",
        "l_max",
        "threshold"
      ],
      "properties": {
        "target": {
          "type": "string",
          "description": "age:NODE | branch:P->C | order_gap:A,B"
        },
        "lower": { "type": "number" },
        "upper": { "type": "number" },
        "mle": { "type": "number" },
        "status": { "$ref": "#/$defs/identification_status" },
        "width_frac": {
          "type": "number",
          "description": "(upper-lower)/total_time for ages/branches; for gaps use (upper-lower)/total_time as scale proxy"
        },
        "l_max": { "type": "number" },
        "threshold": { "type": "number" },
        "n_feasible_grid": { "type": "integer" },
        "grid_size": { "type": "integer" },
        "grid_stable": { "type": "boolean" },
        "ensemble": {
          "type": "object",
          "properties": {
            "union_lower": { "type": "number" },
            "union_upper": { "type": "number" },
            "robust_order": { "type": "boolean" },
            "per_topology": {
              "type": "array",
              "items": {
                "type": "object",
                "required": ["topology", "lower", "upper", "status"],
                "properties": {
                  "topology": { "type": "string" },
                  "lower": { "type": "number" },
                  "upper": { "type": "number" },
                  "status": { "$ref": "#/$defs/identification_status" }
                }
              }
            }
          }
        },
        "certificate_statement": {
          "type": "string",
          "description": "Human-readable profile LR set description"
        }
      }
    },
    "diagnostics": {
      "type": "object",
      "required": ["invariance", "rate_time_guardrail"],
      "properties": {
        "invariance": {
          "type": "object",
          "required": ["ok", "scale"],
          "properties": {
            "ok": { "type": "boolean" },
            "scale": { "type": "number" },
            "ll_base": { "type": "number" },
            "ll_scaled": { "type": "number" },
            "abs_diff": { "type": "number" }
          }
        },
        "rate_time_guardrail": {
          "type": "object",
          "required": ["scale", "feasible_with_anchor_and_bounds"],
          "properties": {
            "scale": { "type": "number" },
            "rates_in_bounds": { "type": "boolean" },
            "preserves_dated_anchor": { "type": "boolean" },
            "feasible_with_anchor_and_bounds": { "type": "boolean" },
            "feasible_without_anchor": { "type": "boolean" },
            "guardrail": { "type": "string" }
          }
        },
        "coverage_ref": {
          "type": "string",
          "description": "Path or ID of coverage artifact used for calibration claim"
        },
        "false_identification_rate": { "type": "number", "minimum": 0, "maximum": 1 },
        "scout_verdict": { "type": "string", "enum": ["PASS", "KILL_OR_REFRAME", "N/A"] }
      }
    }
  }
}
```

---

## 3. Minimal valid example

```json
{
  "schema_version": "b01-certificate-1.0",
  "method": "exact_profile_grid",
  "created_utc": "2026-09-08T00:00:00Z",
  "seed": 7,
  "total_time": 1.0,
  "confidence": 0.95,
  "topology_id": "small_scout_tree",
  "dataset": { "source_id": "synthetic" },
  "rate_bounds": {
    "mode": "global",
    "global": { "min": 0.5, "max": 2.0 },
    "source": "oracle_sim"
  },
  "constraints": {
    "ultrametric": true,
    "dated_anchor": true,
    "irreversible": true,
    "rate_time_guardrail": true
  },
  "split": {
    "locked_before_headline": true,
    "calib_frac": 0.4,
    "calib_sites": [0, 1, 2],
    "analysis_sites": [3, 4, 5],
    "seed": 106
  },
  "targets": [
    {
      "target": "age:A",
      "lower": 0.0,
      "upper": 1.0,
      "mle": 0.5,
      "status": "partially_identified",
      "width_frac": 1.0,
      "l_max": 0.0,
      "threshold": -1.92,
      "certificate_statement": "Profile LR set with rates profiled in locked bounds; replace numeric fields from solver."
    }
  ],
  "diagnostics": {
    "invariance": { "ok": true, "scale": 2.0 },
    "rate_time_guardrail": {
      "scale": 2.0,
      "feasible_with_anchor_and_bounds": false
    },
    "scout_verdict": "N/A"
  },
  "notes": [
    "Example uses placeholder target numbers; emitters must write solver outputs, never invented intervals."
  ]
}
```

---

## 4. Emission rules

1. **No silent defaults:** `rate_bounds`, `split`, `confidence`, and `total_time` must appear explicitly.  
2. **Headline certificates** require `split.locked_before_headline === true` and `rate_bounds.source === "calibration_holdout"` (or documented `config_lock` for pure synthetic oracle bounds).  
3. **Status** must match `classify_identification(lower, upper, total_time)`.  
4. If invariance `ok !== true`, status of all targets becomes `abstain` and `notes` must cite failure.  
5. Ensemble method must fill `targets[].ensemble` or use top-level `method: topology_ensemble_union`.  
6. Do not put biological interpretation in `certificate_statement`; keep to statistical set description.

---

## 5. Mapping from Python

| Schema field | Source |
|--------------|--------|
| `targets[]` | `ExactSolverResult` / ensemble rows |
| `diagnostics.invariance` | `check_rate_time_invariance` |
| `diagnostics.rate_time_guardrail` | `rate_time_transform_feasible` |
| `split` | `lock_split_record` |
| `rate_bounds.global` | `held_out_rate_bounds` or config |
| top-level wrapper | `CertificateReport.to_dict` (extend to full schema in emitters) |
