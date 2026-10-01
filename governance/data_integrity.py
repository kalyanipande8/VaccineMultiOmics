"""
data_integrity.py

Computes a real, auditable overall data-integrity percentage across the
consolidated multi-omics SQLite store. This is not a single arbitrary
number -- it's the mean pass-rate of four concrete, independently
checkable dimensions, each scored on the full dataset:

  1. Completeness       - fraction of all cell values that are non-null
                           across every modality + clinical table.
  2. Range validity      - fraction of clinical values that fall inside a
                           physiologically plausible range (catches
                           data-entry defects like age=999).
  3. Uniqueness          - fraction of sample_id values that are not
                           duplicated within samples.
  4. Referential integrity - fraction of samples that have a matching row
                           in every one of the 4 modality tables (no
                           orphaned / missing modality data for a sample).

Each dimension is reported individually in the compliance report so the
headline percentage is traceable back to specific, reproducible checks
(an ALCOA+ "accurate" + "complete" control), rather than asserted.
"""

import sqlite3
from pathlib import Path

import pandas as pd

DB_PATH = Path("data") / "vaccine_omics.db"

CLINICAL_RANGES = {
    "age": (0, 110), "bmi": (10, 60), "baseline_titer": (0, 200),
    "day7_AE_score": (0, 10), "day7_fever_c": (34, 42),
    "injection_site_reaction_score": (0, 10), "comorbidity_count": (0, 15),
    "prior_infection": (0, 1), "days_since_prior_vaccination": (0, 5000),
    "lymphocyte_count_baseline": (0, 10000),
}

MODALITY_TABLES = ["modality_scrna", "modality_spatial", "modality_tcr_bcr", "modality_clinical"]


def _load_tables(db_path=DB_PATH):
    conn = sqlite3.connect(db_path)
    tables = {}
    for name in ["studies", "samples", *MODALITY_TABLES]:
        tables[name] = pd.read_sql(f"SELECT * FROM {name}", conn)
    conn.close()
    return tables


def completeness_score(tables):
    total_cells, non_null = 0, 0
    for name in ["samples", *MODALITY_TABLES]:
        df = tables[name]
        total_cells += df.size
        non_null += df.notna().sum().sum()
    return non_null / total_cells if total_cells else 1.0, int(total_cells - non_null)


def range_validity_score(tables):
    df = tables["modality_clinical"]
    total, valid = 0, 0
    for col, (lo, hi) in CLINICAL_RANGES.items():
        if col not in df.columns:
            continue
        vals = df[col].dropna()
        total += len(vals)
        valid += int(((vals >= lo) & (vals <= hi)).sum())
    return (valid / total if total else 1.0), int(total - valid)


def uniqueness_score(tables):
    samples = tables["samples"]
    n = len(samples)
    n_unique = samples["sample_id"].nunique()
    dupes = n - n_unique
    return (n_unique / n if n else 1.0), int(dupes)


def referential_integrity_score(tables):
    sample_ids = set(tables["samples"]["sample_id"])
    n = len(sample_ids)
    if n == 0:
        return 1.0, 0
    complete = 0
    missing_links = 0
    for sid in sample_ids:
        present_everywhere = all(sid in set(tables[t]["sample_id"]) for t in MODALITY_TABLES)
        if present_everywhere:
            complete += 1
        else:
            missing_links += 1
    return complete / n, missing_links


def compute_integrity_report(db_path=DB_PATH):
    tables = _load_tables(db_path)

    comp_score, missing_cells = completeness_score(tables)
    range_score, bad_values = range_validity_score(tables)
    uniq_score, dupes = uniqueness_score(tables)
    ref_score, missing_links = referential_integrity_score(tables)

    dimensions = {
        "completeness": {"score": comp_score, "defects": missing_cells, "detail": "non-null cell fraction across samples + all modality tables"},
        "range_validity": {"score": range_score, "defects": bad_values, "detail": "clinical values within physiologically plausible bounds"},
        "uniqueness": {"score": uniq_score, "defects": dupes, "detail": "sample_id uniqueness in samples table"},
        "referential_integrity": {"score": ref_score, "defects": missing_links, "detail": "every sample has a row in all 4 modality tables"},
    }
    overall_pct = 100.0 * sum(d["score"] for d in dimensions.values()) / len(dimensions)

    report = {
        "overall_integrity_pct": round(overall_pct, 2),
        "dimensions": {k: {"score_pct": round(v["score"] * 100, 3), "defect_count": v["defects"], "detail": v["detail"]} for k, v in dimensions.items()},
        "n_samples": len(tables["samples"]),
        "n_studies": tables["studies"]["study_id"].nunique(),
    }
    return report


if __name__ == "__main__":
    import json
    print(json.dumps(compute_integrity_report(), indent=2))
