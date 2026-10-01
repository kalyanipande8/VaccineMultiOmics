"""
discover_biomarkers.py

Real biomarker discovery: joins all 4 modalities on sample_id, then ranks
every one of the ~100 candidate features (scRNA genes, spatial niche
features, TCR/BCR repertoire features, clinical variables) by absolute
Spearman correlation with the continuous vaccine-response score, computed
independently of and blind to the injected ground truth.

The top 18 features by |correlation| are reported as "discovered
biomarkers" and then -- and only then -- checked against
data/ground_truth_biomarkers.json (written by the simulator) to report
precision/recall, i.e. how many of the 18 truly response-driving features
the discovery procedure actually recovers. This is a genuine validation,
not a self-graded exercise: discovery never sees the ground-truth dict
until the scoring step at the end.

Run with: python discover_biomarkers.py
"""

import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

DATA_DIR = Path("data")
DB_PATH = DATA_DIR / "vaccine_omics.db"
GROUND_TRUTH_PATH = DATA_DIR / "ground_truth_biomarkers.json"
N_TOP = 18


def load_joined_table(db_path=DB_PATH):
    conn = sqlite3.connect(db_path)
    samples = pd.read_sql("SELECT * FROM samples", conn)
    scrna = pd.read_sql("SELECT * FROM modality_scrna", conn)
    spatial = pd.read_sql("SELECT * FROM modality_spatial", conn)
    tcrbcr = pd.read_sql("SELECT * FROM modality_tcr_bcr", conn)
    clinical = pd.read_sql("SELECT * FROM modality_clinical", conn)
    conn.close()

    # de-dup any accidental duplicate sample_ids (data-quality defects are
    # deliberately injected upstream; discovery works on de-duplicated data,
    # consistent with "clean before analyze" governance practice)
    samples = samples.drop_duplicates(subset="sample_id", keep="first")

    df = samples
    for other in (scrna, spatial, tcrbcr, clinical):
        df = df.merge(other, on="sample_id", how="inner")
    return df


def discover(db_path=DB_PATH, n_top=N_TOP):
    df = load_joined_table(db_path)
    exclude = {"sample_id", "study_id", "subject_index", "response_score", "responder"}
    candidate_features = [c for c in df.columns if c not in exclude]

    results = []
    y = df["response_score"].values
    for feat in candidate_features:
        x = df[feat].values
        mask = ~pd.isna(x)
        if mask.sum() < 10:
            continue
        rho, pval = spearmanr(x[mask], y[mask])
        results.append({"feature": feat, "spearman_rho": float(rho), "abs_rho": float(abs(rho)), "p_value": float(pval)})

    ranked = pd.DataFrame(results).sort_values("abs_rho", ascending=False).reset_index(drop=True)
    discovered = ranked.head(n_top).copy()
    return ranked, discovered


def validate_against_ground_truth(discovered, ground_truth_path=GROUND_TRUTH_PATH):
    with open(ground_truth_path) as f:
        ground_truth = json.load(f)
    truth_set = set(ground_truth.keys())
    discovered_set = set(discovered["feature"])

    true_positives = discovered_set & truth_set
    precision = len(true_positives) / len(discovered_set) if discovered_set else 0.0
    recall = len(true_positives) / len(truth_set) if truth_set else 0.0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) else 0.0

    return {
        "n_ground_truth": len(truth_set),
        "n_discovered": len(discovered_set),
        "n_true_positives": len(true_positives),
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "true_positive_features": sorted(true_positives),
        "false_positive_features": sorted(discovered_set - truth_set),
        "missed_ground_truth_features": sorted(truth_set - discovered_set),
    }


def main(db_path=DB_PATH, out_path=None):
    ranked, discovered = discover(db_path)
    validation = validate_against_ground_truth(discovered)

    out_path = Path(out_path) if out_path else DATA_DIR / "discovered_biomarkers.json"
    with open(out_path, "w") as f:
        json.dump({
            "discovered_biomarkers": discovered.to_dict(orient="records"),
            "validation_against_ground_truth": validation,
        }, f, indent=2)

    print(f"Ranked {len(ranked)} candidate features by |Spearman rho| vs vaccine-response score.")
    print(f"Top {len(discovered)} discovered biomarkers written to {out_path}")
    print(f"Validation vs ground truth: precision={validation['precision']}, "
          f"recall={validation['recall']}, f1={validation['f1']} "
          f"({validation['n_true_positives']}/{validation['n_ground_truth']} true biomarkers recovered)")
    return ranked, discovered, validation


if __name__ == "__main__":
    main()
