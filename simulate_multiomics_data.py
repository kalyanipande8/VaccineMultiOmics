"""
simulate_multiomics_data.py

Simulates 3 molecular immune-profiling modalities (scRNA-seq, spatial
tissue profiling, TCR/BCR repertoire) plus clinical outcome data across
50+ vaccine studies, and consolidates everything into one unified SQLite
store (data/vaccine_omics.db).

Every consolidation step is logged to the data-governance lineage ledger
(governance/lineage.py) with source, transform, row counts, a content
hash, and a timestamp -- this is the backbone of the governance
framework that is the centerpiece of this project (see governance/).

A ground-truth set of exactly 18 features (spread across all 4 data
types) is wired, with known signed weights, into a simulated continuous
vaccine-response score (a neutralizing-antibody titer proxy). All other
features are pure noise, uncorrelated with response. Ground truth is
saved so biomarker discovery is validated against it rather than graded
on its own homework.

Run with: python simulate_multiomics_data.py
"""

import json
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd

from governance.lineage import LineageLedger

RNG_SEED = 42
N_STUDIES = 56
SUBJECTS_PER_STUDY = (8, 16)  # inclusive range, randomized per study

DATA_DIR = Path("data")
DB_PATH = DATA_DIR / "vaccine_omics.db"

# ---------------------------------------------------------------------------
# Feature definitions
# ---------------------------------------------------------------------------

SCRNA_GENES = [f"GENE_{g}" for g in [
    "IFNG", "TNF", "IL2", "IL6", "IL10", "IL21", "TBX21", "GATA3", "FOXP3",
    "GZMB", "PRF1", "CD38", "CD69", "CXCR5", "PDCD1", "ICOS", "BCL6", "AICDA",
    "MKI67", "CCR7", "SELL", "IL7R", "KLRG1", "TOX", "TCF7", "CD27", "CD28",
    "HLA_DR", "STAT1", "STAT3", "IRF4", "BATF", "XBP1", "PRDM1", "MS4A1",
    "CD19", "IGHG1", "IGHM", "JCHAIN", "TNFRSF17",
]]

SPATIAL_FEATURES = [
    "GC_density_LN", "Tfh_GC_colocalization", "Bcell_follicle_area",
    "Tcell_zone_density", "DC_Tcell_contact_score", "plasma_cell_medullary_cord",
    "niche_diversity_index", "macrophage_subcapsular_density",
    "Tfh_Bcell_distance_um", "germinal_center_count", "stromal_FDC_density",
    "CD8_infiltration_score", "neutrophil_margin_density",
    "vascular_density_LN", "niche_entropy",
]

TCR_BCR_FEATURES = [
    "TCR_clonality", "BCR_clonality", "TCR_shannon_diversity",
    "BCR_shannon_diversity", "top10_TCR_clone_freq", "top10_BCR_clone_freq",
    "vaccine_specific_clonotype_freq", "somatic_hypermutation_rate",
    "class_switch_fraction_IgG", "BCR_clone_expansion_index",
    "TCR_public_clonotype_count", "CDR3_length_mean_BCR",
    "convergent_recombination_score", "TCR_Vgene_skew", "BCR_Vgene_skew",
]

CLINICAL_FEATURES = [
    "age", "bmi", "baseline_titer", "day7_AE_score", "day7_fever_c",
    "injection_site_reaction_score", "comorbidity_count", "prior_infection",
    "days_since_prior_vaccination", "lymphocyte_count_baseline",
]

ALL_FEATURES = {
    "scrna": SCRNA_GENES,
    "spatial": SPATIAL_FEATURES,
    "tcr_bcr": TCR_BCR_FEATURES,
    "clinical": CLINICAL_FEATURES,
}

# ---------------------------------------------------------------------------
# Ground truth: exactly 18 predictive biomarkers, spread across modalities,
# each with a signed weight used to construct the response score.
# ---------------------------------------------------------------------------

GROUND_TRUTH_BIOMARKERS = {
    # scRNA-seq (6)
    "GENE_IFNG": 1.4, "GENE_TBX21": 1.1, "GENE_CXCR5": 0.9,
    "GENE_BCL6": 0.8, "GENE_MKI67": 0.7, "GENE_XBP1": 0.6,
    # spatial (4)
    "GC_density_LN": 1.3, "Tfh_GC_colocalization": 1.2,
    "germinal_center_count": 1.0, "stromal_FDC_density": 0.6,
    # TCR/BCR (4)
    "vaccine_specific_clonotype_freq": 1.6, "BCR_clonality": 1.0,
    "somatic_hypermutation_rate": 0.9, "class_switch_fraction_IgG": 0.7,
    # clinical (4)
    "baseline_titer": 1.1, "age": -0.8, "comorbidity_count": -0.6,
    "lymphocyte_count_baseline": 0.5,
}
assert len(GROUND_TRUTH_BIOMARKERS) == 18

VACCINE_TYPES = ["mRNA", "protein-subunit", "viral-vector", "inactivated", "live-attenuated"]
PLATFORMS_NOTE = "10x Genomics (scRNA), CODEX/Xenium (spatial), Adaptive immunoSEQ (TCR/BCR), EDC clinical capture"


def _feature_mu_sigma(feat, rng):
    """Approximate per-feature center/scale used only to standardize the
    ground-truth contribution to the response score at generation time."""
    if feat in SCRNA_GENES:
        return 5.0, 1.0
    if feat in SPATIAL_FEATURES:
        return 0.5, 0.15
    if feat in TCR_BCR_FEATURES:
        return 0.3, 0.1
    # clinical
    return {
        "age": (45.0, 15.0), "bmi": (26.0, 4.0), "baseline_titer": (20.0, 10.0),
        "day7_AE_score": (2.0, 1.0), "day7_fever_c": (37.0, 0.5),
        "injection_site_reaction_score": (1.5, 1.0), "comorbidity_count": (1.0, 1.0),
        "prior_infection": (0.3, 0.46), "days_since_prior_vaccination": (400, 200),
        "lymphocyte_count_baseline": (2000, 500),
    }[feat]


def simulate(seed=RNG_SEED, inject_missingness=True, db_path=None, ledger_path=None, run_id=None):
    db_path = Path(db_path) if db_path else DB_PATH
    rng = np.random.default_rng(seed)
    DATA_DIR.mkdir(exist_ok=True)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    ledger = LineageLedger(run_id=run_id)

    studies, samples = [], []
    scrna_rows, spatial_rows, tcrbcr_rows, clinical_rows = [], [], [], []

    for s in range(1, N_STUDIES + 1):
        study_id = f"VAX-{s:03d}"
        vaccine_type = rng.choice(VACCINE_TYPES)
        n_subjects = int(rng.integers(SUBJECTS_PER_STUDY[0], SUBJECTS_PER_STUDY[1] + 1))
        studies.append({
            "study_id": study_id,
            "vaccine_type": str(vaccine_type),
            "n_subjects": n_subjects,
            "phase": str(rng.choice(["Phase1", "Phase1/2", "Phase2"])),
            "year": int(rng.integers(2016, 2026)),
        })

        batch_shift = rng.normal(0, 0.15)  # study-level batch effect

        for subj in range(n_subjects):
            sample_id = f"{study_id}-S{subj+1:03d}"
            feature_values = {}

            scrna_vals = {g: rng.normal(5 + batch_shift, 1.0) for g in SCRNA_GENES}
            spatial_vals = {f: rng.normal(0.5 + batch_shift * 0.1, 0.15) for f in SPATIAL_FEATURES}
            tcrbcr_vals = {f: abs(rng.normal(0.3 + batch_shift * 0.05, 0.1)) for f in TCR_BCR_FEATURES}

            clinical_vals = {}
            for f in CLINICAL_FEATURES:
                mu, sigma = _feature_mu_sigma(f, rng)
                if f == "prior_infection":
                    clinical_vals[f] = float(rng.random() < mu)
                else:
                    clinical_vals[f] = abs(rng.normal(mu, sigma)) if f != "age" else rng.normal(mu, sigma)

            feature_values.update(scrna_vals)
            feature_values.update(spatial_vals)
            feature_values.update(tcrbcr_vals)
            feature_values.update(clinical_vals)

            z = 0.0
            for feat, weight in GROUND_TRUTH_BIOMARKERS.items():
                val = feature_values[feat]
                mu, sigma = _feature_mu_sigma(feat, rng)
                z += weight * (val - mu) / sigma
            noise = rng.normal(0, 1.8)
            response_score = z + noise
            response_pct = 100 / (1 + np.exp(-response_score / 3))
            responder = int(response_pct >= 50)

            samples.append({
                "sample_id": sample_id, "study_id": study_id,
                "subject_index": subj + 1,
                "response_score": round(float(response_pct), 3),
                "responder": responder,
            })
            scrna_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in scrna_vals.items()}})
            spatial_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in spatial_vals.items()}})
            tcrbcr_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in tcrbcr_vals.items()}})
            clinical_rows.append({"sample_id": sample_id, **{k: round(v, 4) for k, v in clinical_vals.items()}})

    studies_df = pd.DataFrame(studies)
    samples_df = pd.DataFrame(samples)
    scrna_df = pd.DataFrame(scrna_rows)
    spatial_df = pd.DataFrame(spatial_rows)
    tcrbcr_df = pd.DataFrame(tcrbcr_rows)
    clinical_df = pd.DataFrame(clinical_rows)

    # ------------------------------------------------------------------
    # Inject a small, known amount of realistic data-quality noise
    # (missing values + a few out-of-range/duplicate rows) so the
    # integrity check has real defects to find instead of a trivial 100%.
    # ------------------------------------------------------------------
    n_defects = 0
    if inject_missingness:
        for df, cols in [(clinical_df, CLINICAL_FEATURES), (scrna_df, SCRNA_GENES),
                          (spatial_df, SPATIAL_FEATURES), (tcrbcr_df, TCR_BCR_FEATURES)]:
            n_missing = max(1, int(0.004 * len(df) * len(cols)))
            for _ in range(n_missing):
                r = rng.integers(0, len(df))
                c = rng.choice(cols)
                df.loc[r, c] = np.nan
                n_defects += 1
        # out-of-range clinical values (data-entry errors)
        n_bad = max(1, int(0.01 * len(clinical_df)))
        for _ in range(n_bad):
            r = rng.integers(0, len(clinical_df))
            clinical_df.loc[r, "age"] = 999.0  # impossible age -> caught by range check
            n_defects += 1
        # a few duplicated sample rows (common real-world ingestion defect)
        n_dup = max(1, int(0.004 * len(samples_df)))
        if len(samples_df) > 1:
            dup = samples_df.iloc[rng.integers(0, len(samples_df), size=n_dup)].copy()
            samples_df = pd.concat([samples_df, dup], ignore_index=True)
            n_defects += n_dup
        # a few samples missing an entire modality row (ingestion dropout) ->
        # exercises the referential-integrity dimension honestly.
        n_orphan = max(1, int(0.004 * len(samples_df)))
        orphan_ids = rng.choice(samples_df["sample_id"].unique(), size=n_orphan, replace=False)
        tcrbcr_df = tcrbcr_df[~tcrbcr_df["sample_id"].isin(orphan_ids)].reset_index(drop=True)
        n_defects += n_orphan

    tables = {
        "studies": studies_df, "samples": samples_df, "modality_scrna": scrna_df,
        "modality_spatial": spatial_df, "modality_tcr_bcr": tcrbcr_df,
        "modality_clinical": clinical_df,
    }

    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(db_path)
    for name, df in tables.items():
        df.to_sql(name, conn, index=False)
        ledger.log_step(
            step="consolidate_table",
            source=f"simulate_multiomics_data.py:{name}",
            transform=f"write {len(df)} rows x {len(df.columns)} cols -> sqlite table '{name}'",
            table=name, row_count=len(df), df=df,
        )

    for idx_sql in [
        "CREATE INDEX idx_samples_study ON samples(study_id)",
        "CREATE INDEX idx_scrna_sample ON modality_scrna(sample_id)",
        "CREATE INDEX idx_spatial_sample ON modality_spatial(sample_id)",
        "CREATE INDEX idx_tcrbcr_sample ON modality_tcr_bcr(sample_id)",
        "CREATE INDEX idx_clinical_sample ON modality_clinical(sample_id)",
    ]:
        conn.execute(idx_sql)
    conn.commit()
    conn.close()

    with open(DATA_DIR / "ground_truth_biomarkers.json", "w") as f:
        json.dump(GROUND_TRUTH_BIOMARKERS, f, indent=2)

    ledger.log_step(
        step="inject_known_defects", source="simulate_multiomics_data.py",
        transform="synthetic missingness/out-of-range/duplicate injection for integrity testing",
        table="all", row_count=n_defects, df=None,
    )
    from governance.lineage import LINEAGE_PATH
    ledger.flush(path=Path(ledger_path) if ledger_path else LINEAGE_PATH)

    print(f"Simulated {len(studies_df)} studies, {len(samples_df)} samples.")
    print(f"Modalities: scRNA ({len(SCRNA_GENES)} genes), spatial ({len(SPATIAL_FEATURES)} features), "
          f"TCR/BCR ({len(TCR_BCR_FEATURES)} features), clinical ({len(CLINICAL_FEATURES)} features).")
    print(f"Injected {n_defects} known data-quality defects for integrity validation.")
    print(f"Unified SQLite store written to {db_path}")
    return studies_df, samples_df, ledger


if __name__ == "__main__":
    simulate()
