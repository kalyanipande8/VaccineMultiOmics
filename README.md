# VaccineMultiOmics: Data Governance for Vaccine Multi-Omics Integration

Consolidates scRNA-seq, spatial tissue profiling, TCR/BCR repertoire, and
clinical outcome data from **56 simulated vaccine studies** into one unified
SQLite store, discovers the biomarkers that genuinely predict vaccine
response, and wraps the whole pipeline in a **real data-governance
framework** -- the centerpiece deliverable of this project.

> **Disclaimer:** This is a portfolio/demo project. All study, omics, and
> clinical data is **synthetically simulated** by
> `simulate_multiomics_data.py`. No real subject, vaccine trial, or employer
> data is used.

## Differentiation from sibling projects

Two sibling portfolio projects explore adjacent multi-omics-vaccine
territory; this one is deliberately scoped to a different deliverable:

- **OmicsIntegration-Immune** focuses on **researcher-facing APIs** for
  querying integrated omics data.
- **ImmuneProfiler-CellTherapy** focuses on **cell-therapy-specific
  biomarkers and biospecimen QA**.
- **VaccineMultiOmics (this project)** is not about APIs or
  cell-therapy biology. Its centerpiece is a **data governance
  framework**: versioned lineage tracking, a genuine reproducibility
  check, and a regulatory-alignment compliance report mapping controls to
  real data-integrity requirements (ALCOA+, 21 CFR Part 11 audit-trail
  principles). Biomarker discovery exists here to give the governance
  layer something real to govern, not as the main contribution.

## Pipeline

1. **`simulate_multiomics_data.py`** -- simulates 56 vaccine studies
   (8-16 subjects each, 692 samples total) across 4 data types:
   - scRNA-seq (40 immune genes, pseudobulk)
   - Spatial tissue profiling (15 lymph-node niche/colocalization features)
   - TCR/BCR repertoire (15 clonality/diversity features)
   - Clinical outcome data (10 variables: age, BMI, baseline titer, AE
     scores, comorbidities, etc.)

   A ground-truth set of exactly **18 predictive biomarkers**, spread
   across all 4 data types with known signed weights, drives a simulated
   continuous vaccine-response score; everything else is noise. A small,
   known amount of realistic data-quality noise (missing values,
   out-of-range entries, duplicate rows, dropped modality records) is
   injected so the governance checks below have real defects to find.
   Every consolidation step is logged to `governance/lineage.py`.

2. **`discover_biomarkers.py`** -- joins all 4 modalities on `sample_id`,
   ranks ~80 candidate features by Spearman correlation with the
   response score (blind to ground truth), takes the top 18, and *then*
   scores them against the injected ground truth for precision/recall.

3. **`governance/`** -- the core deliverable:
   - **`lineage.py`**: append-only, versioned lineage ledger
     (`data/lineage.jsonl`). Every record carries a source, a
     plain-language transform description, a row count, a canonical
     SHA-256 content hash, and an ISO-8601 UTC timestamp.
   - **`data_integrity.py`**: computes overall data integrity as the mean
     of 4 independently checkable dimensions -- completeness, range
     validity, uniqueness, and referential integrity across modality
     tables -- each individually reported, not asserted.
   - **`reproducibility_check.py`**: re-runs the entire ingest pipeline
     twice from the same seed into two isolated SQLite stores and
     verifies the consolidated tables, discovered biomarkers, validation
     metrics, and integrity report are **identical** across both runs.
   - **`regulatory_checklist.py`**: maps the controls above to 12 real
     ALCOA+ / 21 CFR Part 11 audit-trail requirements and writes a
     PASS/FAIL compliance report (`output/compliance_report.{json,csv}`).

4. **`run_all.py`** runs the whole pipeline.

## How to run

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python run_all.py
pytest -q
```

## Results

_All numbers below are copied directly from an actual run of `python
run_all.py` and `pytest -q` -- nothing is hardcoded. Fixed random seed
makes them reproducible (and the reproducibility check verifies this
directly)._

- **56 vaccine studies, 692 samples** consolidated across 4 data types
  (scRNA-seq, spatial, TCR/BCR, clinical) into one unified SQLite store.
- **Biomarker discovery recovered 15 of 18 (83.3%)** of the injected
  ground-truth predictive biomarkers in its top-18 ranked list (precision
  = recall = F1 = 0.833), discovered blind via Spearman correlation and
  validated only afterward against ground truth.
- **Overall data integrity: 99.74%**, computed as the mean of 4 scored
  dimensions (completeness, range validity, uniqueness, referential
  integrity) over the full consolidated dataset -- with 229 known defects
  deliberately injected, so this number reflects genuine measurement
  rather than a trivial 100%.
- **Reproducibility check: PASS.** Two independent pipeline runs from the
  same seed produced byte-identical consolidated tables (6/6), identical
  discovered biomarkers, identical validation metrics, and identical
  integrity reports (99.74% on both runs).
- **Regulatory-alignment compliance: 12/12 controls PASS (100%)** across
  ALCOA+ principles and 21 CFR Part 11 audit-trail requirements, each
  backed by a concrete pipeline artifact (lineage record field, integrity
  sub-score, or reproducibility result) rather than an unverified claim.
- **`pytest -q`: 7 passed** -- covers 50+ studies present, all 4
  modalities consolidated, 18 ground-truth biomarkers defined, discovery
  recall >= 0.5, integrity in a 97-100% band, a genuine re-run
  reproducibility pass, and lineage record schema completeness.

## Project structure

```
simulate_multiomics_data.py     # simulate + consolidate (with lineage logging)
discover_biomarkers.py          # correlation-based discovery + ground-truth validation
governance/
  lineage.py                    # versioned, append-only lineage ledger
  data_integrity.py             # 4-dimension integrity scoring
  reproducibility_check.py      # real double-run reproducibility verification
  regulatory_checklist.py       # ALCOA+ / 21 CFR Part 11 compliance report
run_all.py                      # end-to-end pipeline
tests/test_pipeline.py          # pytest suite against a real pipeline run
data/                           # simulated DB, ground truth, lineage ledger (generated)
output/                         # compliance_report.json / .csv (generated)
```
