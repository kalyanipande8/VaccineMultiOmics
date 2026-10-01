"""
run_all.py

Runs the full pipeline end-to-end:
  1. Simulate scRNA-seq, spatial, TCR/BCR, and clinical data across 50+
     vaccine studies, consolidated into one unified SQLite store, with
     every step logged to the governance lineage ledger -> data/
  2. Discover ~18 predictive biomarkers and validate against injected
     ground truth -> data/discovered_biomarkers.json
  3. Compute the data-integrity report, run the reproducibility check,
     and produce the regulatory-alignment compliance report -> output/

Run with:  python run_all.py
"""

import json
from pathlib import Path

import simulate_multiomics_data
import discover_biomarkers
from governance import regulatory_checklist


def main():
    print("--- Step 1: Simulating & consolidating multi-omics vaccine data ---")
    simulate_multiomics_data.simulate()

    print("\n--- Step 2: Biomarker discovery & validation ---")
    discover_biomarkers.main()

    print("\n--- Step 3: Data governance (lineage review, reproducibility, compliance) ---")
    report = regulatory_checklist.main()

    print("\n=== Pipeline complete ===")
    print(json.dumps({
        "overall_data_integrity_pct": report["overall_data_integrity_pct"],
        "reproducible": report["reproducibility"]["reproducible"] if report["reproducibility"] else None,
        "compliance_rate_pct": report["compliance_rate_pct"],
    }, indent=2))


if __name__ == "__main__":
    main()
