"""
reproducibility_check.py

A real reproducibility control: re-runs the full ingest pipeline
(simulate_multiomics_data.simulate) twice from the same fixed seed into
two independent, isolated SQLite stores, then verifies that:

  1. Every consolidated table is byte-for-byte content-identical between
     the two runs (via the same canonical SHA-256 hash the lineage ledger
     uses, so "reproducible" means exactly what lineage records claim).
  2. Downstream biomarker discovery run on each store produces identical
     top-18 discovered biomarkers and identical validation metrics.
  3. The data-integrity report computed on each store is identical.

If all three hold, the pipeline is genuinely deterministic/reproducible
given the same input seed -- not merely asserted to be.

Run with: python -m governance.reproducibility_check
"""

import shutil
from pathlib import Path

import simulate_multiomics_data as sim
import discover_biomarkers as disc
from governance import data_integrity
from governance.lineage import _hash_df
import sqlite3
import pandas as pd

TMP_DIR = Path("data") / "_reproducibility_tmp"


def _table_hashes(db_path):
    conn = sqlite3.connect(db_path)
    names = [r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")]
    hashes = {}
    for name in names:
        df = pd.read_sql(f"SELECT * FROM {name}", conn)
        hashes[name] = _hash_df(df)
    conn.close()
    return hashes


def run_check(seed=sim.RNG_SEED):
    if TMP_DIR.exists():
        shutil.rmtree(TMP_DIR)
    TMP_DIR.mkdir(parents=True)

    db_a = TMP_DIR / "run_a.db"
    db_b = TMP_DIR / "run_b.db"
    ledger_a = TMP_DIR / "lineage_a.jsonl"
    ledger_b = TMP_DIR / "lineage_b.jsonl"

    sim.simulate(seed=seed, db_path=db_a, ledger_path=ledger_a, run_id="reproA")
    sim.simulate(seed=seed, db_path=db_b, ledger_path=ledger_b, run_id="reproB")

    hashes_a = _table_hashes(db_a)
    hashes_b = _table_hashes(db_b)
    table_match = hashes_a == hashes_b

    ranked_a, discovered_a, val_a = disc.main(db_path=db_a, out_path=TMP_DIR / "discovered_a.json")
    ranked_b, discovered_b, val_b = disc.main(db_path=db_b, out_path=TMP_DIR / "discovered_b.json")
    biomarkers_match = list(discovered_a["feature"]) == list(discovered_b["feature"])
    validation_match = val_a == val_b

    integrity_a = data_integrity.compute_integrity_report(db_a)
    integrity_b = data_integrity.compute_integrity_report(db_b)
    integrity_match = integrity_a == integrity_b

    passed = table_match and biomarkers_match and validation_match and integrity_match

    shutil.rmtree(TMP_DIR)

    return {
        "reproducible": bool(passed),
        "checks": {
            "consolidated_tables_identical": bool(table_match),
            "discovered_biomarkers_identical": bool(biomarkers_match),
            "validation_metrics_identical": bool(validation_match),
            "integrity_report_identical": bool(integrity_match),
        },
        "n_tables_compared": len(hashes_a),
        "integrity_pct_run_a": integrity_a["overall_integrity_pct"],
        "integrity_pct_run_b": integrity_b["overall_integrity_pct"],
    }


if __name__ == "__main__":
    import json
    print(json.dumps(run_check(), indent=2))
