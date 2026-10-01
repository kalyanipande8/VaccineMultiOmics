"""
test_pipeline.py

Sanity checks run against a real pipeline execution (module-scoped fixture
runs simulate -> discover -> governance once; tests just assert on the
real artifacts produced, nothing is mocked or hardcoded).
"""

import json
import sqlite3
from pathlib import Path

import pytest

import simulate_multiomics_data as sim
import discover_biomarkers as disc
from governance import data_integrity, reproducibility_check


@pytest.fixture(scope="module")
def pipeline_run(tmp_path_factory):
    tmp_dir = tmp_path_factory.mktemp("pipeline")
    db_path = tmp_dir / "vaccine_omics.db"
    ledger_path = tmp_dir / "lineage.jsonl"
    sim.simulate(db_path=db_path, ledger_path=ledger_path, run_id="test_run")
    ranked, discovered, validation = disc.main(db_path=db_path, out_path=tmp_dir / "discovered.json")
    integrity_report = data_integrity.compute_integrity_report(db_path)
    return {
        "db_path": db_path, "ranked": ranked, "discovered": discovered,
        "validation": validation, "integrity_report": integrity_report,
    }


def test_50_plus_studies_present(pipeline_run):
    conn = sqlite3.connect(pipeline_run["db_path"])
    n_studies = conn.execute("SELECT COUNT(DISTINCT study_id) FROM studies").fetchone()[0]
    conn.close()
    assert n_studies >= 50


def test_all_four_modalities_consolidated(pipeline_run):
    conn = sqlite3.connect(pipeline_run["db_path"])
    tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
    conn.close()
    for expected in ["studies", "samples", "modality_scrna", "modality_spatial", "modality_tcr_bcr", "modality_clinical"]:
        assert expected in tables


def test_18_ground_truth_biomarkers_defined():
    assert len(sim.GROUND_TRUTH_BIOMARKERS) == 18


def test_biomarker_discovery_recovers_most_ground_truth(pipeline_run):
    validation = pipeline_run["validation"]
    assert validation["n_discovered"] == 18
    # discovery should meaningfully outperform random guessing (18/~100
    # candidates at random would recover ~3); require at least half back.
    assert validation["n_true_positives"] >= 9
    assert validation["recall"] >= 0.5


def test_data_integrity_in_reasonable_band(pipeline_run):
    pct = pipeline_run["integrity_report"]["overall_integrity_pct"]
    # Defects are deliberately injected, so integrity should be high but
    # not a suspicious, untestable 100%.
    assert 97.0 <= pct <= 100.0


def test_reproducibility_check_genuinely_passes():
    result = reproducibility_check.run_check()
    assert result["reproducible"] is True
    for check_name, passed in result["checks"].items():
        assert passed is True, f"reproducibility sub-check failed: {check_name}"


def test_lineage_ledger_has_required_fields(pipeline_run):
    with open(pipeline_run["db_path"].parent / "lineage.jsonl") as f:
        records = [json.loads(line) for line in f if line.strip()]
    assert len(records) > 0
    for r in records:
        assert {"lineage_version", "source", "transform", "timestamp_utc", "run_id"} <= r.keys()
