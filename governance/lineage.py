"""
lineage.py

Versioned data-lineage ledger: every consolidation step performed anywhere
in the pipeline (simulation ingest, biomarker discovery, integrity checks)
is appended as one record to data/lineage.jsonl with:
  - a monotonically increasing lineage version number
  - source (which script/table produced the step)
  - transform (what was done, in plain language)
  - row_count affected
  - a SHA-256 content hash of the resulting table/data (for tamper-evidence
    and for the reproducibility check to compare against)
  - an ISO-8601 UTC timestamp

This is intentionally append-only (never overwritten) so the full history
of how the consolidated dataset came to be is always reconstructable --
the core requirement of a data-governance / lineage-tracking control.
"""

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

LINEAGE_PATH = Path("data") / "lineage.jsonl"


def _hash_df(df):
    if df is None:
        return None
    # Hash on a canonical, order-independent representation so the hash
    # reflects *content*, not incidental row/column ordering.
    canon = df.fillna("__NA__").astype(str)
    canon = canon.reindex(sorted(canon.columns), axis=1)
    canon = canon.sort_values(by=list(canon.columns)).reset_index(drop=True)
    payload = canon.to_csv(index=False).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


class LineageLedger:
    """Buffers lineage records in memory, flushes them to the append-only
    JSONL ledger. A fresh ledger file is started each pipeline run (the
    run's full lineage is self-contained), but every record's `run_id`
    distinguishes separate runs if the file is kept across runs."""

    def __init__(self, run_id=None):
        self._records = []
        self.run_id = run_id or datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        self._next_version = 1

    def log_step(self, step, source, transform, table, row_count, df=None):
        record = {
            "lineage_version": self._next_version,
            "run_id": self.run_id,
            "step": step,
            "source": source,
            "transform": transform,
            "table": table,
            "row_count": int(row_count),
            "content_hash": _hash_df(df),
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        }
        self._records.append(record)
        self._next_version += 1
        return record

    def flush(self, path=LINEAGE_PATH):
        """Append this run's records to the ledger file (new file per
        pipeline run so `run_all.py` always produces a clean, complete
        lineage trail for its own run; history across runs accumulates
        if the file isn't deleted between runs)."""
        path.parent.mkdir(exist_ok=True)
        with open(path, "a") as f:
            for rec in self._records:
                f.write(json.dumps(rec) + "\n")
        return len(self._records)

    @property
    def records(self):
        return list(self._records)


def read_ledger(path=LINEAGE_PATH):
    if not path.exists():
        return []
    with open(path) as f:
        return [json.loads(line) for line in f if line.strip()]
