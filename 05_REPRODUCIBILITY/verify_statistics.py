"""Verify recomputed outputs against the frozen release targets.

This validator is deliberately separate from ``recompute_statistics.py``.  It
may read the frozen expected values and generated summaries, but those files
are never inputs to the statistical calculation itself.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PACKAGE / "05_REPRODUCIBILITY" / "generated"
EXPECTED = PACKAGE / "05_REPRODUCIBILITY" / "frozen_expected.json"


def verify(output_root: Path, expected_path: Path = EXPECTED, tolerance: float = 1e-10) -> dict:
    expected = json.loads(expected_path.read_text(encoding="utf-8"))
    numbers = json.loads((output_root / "FINAL_NUMBERS.json").read_text(encoding="utf-8"))
    effects = pd.read_csv(output_root / "Fig2_effects.csv")
    draws = pd.read_csv(output_root / "bootstrap_full.csv")
    summary = pd.read_csv(output_root / "bootstrap_summary.csv")

    audit = numbers["audit"]
    for key, value in expected["audit"].items():
        if audit.get(key) != value:
            raise AssertionError(f"audit mismatch for {key}: {audit.get(key)!r} != {value!r}")
    if len(draws) != 1000 or draws.groupby("analysis_set").size().to_dict() != {"all_39_series_adequate": 500, "exact31_30_series_adequate": 500}:
        raise AssertionError("bootstrap_full.csv does not contain two complete 500-row blocks")
    if not np.allclose(draws["D"], draws["delta_HCA"] - draws["delta_coumarin"], atol=tolerance, rtol=0):
        raise AssertionError("D identity failed in the saved bootstrap draws")
    if len(summary) != 6 or len(effects) != 6:
        raise AssertionError("six-row summary/effect outputs are required")
    for row in expected["effects"]:
        match = effects[(effects.analysis_set == row["analysis_set"]) & (effects.metric_id == row["metric_id"])]
        if len(match) != 1:
            raise AssertionError(f"missing effect row {row['analysis_set']} / {row['metric_id']}")
        got = match.iloc[0]
        for field in ["estimate", "lower", "upper"]:
            if not np.isclose(float(got[field]), float(row[field]), atol=tolerance, rtol=0):
                raise AssertionError(f"{field} mismatch for {row['analysis_set']} / {row['metric_id']}: {got[field]} != {row[field]}")
    return {"status": "passed", "tolerance": tolerance, "audit": audit, "bootstrap_rows": len(draws), "effect_rows": len(effects)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected", type=Path, default=EXPECTED)
    parser.add_argument("--tolerance", type=float, default=1e-10)
    args = parser.parse_args()
    result = verify(args.output_root.resolve(), args.expected.resolve(), args.tolerance)
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
