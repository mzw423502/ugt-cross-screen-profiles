"""Recompute the UGT series statistics from reaction-level inputs.

This module is intentionally independent of frozen summaries and figure tables.
It reads only S3 (identity-included reaction units), S1 (sequence identity),
S2 (acceptor identity, for the input audit) and the local configuration.  The
generated summaries are verification targets for the scientific manuscript;
they are not read while computing the estimates.
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


PACKAGE = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = PACKAGE / "05_REPRODUCIBILITY" / "input"
DEFAULT_OUTPUT = PACKAGE / "05_REPRODUCIBILITY" / "generated"


def _bool(value: object, field: str, row: int) -> bool:
    """Parse a strict boolean field without treating arbitrary text as False."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    text = str(value).strip().upper()
    if text in {"TRUE", "1", "YES", "Y"}:
        return True
    if text in {"FALSE", "0", "NO", "N"}:
        return False
    raise ValueError(f"{field} row {row}: expected an explicit boolean, got {value!r}")


def _status_active(value: object, field: str, row: int) -> int | None:
    """Return 1/0 for tested calls and None for non-comparable states."""
    status = str(value).strip().upper()
    if status == "TESTED_ACTIVE":
        return 1
    if status == "TESTED_INACTIVE":
        return 0
    if status in {"NOT_TESTED", "AMBIGUOUS", "EXCLUDED", "NOT_COMPARABLE"}:
        return None
    raise ValueError(f"{field} row {row}: unknown activity status {value!r}")


def load_inputs(input_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    config = json.loads((input_root / "recompute_config.json").read_text(encoding="utf-8"))
    s3 = pd.read_csv(input_root / config["reaction_units_file"])
    s1 = pd.read_csv(input_root / config["enzyme_identity_file"])
    s2 = pd.read_csv(input_root / config["acceptor_identity_file"])

    required_s3 = {
        "enzyme_id", "key14", "y_status", "m_status", "identity_include",
        "focal_series", "acceptor_label", "comparable_status",
    }
    missing = sorted(required_s3.difference(s3.columns))
    if missing:
        raise ValueError(f"S3 is missing required columns: {missing}")
    required_s1 = {"enzyme_id", "y18_sequence_sha256", "m25_sequence_sha256"}
    missing = sorted(required_s1.difference(s1.columns))
    if missing:
        raise ValueError(f"S1 is missing required columns: {missing}")
    if not {"key14", "m_names"}.issubset(s2.columns):
        raise ValueError("S2 must contain key14 and m_names")

    # The primary key is an input invariant.  A duplicate is an ambiguity, not
    # a row that may be silently dropped.
    duplicate_mask = s3.duplicated(["enzyme_id", "key14"], keep=False)
    if duplicate_mask.any():
        examples = s3.loc[duplicate_mask, ["enzyme_id", "key14"]].head(10).to_dict("records")
        raise ValueError(f"duplicate enzyme–acceptor primary keys in S3: {examples}")
    if s3["enzyme_id"].nunique() != 40 or s3["key14"].nunique() != 15 or len(s3) != 600:
        raise ValueError(
            f"S3 shape invariant failed: rows={len(s3)}, enzymes={s3.enzyme_id.nunique()}, "
            f"acceptors={s3.key14.nunique()}"
        )

    s3["identity_include_bool"] = [
        _bool(v, "identity_include", i) for i, v in enumerate(s3["identity_include"], start=2)
    ]
    if not s3["identity_include_bool"].all():
        raise ValueError("this release expects all 600 S3 rows to be identity-included")
    if s1["enzyme_id"].duplicated().any() or s1.enzyme_id.nunique() != 40:
        raise ValueError("S1 must contain one identity row for each of the 40 enzymes")
    if s2["key14"].duplicated().any() or s2.key14.nunique() < 15:
        raise ValueError("S2 must contain one identity row for every acceptor used by S3")
    if not set(s3["key14"]).issubset(set(s2["key14"])):
        raise ValueError("S3 contains an acceptor key absent from S2")

    # S3 carries an explicit focal-series field.  The six focal keys are also
    # pinned in the config so a fuzzy name match cannot expand the analysis.
    focal = config["focal_acceptors"]
    focal_keys = set(focal)
    if set(s3.loc[s3.focal_series.isin(["HCA", "Coumarins"]), "key14"]) != focal_keys:
        raise ValueError("S3 focal_series does not match the six configured focal keys")
    s2_labels = dict(zip(s2["key14"], s2["m_names"]))
    for key, spec in focal.items():
        if key not in s2_labels:
            raise ValueError(f"configured focal key {key} is absent from S2")
        expected_s2_label = str(spec.get("s2_label", spec["label"])).strip().lower()
        if expected_s2_label != str(s2_labels[key]).strip().lower():
            # S2 keeps the source label, so this check is exact after the
            # key14 lookup; it never fuzzy-matches acceptor names.
            raise ValueError(f"S2 label mismatch for configured focal key {key}")
        observed = set(s3.loc[s3.key14.eq(key), "focal_series"].dropna())
        if observed != {spec["series"]}:
            raise ValueError(f"S3 series mismatch for configured focal key {key}: {observed}")

    # Status fields, rather than the convenience y_active/m_active columns,
    # define the two binary calls.  NOT_TESTED and AMBIGUOUS stay outside the
    # denominator.
    s3["y_call"] = [
        _status_active(v, "y_status", i) for i, v in enumerate(s3["y_status"], start=2)
    ]
    s3["m_call"] = [
        _status_active(v, "m_status", i) for i, v in enumerate(s3["m_status"], start=2)
    ]
    s3["comparable"] = s3["y_call"].notna() & s3["m_call"].notna()
    declared_comparable = s3["comparable_status"].eq("comparable_tested")
    if not (declared_comparable == s3["comparable"]).all():
        bad = s3.loc[declared_comparable != s3["comparable"], ["enzyme_id", "key14", "y_status", "m_status", "comparable_status"]].head(10)
        raise ValueError(f"S3 comparable_status disagrees with explicit statuses:\n{bad}")
    return s3, s1, s2, config


def family_prefix(value: object, pattern: str) -> str:
    match = re.match(pattern, str(value))
    return match.group(1) if match else str(value)


def _pooled(frame: pd.DataFrame, series: str) -> dict:
    subset = frame[frame["focal_series"].eq(series)]
    n = int(len(subset))
    y = int(subset["y_call"].sum())
    m = int(subset["m_call"].sum())
    return {
        "tested_units": n,
        "y18_positive": y,
        "m25_positive": m,
        "pooled_delta": float((m - y) / n) if n else None,
    }


def _per_enzyme(s3: pd.DataFrame, exact_ids: set[str], family_pattern: str) -> pd.DataFrame:
    focal = s3[s3["comparable"] & s3["focal_series"].isin(["HCA", "Coumarins"])].copy()
    rows: list[dict] = []
    for enzyme_id in sorted(s3["enzyme_id"].unique()):
        subset = focal[focal["enzyme_id"].eq(enzyme_id)]
        row = {
            "enzyme_id": enzyme_id,
            "sequence_exact_status": "exact" if enzyme_id in exact_ids else "high_identity_core",
            "family_prefix": family_prefix(enzyme_id, family_pattern),
        }
        for series in ["HCA", "Coumarins"]:
            current = subset[subset["focal_series"].eq(series)]
            n = int(len(current))
            y = int(current["y_call"].sum())
            m = int(current["m_call"].sum())
            row[f"n_{series}"] = n
            row[f"y_{series}"] = y
            row[f"m_{series}"] = m
            row[f"delta_{series}"] = float((m - y) / n) if n else np.nan
        row["D_HCA_minus_Coumarins"] = (
            row["delta_HCA"] - row["delta_Coumarins"]
            if row["n_HCA"] and row["n_Coumarins"]
            else np.nan
        )
        row["included_in_series_analysis"] = bool(row["n_HCA"] >= 2 and row["n_Coumarins"] >= 2)
        rows.append(row)
    return pd.DataFrame(rows)


def _bootstrap(per_enzyme: pd.DataFrame, seed: int, n_resamples: int, scope: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Shared family-prefix draws, preserving repeated cluster multiplicity."""
    clusters = sorted(per_enzyme["family_prefix"].unique())
    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    for replicate in range(1, n_resamples + 1):
        chosen = rng.choice(clusters, size=len(clusters), replace=True)
        sampled = pd.concat(
            [per_enzyme[per_enzyme["family_prefix"].eq(cluster)] for cluster in chosen],
            ignore_index=True,
        )
        hca = float(sampled["delta_HCA"].mean())
        coumarin = float(sampled["delta_Coumarins"].mean())
        d = float(sampled["D_HCA_minus_Coumarins"].mean())
        if not np.isclose(d, hca - coumarin, atol=1e-12, rtol=0):
            raise AssertionError(f"D identity failed in {scope} replicate {replicate}")
        rows.append({"analysis_set": scope, "replicate": replicate, "delta_HCA": hca, "delta_coumarin": coumarin, "D": d})
    draws = pd.DataFrame(rows)
    summary_rows = []
    for metric, column, label in [("HCA", "delta_HCA", "HCA"), ("coumarins", "delta_coumarin", "coumarins"), ("D", "D", "D")]:
        values = draws[column].to_numpy(dtype=float)
        summary_rows.append({
            "analysis_set": scope,
            "metric_id": metric,
            "display_label": "ΔHCA" if label == "HCA" else "Δcoumarin" if label == "coumarins" else "D = ΔHCA − Δcoumarin",
            "estimate": float(per_enzyme[column if column != "delta_coumarin" else "delta_Coumarins"].mean()) if column != "D" else float(per_enzyme["D_HCA_minus_Coumarins"].mean()),
            "lower": float(np.quantile(values, 0.025)),
            "upper": float(np.quantile(values, 0.975)),
            "bootstrap_mean": float(values.mean()),
            "n_enzymes": int(len(per_enzyme)),
            "n_groups": int(len(clusters)),
            "seed": int(seed),
            "n_resamples": int(n_resamples),
        })
    return draws, pd.DataFrame(summary_rows)


def recompute(input_root: Path, output_root: Path) -> dict:
    s3, s1, s2, config = load_inputs(input_root)
    output_root.mkdir(parents=True, exist_ok=True)

    exact_ids = set(s1.loc[s1["y18_sequence_sha256"].eq(s1["m25_sequence_sha256"]), "enzyme_id"])
    family_pattern = config["family_prefix_regex"]
    per = _per_enzyme(s3, exact_ids, family_pattern)
    analysis = per[per["included_in_series_analysis"]].copy()
    exact_analysis = analysis[analysis["sequence_exact_status"].eq("exact")].copy()
    if len(exact_ids) != 31 or len(analysis) != 39 or len(exact_analysis) != 30:
        raise AssertionError(f"sequence/coverage invariants failed: exact={len(exact_ids)}, all={len(analysis)}, exact_eligible={len(exact_analysis)}")
    if analysis["family_prefix"].nunique() != 17 or exact_analysis["family_prefix"].nunique() != 13:
        raise AssertionError("family-prefix group counts do not match the frozen analysis")
    for frame in [analysis, exact_analysis]:
        if not np.isclose(frame["D_HCA_minus_Coumarins"].mean(), frame["delta_HCA"].mean() - frame["delta_Coumarins"].mean(), atol=1e-12, rtol=0):
            raise AssertionError("raw D mean identity failed")

    seed = int(config["bootstrap"]["seed"])
    n_resamples = int(config["bootstrap"]["n_resamples"])
    all_draws, all_summary = _bootstrap(analysis, seed, n_resamples, "all_39_series_adequate")
    exact_draws, exact_summary = _bootstrap(exact_analysis, seed, n_resamples, "exact31_30_series_adequate")
    bootstrap_full = pd.concat([all_draws, exact_draws], ignore_index=True)
    bootstrap_summary = pd.concat([all_summary, exact_summary], ignore_index=True)

    # All 600 rows and all tested states are recomputed here, not read from a
    # summary file.  Pair states are only assigned after both calls are tested.
    tested = s3[s3["comparable"]].copy()
    pair_state = np.select(
        [tested.y_call.eq(1) & tested.m_call.eq(1), tested.y_call.eq(0) & tested.m_call.eq(0), tested.y_call.eq(1) & tested.m_call.eq(0), tested.y_call.eq(0) & tested.m_call.eq(1)],
        ["double_positive", "double_not_detected", "y18_only", "m25_only"],
        default="invalid",
    )
    if (pair_state == "invalid").any():
        raise AssertionError("invalid pair state encountered")
    pair_counts = pd.Series(pair_state).value_counts().to_dict()
    audit = {
        "identity_included_units": int(len(s3)),
        "enzyme_n": int(s3.enzyme_id.nunique()),
        "acceptor_n": int(s3.key14.nunique()),
        "y18_not_tested_units": int(s3.y_call.isna().sum()),
        "tested_pairs": int(len(tested)),
        "four_state_counts": {str(k): int(v) for k, v in pair_counts.items()},
        "exact_sequence_enzyme_n": int(len(exact_ids)),
        "exact_sequence_comparable_units": int(len(s3[s3.enzyme_id.isin(exact_ids) & s3.comparable])),
        "series_eligible_enzyme_n": int(len(analysis)),
        "exact_series_eligible_enzyme_n": int(len(exact_analysis)),
    }

    molecule_rows = []
    focal = s3[s3["comparable"] & s3["focal_series"].isin(["HCA", "Coumarins"])].copy()
    for key, spec in config["focal_acceptors"].items():
        x = focal[focal.key14.eq(key)]
        molecule_rows.append({
            "key14": key, "acceptor_label": spec["label"], "series": spec["series"],
            "tested_units": int(len(x)), "y18_positive": int(x.y_call.sum()), "m25_positive": int(x.m_call.sum()),
            "y18_fraction": float(x.y_call.mean()) if len(x) else None, "m25_fraction": float(x.m_call.mean()) if len(x) else None,
        })
    molecule_counts = pd.DataFrame(molecule_rows)

    exact_comp = s3[s3.enzyme_id.isin(exact_ids) & s3.comparable].copy()
    exact_series_rows = []
    for series in ["HCA", "Coumarins"]:
        exact_series_rows.append({"series": series, **_pooled(exact_comp, series)})
    exact_series_summary = pd.DataFrame(exact_series_rows)
    exact_series_summary["focal_series"] = exact_series_summary["series"]
    exact_series_summary["enzyme_n"] = exact_series_summary["series"].map({s: int(exact_comp[exact_comp.focal_series.eq(s)].enzyme_id.nunique()) for s in ["HCA", "Coumarins"]})
    exact_series_summary["acceptor_n"] = exact_series_summary["series"].map({s: int(exact_comp[exact_comp.focal_series.eq(s)].key14.nunique()) for s in ["HCA", "Coumarins"]})
    exact_series_summary["y18_active_fraction"] = exact_series_summary["y18_positive"] / exact_series_summary["tested_units"]
    exact_series_summary["m25_active_fraction"] = exact_series_summary["m25_positive"] / exact_series_summary["tested_units"]
    exact_series_summary["delta_m25_minus_y18"] = exact_series_summary["m25_active_fraction"] - exact_series_summary["y18_active_fraction"]

    per.to_csv(output_root / "per_enzyme_series.csv", index=False, encoding="utf-8-sig")
    molecule_counts.to_csv(output_root / "molecule_counts.csv", index=False, encoding="utf-8-sig")
    exact_series_summary.to_csv(output_root / "exact31_series_summary.csv", index=False, encoding="utf-8-sig")
    bootstrap_full.to_csv(output_root / "bootstrap_full.csv", index=False, encoding="utf-8-sig")
    bootstrap_summary.to_csv(output_root / "bootstrap_summary.csv", index=False, encoding="utf-8-sig")
    compatibility = bootstrap_summary.rename(columns={"analysis_set": "scope", "metric_id": "metric", "estimate": "raw_point_estimate", "lower": "bootstrap_low", "upper": "bootstrap_high", "n_groups": "n_clusters"}).copy()
    compatibility["weighting"] = "equal enzyme weight; UGT family-prefix cluster bootstrap"
    compatibility["bootstrap_mean"] = bootstrap_summary["bootstrap_mean"]
    compatibility["seed"] = seed
    compatibility["n_resamples"] = n_resamples
    compatibility[["scope", "metric", "raw_point_estimate", "bootstrap_low", "bootstrap_high", "bootstrap_mean", "n_enzymes", "n_clusters", "weighting", "seed", "n_resamples"]].to_csv(output_root / "bootstrap_series_final.csv", index=False, encoding="utf-8-sig")
    # The audit has a derived pair state only where both source calls were
    # actually tested; non-comparable rows remain explicit.
    audit_frame = s3[["enzyme_id", "key14", "y_status", "m_status", "comparable"]].copy()
    audit_frame["derived_pair_status"] = "not_comparable"
    audit_frame.loc[audit_frame.comparable, "derived_pair_status"] = pair_state
    audit_frame.to_csv(output_root / "reaction_audit.csv", index=False, encoding="utf-8-sig")
    enzyme_paired = s3[s3["comparable"]].groupby("enzyme_id", as_index=False).agg(n_comparable=("key14", "size"))
    enzyme_paired.to_csv(output_root / "enzyme_paired_summary.csv", index=False, encoding="utf-8-sig")
    per.to_csv(output_root / "per_enzyme_series_final.csv", index=False, encoding="utf-8-sig")
    exact_series_summary.to_csv(output_root / "exact31_series_summary.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame({"enzyme_id": sorted(s3.enzyme_id.unique())}).assign(
        family_prefix=[family_prefix(e, family_pattern) for e in sorted(s3.enzyme_id.unique())],
        sequence_exact_status=["exact" if e in exact_ids else "high_identity_core" for e in sorted(s3.enzyme_id.unique())],
    ).to_csv(output_root / "family_prefix_map.csv", index=False, encoding="utf-8-sig")
    bootstrap_summary[["analysis_set", "metric_id", "display_label", "estimate", "lower", "upper", "n_enzymes", "n_groups"]].assign(
        display_order=range(1, len(bootstrap_summary) + 1)
    ).to_csv(output_root / "Fig2_effects.csv", index=False, encoding="utf-8-sig")

    final_numbers = {
        "audit": audit,
        "focal_molecule_counts": molecule_counts.to_dict("records"),
        "all_39_series_adequate": {"n_enzymes": len(analysis), "n_groups": int(analysis.family_prefix.nunique()), "pooled_descriptive": {s: _pooled(focal[focal.enzyme_id.isin(analysis.enzyme_id)], s) for s in ["HCA", "Coumarins"]}, "summary": all_summary.to_dict("records")},
        "exact31_30_series_adequate": {"exact_sequence_enzyme_n": len(exact_ids), "n_enzymes": len(exact_analysis), "n_groups": int(exact_analysis.family_prefix.nunique()), "pooled_descriptive": {s: _pooled(exact_comp, s) for s in ["HCA", "Coumarins"]}, "summary": exact_summary.to_dict("records")},
        "bootstrap": {"seed": seed, "n_resamples": n_resamples, "rng": "numpy.random.default_rng", "cluster_order": sorted(analysis.family_prefix.unique()), "reset_seed_per_analysis_set": True, "quantile_method": "numpy.quantile default linear", "full_output": "bootstrap_full.csv", "D_recomputed_each_replicate": True},
    }
    (output_root / "FINAL_NUMBERS.json").write_text(json.dumps(final_numbers, ensure_ascii=False, indent=2), encoding="utf-8")
    (output_root / "recompute_config_used.json").write_text(json.dumps(config, ensure_ascii=False, indent=2), encoding="utf-8")
    return {"audit": audit, "summary_rows": int(len(bootstrap_summary)), "bootstrap_rows": int(len(bootstrap_full)), "output": str(output_root)}


def main(argv: Iterable[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(list(argv) if argv is not None else None)
    result = recompute(args.input_root.resolve(), args.output_root.resolve())
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
