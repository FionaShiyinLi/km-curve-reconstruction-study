#!/usr/bin/env python3
"""Recalculate S12 from saved IPD; no model calls or changes to primary scoring.

Run: python scripts/analyze_followup_windows.py (after environment activation).
Writes dataset_results.csv and summary.json under ignored build/followup_windows/.
Late TAAE uses the same [24, U] for all ten arms per dataset, with no extrapolation.
Bootstrap: 10,000 paired draws, 8 datasets per scenario, PCG64 seed 2026100406,
type-1 percentile endpoints. These are the original post hoc analysis settings.
"""
import csv
import json
from pathlib import Path

import numpy as np

from calculate_uncertainty_and_events import arm_arrays, read_records

ROOT = Path(__file__).resolve().parents[1]
CONDITIONS = ("PDF", "NOMARK_PDF", "PNG300", "NOMARK_PNG300")
SEED = 2026100406
DRAWS = 10000


def load_arms(path, n):
    rows = read_records(ROOT / path)
    assert len(rows) == 2 * n and {r[0] for r in rows} == {"A", "B"}, path
    arms = {g: arm_arrays(rows, g) for g in "AB"}
    for values in arms.values():
        assert values is not None and len(values[0]) == n, path
        assert np.all(values[0] <= 36), path
    return arms


def km(time, event):
    knots = np.unique(time[event == 1])
    # Censoring tied to an event remains in the risk set for that event.
    risk = len(time) - np.searchsorted(np.sort(time), knots, side="left")
    deaths = np.array([np.sum((time == t) & (event == 1)) for t in knots])
    return knots, np.cumprod(1 - deaths / risk)


def taae(reference, reconstructed, start, end):
    assert end > start
    kr, sr = km(*reference)
    kx, sx = km(*reconstructed)
    grid = np.unique(np.concatenate(([start, end], kr[(kr > start) & (kr < end)],
                                     kx[(kx > start) & (kx < end)])))
    r = np.r_[1.0, sr][np.searchsorted(kr, grid[:-1], side="right")]
    x = np.r_[1.0, sx][np.searchsorted(kx, grid[:-1], side="right")]
    return float(100 * np.sum(np.diff(grid) * np.abs(x - r)) / (end - start))


def event_counts(time, event):
    return {"early": int(event[(time > 0) & (time <= 24)].sum()),
            "later": int(event[time > 24].sum()), "full": int(event.sum()),
            "at_zero": int(event[time == 0].sum())}


def main():
    with (ROOT / "data/run_manifest.csv").open(newline="") as stream:
        runs = [r for r in csv.DictReader(stream) if r["batch"] == "first"
                and r["complete"] == "True" and r["condition"] in CONDITIONS
                and (r["primary64"] == "True" or r["experiment"] == "followup")]
    lookup = {(r["dataset_id"], r["condition"]): r for r in runs}
    datasets = sorted({r["dataset_id"] for r in runs})
    assert len(runs) == len(lookup) == 128 and len(datasets) == 32
    rows, endpoints = [], []
    for ds in datasets:
        assigned = [lookup[ds, c] for c in CONDITIONS]
        assert len({r["reference_ipd"] for r in assigned}) == 1
        assert len({(r["scenario"], r["n_per_arm"]) for r in assigned}) == 1
        n = int(assigned[0]["n_per_arm"])
        reference = load_arms(assigned[0]["reference_ipd"], n)
        reconstructed = {c: load_arms(lookup[ds, c]["reconstruction_output"], n)
                         for c in CONDITIONS}
        # Final censoring limits support even when a curve has already reached zero.
        end = min(36., *(a[g][0].max() for a in [reference, *reconstructed.values()]
                        for g in "AB"))
        assert end > 24, f"{ds}: no common late-follow-up support"
        endpoints.append(float(end))
        for condition, rec in reconstructed.items():
            row = dict(dataset=ds, scenario=assigned[0]["scenario"], condition=condition,
                       run_id=lookup[ds, condition]["run_id"], endpoint_months=float(end))
            for arm in "AB":
                for label, data in (("reference", reference[arm]), ("reconstructed", rec[arm])):
                    counts = event_counts(*data)
                    assert counts["early"] + counts["later"] + counts["at_zero"] == counts["full"]
                    row.update({f"{label}_{w}_events_{arm}": v for w, v in counts.items()})
                row[f"early_taae_pp_{arm}"] = taae(reference[arm], rec[arm], 0, 24)
                row[f"late_taae_pp_{arm}"] = taae(reference[arm], rec[arm], 24, end)
            for w in ("early", "late"):
                row[f"{w}_taae_pp"] = (row[f"{w}_taae_pp_A"] + row[f"{w}_taae_pp_B"]) / 2
            rows.append(row)

    by_pair = {(r["dataset"], r["condition"]): r for r in rows}
    summary = dict(datasets=32, units="percentage points", bootstrap_draws=DRAWS,
                   seed=SEED, rng="PCG64", percentile_type=1, conditions=[], contrasts=[])
    for c in CONDITIONS:
        subset = [r for r in rows if r["condition"] == c]
        summary["conditions"].append(dict(condition=c,
            mean_early_taae_pp=float(np.mean([r["early_taae_pp"] for r in subset])),
            mean_late_taae_pp=float(np.mean([r["late_taae_pp"] for r in subset])),
            exact_event_counts_both_arms={w: sum(all(r[f"reference_{w}_events_{g}"] ==
                r[f"reconstructed_{w}_events_{g}"] for g in "AB") for r in subset)
                for w in ("early", "later", "full", "at_zero")}))
    rng = np.random.Generator(np.random.PCG64(SEED))
    picks = {f"S{i}": rng.integers(0, 8, size=(DRAWS, 8)) for i in range(1, 5)}
    for window in ("early", "late"):
        for fmt in ("PDF", "PNG300"):
            means, draws = [], []
            for scenario, indices in picks.items():
                ids = [ds for ds in datasets if lookup[ds, fmt]["scenario"] == scenario]
                assert len(ids) == 8
                delta = np.array([by_pair[ds, "NOMARK_" + fmt][window + "_taae_pp"] -
                                  by_pair[ds, fmt][window + "_taae_pp"] for ds in ids])
                means.append(delta.mean())
                draws.append(delta[indices].mean(axis=1))
            lo, hi = np.quantile(np.mean(draws, axis=0), [.025, .975], method="inverted_cdf")
            summary["contrasts"].append(dict(format=fmt, window=window,
                contrast="without marks minus with marks", mean_difference_pp=float(np.mean(means)),
                lower=float(lo), upper=float(hi)))
    summary["support_endpoint_months"] = dict(zip(
        ("minimum", "q1", "median", "q3", "maximum"),
        np.quantile(endpoints, [0, .25, .5, .75, 1]).tolist()))
    target = ROOT / "build/followup_windows"
    target.mkdir(parents=True, exist_ok=True)
    with (target / "dataset_results.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (target / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    print("Calculated S12 for 32 datasets and 128 reconstructions: build/followup_windows/")


if __name__ == "__main__":
    main()
