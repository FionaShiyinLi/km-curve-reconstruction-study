#!/usr/bin/env python3
"""Reproduce Figure 3: original image and signed reconstruction errors over time.
Inputs: build/scores.csv, data/run_manifest.csv and saved reference/LLM IPDs.
Outputs: build/figures/Figure3_reconstruction_example.{pdf,csv,json}.
Run: python scripts/plot_reconstruction_example.py (repository root).
Retains the original lower-median selection (first-batch PNG300 rank 16/32),
right-continuous curves, exact union of knots and pandas round-trip parsing.
"""
from pathlib import Path
import hashlib
import json
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
matplotlib.use('Agg')
ROOT=Path(__file__).resolve().parents[1]
FIGURES=ROOT/'build/figures'

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def km(records, group):
    arm = records[records.group == group]
    counts = arm.groupby("observed_time").event.agg(["count", "sum"]).sort_index()
    at_risk = len(arm) - counts["count"].cumsum().shift(fill_value=0)
    survival = (1 - counts["sum"] / at_risk).cumprod()
    return np.r_[0.0, counts.index], np.r_[1.0, survival]

def prepare():
    FIGURES.mkdir(parents=True,exist_ok=True)
    scores=pd.read_csv(ROOT/'build/scores.csv',float_precision='round_trip')
    ranked=scores[scores.primary64 & (scores.batch=='first') & (scores.condition=='PNG300')].sort_values(['taae_mean','dataset_id'])
    assert len(ranked)==32
    case=dict(rank=16,denominator=len(ranked),png_taae=float(ranked.iloc[15].taae_mean))
    dataset=ranked.iloc[15].dataset_id
    assert dataset=='PP32_S3_D03'
    subset=scores[scores.primary64 & (scores.dataset_id==dataset)]
    input_png=ROOT/subset[subset.condition=='PNG300'].iloc[0].input_figure
    reference_path=ROOT/subset.iloc[0].reference_ipd
    records={'Reference':pd.read_csv(reference_path,float_precision='round_trip')}
    sources={str(reference_path.relative_to(ROOT)):sha(reference_path),str(input_png.relative_to(ROOT)):sha(input_png)}
    sessions={}
    for fmt in ('PDF','PNG300'):
        row=subset[subset.condition==fmt].iloc[0]
        path=ROOT/row.reconstruction_output
        records[fmt]=pd.read_csv(path,float_precision='round_trip');sessions[fmt]=row.run_id
        sources[str(path.relative_to(ROOT))]=sha(path)
    tables, metrics, time_shifts = [], {}, {}
    for group in "AB":
        original = records["Reference"][records["Reference"].group == group].sort_values("observed_time")
        native = records["PDF"][records["PDF"].group == group].sort_values("observed_time")
        np.testing.assert_array_equal(original.event, native.event)
        time_shifts[group] = float(np.max(np.abs(original.observed_time.to_numpy() - native.observed_time.to_numpy())))
        assert time_shifts[group] <= 5.21e-6
        curves = {fmt: km(data, group) for fmt, data in records.items()}
        # The union of actual step times preserves all interval widths; no grid or smoothing.
        knots = np.unique(np.concatenate([np.array([0.0, 24.0])] +
                          [t[(t > 0) & (t < 24)] for t, _ in curves.values()]))
        values = {fmt: s[np.searchsorted(t, knots, side="right") - 1]
                  for fmt, (t, s) in curves.items()}
        table = pd.DataFrame({"group": group, "time_months": knots,
                              "interval_end_months": np.r_[knots[1:], 24.0],
                              "reference_survival": values["Reference"]})
        for fmt, column in (("PDF", "pdf"), ("PNG300", "png300")):
            delta = 100 * (values[fmt] - values["Reference"])
            table[f"{column}_survival"] = values[fmt]
            table[f"{column}_error_pp"] = delta
            taae = float(np.dot(np.diff(knots), np.abs(delta[:-1])) / 24)
            maximum = float(np.max(np.abs(delta)))
            scored = scores[scores.run_id == sessions[fmt]].iloc[0]
            np.testing.assert_allclose(taae, 100 * scored[f"taae_{group}"], rtol=0, atol=1e-10)
            np.testing.assert_allclose(maximum, 100 * scored[f"maximum_{group}"], rtol=0, atol=1e-10)
            metrics[f"{fmt}_{group}"] = {"taae_pp": taae, "maximum_absolute_error_pp": maximum}
        tables.append(table)
    data = pd.concat(tables, ignore_index=True)
    means = {fmt: float(np.mean([metrics[f"{fmt}_{group}"]["taae_pp"] for group in "AB"]))
             for fmt in ("PDF", "PNG300")}
    np.testing.assert_allclose(means["PNG300"], 100 * case["png_taae"], rtol=0, atol=1e-10)
    assert data[["pdf_error_pp", "png300_error_pp"]].abs().to_numpy().max() < 1.5
    pdf_mantissa, pdf_exponent = f"{means['PDF']:.2e}".split("e")

    colors = {"PDF": "#0072B2", "PNG300": "#D55E00"}
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 9,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "axes.linewidth": 0.6, "path.simplify": False, "pdf.fonttype": 42}):
        fig = plt.figure(figsize=(6.5, 6.9))
        fig.text(0.035, 0.98, "A  Original Kaplan-Meier figure", fontsize=12, weight="bold", va="top")
        fig.text(0.035, 0.948, "Dataset S3 D03  |  1,000 participants per group", fontsize=9, va="top")
        ax_input = fig.add_axes([0.08, 0.475, 0.87, 0.455])
        ax_input.imshow(plt.imread(input_png), interpolation="none")
        ax_input.set_axis_off()

        fig.text(0.035, 0.444, "B  Reconstruction errors over time", fontsize=12, weight="bold", va="top")
        handles = [Line2D([0], [0], color=colors[fmt], lw=1.5,
                          label="Native PDF" if fmt == "PDF" else "PNG300") for fmt in colors]
        fig.legend(handles=handles, loc="center left", bbox_to_anchor=(0.145, 0.402),
                   frameon=False, ncol=2, fontsize=9, borderaxespad=0)
        fig.text(0.95, 0.402, "Reconstructed minus reference", ha="right", va="center", fontsize=8)
        error_label = (rf"Time-averaged absolute error: PDF ${pdf_mantissa}\times10^{{{int(pdf_exponent)}}}$; "
                       f"PNG300 {means['PNG300']:.3g} pp")
        fig.text(0.15, 0.369, error_label, fontsize=8, va="center")
        for group, bottom in (("A", 0.220), ("B", 0.055)):
            ax = fig.add_axes([0.15, bottom, 0.80, 0.125])
            arm = data[data.group == group]
            ax.axhline(0, color="#9a9a9a", lw=0.6, zorder=1)
            for fmt, column, width in (("PNG300", "png300", 0.75), ("PDF", "pdf", 0.55)):
                ax.step(arm.time_months, arm[f"{column}_error_pp"], where="post",
                        color=colors[fmt], lw=width, zorder=2)
            ax.set(xlim=(0, 24), ylim=(-1.5, 1.5), yticks=[-1, 0, 1], xticks=[0, 6, 12, 18, 24])
            ax.tick_params(labelsize=8, length=3, width=0.6)
            ax.text(0.985, 0.90, f"Group {group}", transform=ax.transAxes, fontsize=9, va="top", ha="right")
            if group == "A":
                ax.tick_params(labelbottom=False)
            else:
                ax.set_xlabel("Time (months)", fontsize=9, labelpad=3)
        fig.text(0.035, 0.20, "Survival probability difference\n(percentage points)",
                 rotation=90, ha="center", va="center", fontsize=9)
        output = FIGURES / "Figure3_reconstruction_example.pdf"
        fig.savefig(output, metadata={"CreationDate": None, "ModDate": None})
        plt.close(fig)

    data.to_csv(FIGURES / "Figure3_reconstruction_example.csv", index=False, float_format="%.17g")
    audit = {"dataset": dataset, "selection": "Existing post-result lower-median PNG300 case",
             "png300_error_rank": case["rank"], "selection_denominator": case["denominator"],
             "session_ids": sessions, "evaluation_months": [0, 24],
             "difference": "100 * (reconstructed survival - reference survival)",
             "display": "All right-continuous KM steps; no smoothing, downsampling, or clipping; identical axes for both groups",
             "input_display": "Unaltered original PNG300 input, with its original 0-36-month axis",
             "metrics": metrics, "mean_taae_pp": means,
             "native_pdf_event_censor_order_preserved": True,
             "native_pdf_maximum_observed_time_shift_months": time_shifts,
             "frozen_score_and_reference_checks": "passed", "sources_sha256": sources,
             "outputs_sha256": {name: sha(FIGURES / name) for name in
                                ("Figure3_reconstruction_example.pdf", "Figure3_reconstruction_example.csv")}}
    (FIGURES / "Figure3_reconstruction_example.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps({key: audit[key] for key in
                     ("dataset", "png300_error_rank", "selection_denominator", "mean_taae_pp", "frozen_score_and_reference_checks")}, indent=2))

if __name__=="__main__":prepare()
