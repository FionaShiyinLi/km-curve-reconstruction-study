#!/usr/bin/env python3
"""Plot paired accuracy comparisons: manuscript Figures 2, 4, S1 and S2.
Inputs: build/{primary_pairs,scores,reconstruction_audit,followup_contrasts,
image_only_pairs,image_only_bootstrap,survminer_pairs}.csv.
Outputs: numbered PDFs and supporting plot data in build/figures/.
Run: python scripts/plot_accuracy_comparisons.py (repository root).
Original figure calculations and layout are retained; S2 has fewer x-axis ticks
to prevent label collisions. Data loading/paths are refactored.
Probability and percentage-point axes remain distinct, as in the manuscript.
"""
from pathlib import Path
import csv
import json
import math
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator
from matplotlib.transforms import blended_transform_factory
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
FIGURES=ROOT/'build/figures'
COLORS={'PDF':'#0072B2','RASTER_PDF':'#009E73','PNG300':'#D55E00'}

def read_csv(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))

def savefig(fig,name):
    fig.savefig(FIGURES/(name+'.pdf'),bbox_inches='tight',metadata={'CreationDate':None,'ModDate':None})
    plt.close(fig)


def primary_figure():
    scores=read_csv('build/scores.csv')
    audit={r['run_id']:r for r in read_csv('build/reconstruction_audit.csv')}
    png={r['dataset_id']:r for r in scores if r['primary64']=='TRUE' and r['condition']=='PNG300'}
    rows=[]
    for pair in read_csv('build/primary_pairs.csv'):
        run=png[pair['dataset']]
        rows.append(dict(dataset=pair['dataset'],scenario=pair['scenario'],native_pdf_taae_pp=100*float(pair['taae_mean_PDF']),
                         png300_taae_pp=100*float(pair['taae_mean_PNG300']),png300_excess_revisions=audit[run['run_id']]['confirmed_excess_corrections']=='True',png300_session_id=run['run_id']))
    rows.sort(key=lambda r:(r['scenario'],r['dataset']))
    assert len(rows)==64 and sum(r['png300_excess_revisions'] for r in rows)==13
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.titleweight": "bold", "legend.frameon": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(9, 4.5))
    fig.subplots_adjust(left=.085, right=.985, bottom=.22, top=.90, wspace=.32)
    for ax in axes:
        for i, row in enumerate(rows):
            offset = ((i * 17) % len(rows) / (len(rows) - 1) - .5) * .13
            values = [row["native_pdf_taae_pp"], row["png300_taae_pp"]]
            ax.plot([offset, 1 + offset], values, color="#b5b5b5", lw=.55, alpha=.65, zorder=1)
            ax.scatter(offset, values[0], color="#0072B2", s=13, linewidths=0, zorder=2)
            marker = "x" if row["png300_excess_revisions"] else "o"
            ax.scatter(1 + offset, values[1], color="#D55E00", marker=marker,
                       s=24 if marker == "x" else 13, linewidths=.9 if marker == "x" else 0, zorder=3)
        ax.set_xticks([0, 1], ["Native PDF\n64 datasets", "PNG300\n64 datasets"])
        ax.set_xlim(-.3, 1.3)
        ax.set_ylabel("Mean TAAE across arms (percentage points)")
        ax.grid(axis="y", alpha=.18)
    axes[0].set(title="A  Absolute magnitude", ylim=(-.012, .36))
    axes[1].set(title="B  Detail on a logarithmic scale", yscale="log", ylim=(3e-6, .65))
    fig.legend(handles=[
        Line2D([], [], marker="x", color="#D55E00", linestyle="none", label="PNG300: >2 confirmed revisions"),
    ], loc="lower center", bbox_to_anchor=(.5, .025), fontsize=9)
    output = FIGURES / "Figure2_primary_accuracy.pdf"
    fig.savefig(output, bbox_inches="tight", metadata={"CreationDate": None, "ModDate": None})
    plt.close(fig)
    with output.with_suffix(".csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def display_number(value):
    if value != 0 and abs(value) < 0.001:
        mantissa, exponent = f"{value:.1e}".split("e")
        return rf"${mantissa}\times10^{{{int(exponent)}}}$"
    return f"{value:.4f}"

def draw_followup_figure(rows):
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8.5,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "mathtext.default": "regular"})
    fig = plt.figure(figsize=(6.5, 6.3), facecolor="white")
    color = "#225C78"
    fig.text(0.04, 0.975, "32 paired datasets per comparison; points = means, lines = 95% CIs",
             fontsize=8.5, color="#424B50")
    fig.text(0.04, 0.92, "Condition\nminus comparator", fontsize=8.5, weight="bold", va="center")
    for suffix, left, title, subtitle, limits, ticks, tick_labels, xlabel in (
        ("pp", 0.375, "A  TAAE", "Primary outcome", (-0.012, 0.165),
         [0, 0.05, 0.10, 0.15], ["0", "0.05", "0.10", "0.15"],
         "Difference in TAAE\n(percentage points)"),
        ("log_hr", 0.71, "B  Log-HR error", "Secondary outcome", (-0.002, 0.0045),
         [-0.002, 0, 0.002, 0.004], ["-0.002", "0", "0.002", "0.004"],
         "Difference in absolute\nlog-HR error"),
    ):
        fig.text(left, 0.93, title, fontsize=10, weight="bold")
        fig.text(left, 0.897, subtitle, fontsize=8, color="#424B50")
        ax = fig.add_axes((left, 0.595, 0.255, 0.28))
        ax.set_xlim(limits)
        ax.set_ylim(-0.6, 4.6)
        ax.set_xticks(ticks, tick_labels)
        ax.set_yticks([])
        ax.grid(axis="x", color="#DFE4E7", linewidth=0.65)
        ax.axvline(0, color="#56616A", linestyle=(0, (3, 3)), linewidth=0.9)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_color("#869199")
        ax.tick_params(axis="x", length=3, color="#869199", labelsize=7.5)
        ax.set_xlabel(xlabel, fontsize=8, labelpad=5)
        label_transform = blended_transform_factory(fig.transFigure, ax.transData)
        for index, row in enumerate(rows):
            y = 4 - index
            estimate, lower, upper = (row[f"{key}_{suffix}"] for key in ("estimate", "lower", "upper"))
            # Keep the native-PDF TAAE interval at its true width, even when it
            # is narrower than the plotted marker. The table preserves its sign.
            ax.errorbar(estimate, y, xerr=[[estimate - lower], [upper - estimate]],
                        fmt="o", color=color, markersize=4.5, elinewidth=1.5,
                        capsize=2.5, capthick=1, zorder=3)
            if suffix == "pp":
                ax.text(0.04, y, row["label"].replace(" minus ", "\nminus "), transform=label_transform,
                        va="center", ha="left", fontsize=8)
    fig.text(0.04, 0.488, "Mean paired differences: positive values indicate greater error in the first condition.",
             fontsize=8, color="#424B50")

    fig.text(0.04, 0.44, "C  Mean paired differences", fontsize=10.5, weight="bold")
    table = fig.add_axes((0.04, 0.088, 0.92, 0.33))
    table.set_xlim(0, 1)
    table.set_ylim(0, 7.4)
    table.axis("off")
    table.hlines([7.35, 6.0, 0.0], 0, 1, color="#46525A", linewidth=[0.8, 0.65, 0.8])
    table.text(0.01, 6.68, "Condition minus comparator", va="center", weight="bold", fontsize=7.8)
    table.text(0.525, 6.68, "TAAE (pp)\nMean [95% CI]", va="center", ha="center", weight="bold", fontsize=7.8)
    table.text(0.84, 6.68, "Absolute log-HR error\nMean [95% CI]", va="center", ha="center", weight="bold", fontsize=7.8)
    for index, row in enumerate(rows):
        y = 5.4 - 1.2 * index
        if index % 2 == 0:
            table.axhspan(y - 0.58, y + 0.58, color="#F2F5F7", zorder=-1)
        table.text(0.01, y, row["label"].replace(" minus ", "\nminus "), va="center", fontsize=7.5)
        taae = (f"{display_number(row['estimate_pp'])}\n"
                f"[{display_number(row['lower_pp'])}, {display_number(row['upper_pp'])}]")
        log_hr = (f"{row['estimate_log_hr']:.5f}\n"
                  f"[{row['lower_log_hr']:.5f}, {row['upper_log_hr']:.5f}]")
        table.text(0.525, y, taae, va="center", ha="center", fontsize=7.6, linespacing=1.35)
        table.text(0.84, y, log_hr, va="center", ha="center", fontsize=7.6, linespacing=1.35)
    fig.text(0.04, 0.052, "Without marks: censoring marks omitted. pp: percentage points. HR: hazard ratio.", fontsize=7.2)
    fig.text(0.04, 0.023, "In A, the native-PDF CI is narrower than the point; the estimate and CI are reported in C.", fontsize=7.2)
    output = FIGURES / "Figure4_figure_characteristics.pdf"
    fig.savefig(output, metadata={"CreationDate": None, "ModDate": None,
                                  "Title": "Paired effects of censoring marks, resolution, and compression"})
    plt.close(fig)
    return output

def followup_figure():
    source=read_csv('build/followup_contrasts.csv')
    contrasts=[('marks_vector','PDF without marks minus native PDF'),('marks_raster','PNG300 without marks minus PNG300'),
               ('resolution_150','PNG150 minus PNG300'),('resolution_100','PNG100 minus PNG300'),('compression_100','JPEG100 minus PNG100')]
    rows=[]
    for contrast,label in contrasts:
        item=dict(contrast=contrast,label=label,pairs=32)
        for metric,suffix,scale in [('taae_mean','pp',100),('absolute_log_hr_error','log_hr',1)]:
            row=next(r for r in source if r['contrast']==contrast and r['metric']==metric and r['scope']=='all')
            item.update({f'{k}_{suffix}':scale*float(row[k]) for k in ('estimate','lower','upper')})
        rows.append(item)
    out=draw_followup_figure(rows)
    with out.with_suffix('.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def supplementary_figures():
    three_pairs=pd.read_csv(ROOT/'build/image_only_pairs.csv',float_precision='round_trip')
    ext=pd.read_csv(ROOT/'build/image_only_bootstrap.csv',float_precision='round_trip')
    extension_primary=ext[(ext.metric=='taae_mean') & (ext.scope=='all')]
    template_pairs=pd.read_csv(ROOT/'build/survminer_pairs.csv',float_precision='round_trip')
    audit=read_csv('build/reconstruction_audit.csv')
    def flagged(experiment):return {r['dataset_id'] for r in audit if r['experiment']==experiment and r['confirmed_excess_corrections']=='True'}
    bad_datasets=flagged('primary');raster_excess=flagged('image_only');template_excess=flagged('survminer')
    new_difference=template_pairs.taae_mean_new_PDF-template_pairs.taae_mean_new_PNG300
    old_difference=template_pairs.taae_mean_old_PDF-template_pairs.taae_mean_old_PNG300
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'axes.titleweight':'bold','axes.labelsize':10,'legend.frameon':False})
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6), constrained_layout=True,
        gridspec_kw={'width_ratios': [1.1, 1]})
    formats = ['PDF', 'RASTER_PDF', 'PNG300']
    labels = ['Native PDF', 'Image-only PDF', 'PNG300']
    for i, row in three_pairs.iterrows():
        offset = (i % 8 - 3.5) * .013
        values = [row['taae_mean_' + fmt] for fmt in formats]
        axes[0].plot(np.arange(3) + offset, values, color='#b8b8b8', lw=.7, alpha=.7, zorder=1)
        for x, fmt, value in zip(np.arange(3) + offset, formats, values):
            flagged = (fmt == 'PNG300' and row.dataset in bad_datasets) or (
                fmt == 'RASTER_PDF' and row.dataset in raster_excess)
            axes[0].scatter(x, value, color=COLORS[fmt], s=23,
                marker='x' if flagged else 'o', zorder=3)
    axes[0].set_xticks(range(3), labels)
    axes[0].set_xlim(-.3, 2.3)
    axes[0].set_yscale('log')
    axes[0].set_ylabel('Mean TAAE across arms (probability scale)')
    axes[0].set_title('A  Same 32 datasets, three inputs')
    axes[0].grid(axis='y', alpha=.18)
    contrasts = ['PDF_minus_RASTER_PDF', 'RASTER_PDF_minus_PNG300']
    comparison_labels = ['Native PDF minus image-only PDF', 'Image-only PDF minus PNG300']
    for i, contrast in enumerate(contrasts):
        row = extension_primary[extension_primary.contrast == contrast].iloc[0]
        axes[1].errorbar(100*row.estimate, i,
            xerr=np.array([[100*(row.estimate-row.lower)], [100*(row.upper-row.estimate)]]),
            color='#333333', fmt='D', capsize=4, markersize=6)
    axes[1].axvline(0, color='#777777', lw=.9, ls='--')
    axes[1].set_yticks(range(2), comparison_labels)
    axes[1].set_ylim(-.5, 1.5)
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Difference in mean TAAE across arms (percentage points)')
    axes[1].set_title('B  Paired differences and 95% CIs')
    axes[1].xaxis.set_major_locator(MaxNLocator(nbins=5))
    axes[1].grid(axis='x', alpha=.18)
    savefig(fig, 'FigureS2_image_only')

    fig, axes = plt.subplots(1, 2, figsize=(10, 4.7), constrained_layout=True)
    for i, row in template_pairs.iterrows():
        offset = (i % 4 - 1.5) * .035
        values = [row.taae_mean_new_PDF, row.taae_mean_new_PNG300]
        axes[0].plot([offset, 1+offset], values, color='#b8b8b8', lw=.8, zorder=1)
        for x, fmt, value in zip([offset, 1+offset], ['PDF', 'PNG300'], values):
            axes[0].scatter(x, value, color=COLORS[fmt], s=29, zorder=3,
                marker='x' if fmt == 'PNG300' and row.dataset in template_excess else 'o')
        values = [100*old_difference.iloc[i], 100*new_difference.iloc[i]]
        axes[1].plot([offset, 1+offset], values, color='#b8b8b8', lw=.8, zorder=1)
        for x, value, flagged in zip([offset, 1+offset], values,
                [row.dataset in bad_datasets, row.dataset in template_excess]):
            axes[1].scatter(x, value, color='#555555', s=26, zorder=3,
                            marker='x' if flagged else 'o')
    axes[0].set_xticks([0,1], ['Native PDF', 'PNG300'])
    axes[0].set_yscale('log')
    axes[0].set(title='A  survminer: 16 paired datasets',
                ylabel='Mean TAAE across arms (probability scale)')
    axes[1].set_xticks([0,1], ['Original ReportLab', 'Later survminer'])
    axes[1].scatter([0,1], [100*old_difference.mean(),100*new_difference.mean()],
                    marker='D', s=65, color='#0072B2', zorder=4, label='Mean of 16 pairs')
    axes[1].axhline(0, color='#777777', lw=.9, ls='--')
    axes[1].set(title='B  Format differences: same 16 datasets',
                ylabel='PDF minus PNG300 TAAE (percentage points)')
    axes[1].legend(loc='lower left', fontsize=9)
    for ax in axes:
        ax.set_xlim(-.3,1.3)
        ax.grid(axis='y', alpha=.18)
    savefig(fig, 'FigureS1_plotting_style')

    three_pairs.to_csv(FIGURES/'FigureS2_plot_data.csv',index=False)
    template_pairs.to_csv(FIGURES/'FigureS1_plot_data.csv',index=False)

if __name__=='__main__':
    FIGURES.mkdir(parents=True,exist_ok=True)
    primary_figure();followup_figure();supplementary_figures()
    print('Generated Figures 2, 4, S1 and S2 from recalculated comparisons.')
