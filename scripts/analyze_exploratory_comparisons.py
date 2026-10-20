#!/usr/bin/env python3
"""Recompute uncertainty, exact-coordinate, figure-characteristic and sensitivity comparisons.

Inputs: build/scores.csv, uncertainty.csv, primary_pairs.csv, historical62_pairs.csv.
Outputs: build/exploratory_contrasts.csv, followup_contrasts.csv,
followup_conditions.csv, revision_sensitivity.json, missing_bounds.json.
Run: python scripts/analyze_exploratory_comparisons.py (repository root).
Preserves NumPy PCG64 seeds, assigned-dataset resampling, percentile type 1,
and historical ordering. Current64 uses four strata of 16. Follow-up uses
four of eight. Revision sensitivity uses eight draws per retained batch/scenario
stratum, even when that stratum contains fewer than eight retained datasets.
"""
from collections import defaultdict
from pathlib import Path
import csv
import json
import math
import statistics as st
import numpy as np
from calculate_uncertainty_and_events import METRICS as SE_METRICS

ROOT = Path(__file__).resolve().parents[1]
B = 10000

def read_csv(path):
    with (ROOT/path).open(newline='') as stream:
        return list(csv.DictReader(stream))

def write_csv(path, rows):
    with (ROOT/path).open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def bootstrap(rows, metric, left, right, seed, scenario_only):
    grouped = defaultdict(list)
    for row in rows:
        key = row['scenario'] if scenario_only else row['stratum']
        grouped[key].append(row[f'{metric}_{left}'] - row[f'{metric}_{right}'])
    size = 16 if scenario_only else 8
    arrays = np.array([grouped[key] for key in sorted(grouped)])
    assert arrays.shape == (64 // size, size) and np.isfinite(arrays).all()
    # The flattened sequence is identical to the historical draw/stratum loops
    # because all strata have the same size and RNG integer bounds.
    rng = np.random.default_rng(seed)
    picks = rng.integers(0, size, size=(B, len(arrays), size))
    draws = arrays[np.arange(len(arrays))[None, :, None], picks].mean(axis=2).mean(axis=1)
    ordered = np.sort(draws)
    return dict(estimate=float(arrays.mean(axis=1).mean()),
                lower=float(ordered[int(np.ceil(B * .025)) - 1]),
                upper=float(ordered[int(np.ceil(B * .975)) - 1]),
                standard_error=float(draws.std(ddof=1)),
                status='available' if ordered[0] != ordered[-1] else 'degenerate_resamples')

CURVE = [("taae_mean", "TAAE, both arms (pp)", 100),
         ("maximum_A", "Maximum curve error, A (pp)", 100),
         ("maximum_B", "Maximum curve error, B (pp)", 100),
         ("duration_above_002_A", "Time above 2 pp, A (months)", 1),
         ("duration_above_002_B", "Time above 2 pp, B (months)", 1),
         ("absolute_rmst_difference_error", "Absolute RMST-contrast error (months)", 1),
         ("absolute_log_hr_error", "Absolute log-HR error", 1)]
EXCEEDED = {"PP32B_R009", "PP32B_R013", "PP32B_R025", "PP32B_R026", "PP32B_R047", "PP32B_R062",
            "RP32_R010", "RP32_R023", "RP32_R031", "RP32_R033", "RP32_R052", "RP32_R061"}
FOLLOW_SCENARIOS = ("S1", "S2", "S3", "S4")
FOLLOW_SCOPES = {"all": FOLLOW_SCENARIOS, "n100": ("S1", "S2"), "n1000": ("S3", "S4"),
          "noncrossing": ("S1", "S3"), "crossing": ("S2", "S4"),
          "S1": ("S1",), "S2": ("S2",), "S3": ("S3",), "S4": ("S4",)}
FOLLOW_CONTRASTS = (("marks_vector", "NOMARK_PDF", "PDF"), ("marks_raster", "NOMARK_PNG300", "PNG300"),
             ("resolution_150", "PNG150", "PNG300"), ("resolution_100", "PNG100", "PNG300"),
             ("compression_100", "JPEG100", "PNG100"))
FOLLOW_CONDITIONS = ("PDF", "PNG300", "NOMARK_PDF", "NOMARK_PNG300", "PNG150", "PNG100", "JPEG100")
FOLLOW_METRICS = ['taae_mean', 'taae_A', 'taae_B', 'maximum_A', 'maximum_B', 'duration_above_002_A', 'duration_above_002_B', 'absolute_rmst_error_A', 'absolute_rmst_error_B', 'absolute_rmst_difference_error', 'absolute_log_hr_error', 'abs_rel_se_log_hr_error', 'signed_events_total_error_A', 'signed_events_total_error_B']
FOLLOW_SEED = 2026100406

def quantile_type1(x, p):
    x = np.sort(x)
    return float(x[min(max(1, math.ceil(len(x) * p)), len(x)) - 1])

def followup_bootstrap(table, datasets):
    """Prespecified contrasts; whole datasets resampled within scenario, indices shared by all outputs."""
    rng = np.random.Generator(np.random.PCG64(FOLLOW_SEED))
    picks = {s: rng.integers(0, 8, size=(B, 8)) for s in FOLLOW_SCENARIOS}
    results = []
    for name, new, old in FOLLOW_CONTRASTS:
        for metric in FOLLOW_METRICS:
            means, draws, complete = {}, {}, {}
            for s in FOLLOW_SCENARIOS:
                d = np.array([table[(ds, new)][metric] - table[(ds, old)][metric] for ds in datasets[s]], float)
                finite = np.isfinite(d)
                complete[s] = int(finite.sum())
                means[s] = float(d[finite].mean()) if finite.any() else math.nan
                sampled = d[picks[s]]
                counts = np.isfinite(sampled).sum(axis=1)
                totals = np.where(np.isfinite(sampled), sampled, 0.0).sum(axis=1)
                draws[s] = np.where(counts > 0, totals / np.maximum(counts, 1), np.nan)
            for scope, members in FOLLOW_SCOPES.items():
                estimate = float(np.mean([means[s] for s in members]))
                x = np.mean([draws[s] for s in members], axis=0)
                undefined = int((~np.isfinite(x)).sum())
                degenerate = not undefined and np.unique(x).size == 1
                status = ("undefined_estimate" if not math.isfinite(estimate) else "undefined_resamples"
                          if undefined else "degenerate_resamples" if degenerate else "available")
                ok = status == "available"
                results.append(dict(contrast=name, new_condition=new, comparator=old, metric=metric, scope=scope,
                    estimate=estimate, complete_pairs=sum(complete[s] for s in members),
                    assigned_pairs=8 * len(members),
                    lower=quantile_type1(x, .025) if ok else math.nan,
                    upper=quantile_type1(x, .975) if ok else math.nan,
                    standard_error=float(np.std(x, ddof=1)) if ok else math.nan,
                    undefined_resamples=undefined, status=status, B=B, seed=FOLLOW_SEED))
    return results

def describe_followup(values):
    x = np.array([v for v in values if math.isfinite(v)], float)
    q = np.percentile(x, [0, 25, 50, 75, 100])
    return dict(n=int(x.size), mean=float(x.mean()), minimum=q[0], q1=q[1], median=q[2], q3=q[3], maximum=q[4])

def q1(x, p):
    x = np.sort(np.asarray(x, float)); n = len(x)
    return float(x[min(max(1, math.ceil(n * p)), n) - 1])

def revision_contrast(rows, seed):
    strata = defaultdict(list)
    for r in rows:
        strata[r["stratum"]].append(r["diff"])
    arrays = [np.array(strata[s]) for s in sorted(strata)]
    if not arrays:
        return dict(estimate=math.nan, status="undefined")
    est = float(np.mean([a.mean() for a in arrays]))
    rng = np.random.default_rng(seed)
    draws = np.array([np.mean([a[rng.integers(0, len(a), size=8)].mean() for a in arrays]) for _ in range(B)])
    return dict(estimate=est, lower=q1(draws, .025), upper=q1(draws, .975), pairs=len(rows), strata=len(arrays),
                min_per_stratum=int(min(len(a) for a in arrays)), status="available")

def pp(value):
    return 100.0 * float(value)

def stratified_mean(rows, fill_pp=None):
    """Equal-scenario-weight mean of PDF-minus-PNG300 TAAE differences (pp).

    Complete pairs use the exported difference. If fill_pp is given, each
    PDF-only pair contributes PDF TAAE minus fill_pp.
    """
    strata = defaultdict(list)
    for row in rows:
        if row["pair_status"] == "both":
            strata[row["scenario"]].append(pp(row["taae_difference"]))
        elif fill_pp is not None and row["status_PDF"] == "complete_scorable":
            strata[row["scenario"]].append(pp(row["taae_mean_PDF"]) - fill_pp)
    return st.mean(st.mean(strata[s]) for s in sorted(strata))

def describe_exclusive(values):
    q = st.quantiles(values, n=4, method="exclusive")  # R type 6
    return dict(n=len(values), mean=st.mean(values), median=st.median(values),
                q1=q[0], q3=q[2], min=min(values), max=max(values))

def missing_bounds():
    current = read_csv('build/primary_pairs.csv')
    historical = read_csv('build/historical62_pairs.csv')
    original = [r for r in current if r['dataset'].startswith('PP32')]
    second = [r for r in historical if r['dataset'].startswith('RP32')]
    transport = [r for r in current if r['dataset'].startswith('RP32')]
    png = lambda rows: [pp(r["taae_mean_PNG300"]) for r in rows if r["pair_status"] == "both"]
    complete_png = png(original) + png(second)
    low, high = min(complete_png), max(complete_png)
    missing = [r for r in second if r["pair_status"] != "both"]
    out = dict(
        units="percentage points of survival probability",
        original=dict(difference=stratified_mean(original), png300=describe_exclusive(png(original))),
        second_main=dict(difference=stratified_mean(second), png300=describe_exclusive(png(second)),
                         missing=[dict(dataset=r["dataset"], scenario=r["scenario"], status_PNG300=r["status_PNG300"],
                                       tipping_point_pdf_taae=pp(r["taae_mean_PDF"])) for r in missing],
                         fill_observed_min=dict(value=low, difference=stratified_mean(second, low)),
                         fill_observed_max=dict(value=high, difference=stratified_mean(second, high)),
                         fill_zero=dict(value=0.0, difference=stratified_mean(second, 0.0))),
        second_transport=dict(difference=stratified_mean(transport), png300=describe_exclusive(png(transport)),
                              recovered=[dict(dataset=r["dataset"], taae_PNG300=pp(r["taae_mean_PNG300"]))
                                         for r in transport if r["dataset"] in {m["dataset"] for m in missing}]),
        pooled62=dict(png300=describe_exclusive(complete_png),
                      difference=(stratified_mean(original) + stratified_mean(second)) / 2,
                      fill_observed_min=(stratified_mean(original) + stratified_mean(second, low)) / 2,
                      fill_observed_max=(stratified_mean(original) + stratified_mean(second, high)) / 2,
                      fill_zero=(stratified_mean(original) + stratified_mean(second, 0.0)) / 2),
        pooled64_transport=dict(png300=describe_exclusive(png(original) + png(transport)),
                                difference=(stratified_mean(original) + stratified_mean(transport)) / 2),
        note="Pooled values combine equal-weight cohort means; the frozen pooled bootstrap used the same weighting.",
    )
    (ROOT/'build/missing_bounds.json').write_text(json.dumps(out,indent=2)+'\n')
    return out


def primary_comparisons(scores, uncertainty):
    primary = [r for r in scores if r['primary64']=='TRUE']
    by_dataset = {}
    # Preserve the original run first-appearance order used for SE bootstrapping.
    for r in primary:
        row=by_dataset.setdefault(r['dataset_id'],dict(dataset=r['dataset_id'],scenario=r['scenario'],stratum=r['batch']+'.'+r['scenario']))
        for metric in SE_METRICS:
            row[f"{metric}_{r['condition']}"]=float(uncertainty[r['run_id']][metric])
    for r in scores:
        if r['method'] in ('Guyot','GuyotExact'):
            for metric in SE_METRICS:
                by_dataset[r['dataset_id']][f"{metric}_{r['method']}"]=float(uncertainty[r['run_id']][metric])
    # The prior exporter consumed four scopes, each with 12*5 contrasts, before
    # the completed64 scope. Preserve that original seed offset exactly.
    seed=2026100401 + 4*len(SE_METRICS)*5
    output=[]
    def add(family, rows, metric, left, right, seed):
        result=bootstrap(rows,metric,left,right,seed,True)
        output.append(dict(family=family,metric=metric,contrast=left+'_minus_'+right,
                           mean_left=float(np.mean([r[f'{metric}_{left}'] for r in rows])),
                           mean_right=float(np.mean([r[f'{metric}_{right}'] for r in rows])),
                           pairs=len(rows),B=B,seed=seed,
                           **{'new_'+k:v for k,v in result.items()}))
    for metric in SE_METRICS:
        for left,right in [('PDF','PNG300'),('Guyot','PNG300'),('Guyot','GuyotExact'),('GuyotExact','PDF'),('GuyotExact','PNG300')]:
            add('precision',list(by_dataset.values()),metric,left,right,seed);seed+=1
    # Curve decomposition originally consumed the dataset-sorted R pair table.
    pairs=read_csv('build/primary_pairs.csv')
    lookup={(r['dataset_id'],r['method']):r for r in scores if r['method']!='LLM'}
    curve_rows=[]
    for p in pairs:
        row=dict(dataset=p['dataset'],scenario=p['scenario'],stratum=p['cohort']+'.'+p['scenario'])
        for metric,_,scale in CURVE:
            for method in ('PDF','PNG300'):
                row[f'{metric}_{method}']=float(p[f'{metric}_{method}'])*scale
            for method in ('Guyot','GuyotExact'):
                row[f'{metric}_{method}']=float(lookup[(p['dataset'],method)][metric])*scale
        curve_rows.append(row)
    seed=2026100402
    for metric,_,scale in CURVE:
        for left,right in [('Guyot','GuyotExact'),('GuyotExact','PDF'),('GuyotExact','PNG300')]:
            add('decomposition',curve_rows,metric,left,right,seed);seed+=1
    write_csv('build/exploratory_contrasts.csv',output)


def figure_conditions(scores,uncertainty):
    chosen=[r for r in scores if r['experiment']=='followup' or
            (r['primary64']=='TRUE' and r['batch']=='first')]
    assert len(chosen)==224
    table={}
    for r in chosen:
        table[(r['dataset_id'],r['condition'])] = dict(
            session_id=r['run_id'],scenario=r['scenario'],
            **{m:float(r[m]) if m in r else float(uncertainty[r['run_id']][m]) for m in FOLLOW_METRICS})
    datasets={s:sorted({ds for ds,c in table if table[(ds,c)]['scenario']==s}) for s in FOLLOW_SCENARIOS}
    assert all(len(v)==8 for v in datasets.values())
    write_csv('build/followup_contrasts.csv',followup_bootstrap(table,datasets))
    summaries=[]
    for condition in FOLLOW_CONDITIONS:
        cells=[table[(ds,condition)] for s in FOLLOW_SCENARIOS for ds in datasets[s]]
        for scope,members in FOLLOW_SCOPES.items():
            c=[v for v in cells if v['scenario'] in members]
            stats=describe_followup([v['taae_mean'] for v in c])
            summaries.append(dict(condition=condition,scope=scope,**{'taae_'+k:v for k,v in stats.items()},
                se_log_hr_within_1pct=sum(v['abs_rel_se_log_hr_error']<=.01 for v in c),
                se_log_hr_within_2pct=sum(v['abs_rel_se_log_hr_error']<=.02 for v in c),
                max_abs_rel_se_log_hr_error=max(v['abs_rel_se_log_hr_error'] for v in c),
                exact_event_totals_both_arms=sum(v['signed_events_total_error_A']==0 and v['signed_events_total_error_B']==0 for v in c),
                datasets=len(c),maximum_curve_error=max(max(v['maximum_A'],v['maximum_B']) for v in c),
                mean_absolute_log_hr_error=float(np.mean([v['absolute_log_hr_error'] for v in c]))))
    write_csv('build/followup_conditions.csv',summaries)
    write_csv('build/followup_dataset_conditions.csv',[
        dict(dataset=ds,condition=condition,**values) for (ds,condition),values in table.items()])


def revision_sensitivity(scores):
    rows=[]
    session_of={(r['dataset_id'],r['condition']):r['run_id'] for r in scores if r['primary64']=='TRUE'}
    for r in read_csv('build/historical62_pairs.csv'):
        if r['pair_status']!='both': continue
        sid=session_of[(r['dataset'],'PNG300')]
        cohort='original' if r['dataset'].startswith('PP32') else 'second'
        rows.append(dict(cohort=cohort,dataset=r['dataset'],stratum=cohort+'.'+r['scenario'],session=sid,
                         exceeded=sid in EXCEEDED,png=100*float(r['taae_mean_PNG300']),
                         pdf=100*float(r['taae_mean_PDF']),diff=100*float(r['taae_difference'])))
    exc=[r for r in rows if r['exceeded']];non=[r for r in rows if not r['exceeded']]
    def desc(v):
        v=np.asarray(v);return dict(n=len(v),mean=float(v.mean()),median=float(np.median(v)),min=float(v.min()),max=float(v.max()))
    out=dict(seed=2026100403,B=B,exceeded_sessions_matched=sorted(r['session'] for r in exc),n_pairs=len(rows),
             png300_taae_pp=dict(exceeded=desc([r['png'] for r in exc]),not_exceeded=desc([r['png'] for r in non])),
             contrast_all=revision_contrast(rows,2026100403),contrast_excluding_exceeded=revision_contrast(non,2026100404),
             contrast_within_exceeded=revision_contrast(exc,2026100405),all_pairs_favor_pdf_excluding=all(r['diff']<0 for r in non))
    (ROOT/'build/revision_sensitivity.json').write_text(json.dumps(out,indent=2)+'\n')


def main():
    scores=read_csv('build/scores.csv')
    uncertainty={r['run_id']:r for r in read_csv('build/uncertainty.csv')}
    primary_comparisons(scores,uncertainty)
    figure_conditions(scores,uncertainty)
    revision_sensitivity(scores)
    missing_bounds()
    print('Computed 81 current64 Python contrasts, figure conditions, revision50 and missing-output bounds.')

if __name__=='__main__': main()
