#!/usr/bin/env python3
"""Calculate Cox/Greenwood uncertainty and event-count errors from saved IPD.

Inputs: data/run_manifest.csv and its original/reconstructed CSV paths.
Outputs: build/uncertainty.csv and build/reference_uncertainty.json.
Run: python scripts/calculate_uncertainty_and_events.py (repository root).
Algorithms: unchanged prepare_se_estimands.py kernels; Efron ties without
near-time snapping, Newton tolerance 1e-12, 100 iterations, tau=24 months.
Relative SE errors are fractions; Greenwood errors are probabilities.
"""
from pathlib import Path
import csv
import json
import math
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
TAU = 24.0
METRICS = ["abs_log_hr_error", "abs_se_log_hr_error", "rel_se_log_hr_error", "abs_rel_se_log_hr_error",
           "abs_greenwood_se24_error_A", "abs_greenwood_se24_error_B",
           "signed_events_total_error_A", "signed_events_total_error_B",
           "abs_events_total_error_A", "abs_events_total_error_B",
           "abs_events_24_error_A", "abs_events_24_error_B"]

def read_records(path):
    rows = []
    with open(path, newline="") as stream:
        for r in csv.DictReader(stream):
            rows.append((r["group"], r["observed_time"], r["event"]))
    return rows

def arm_arrays(rows, group):
    t, e = [], []
    for g, time, event in rows:
        if g != group:
            continue
        try:
            tv, ev = float(time), float(event)
        except (TypeError, ValueError):
            return None
        if not (math.isfinite(tv) and tv >= 0 and ev in (0.0, 1.0)):
            return None
        t.append(tv); e.append(int(ev))
    if not t:
        return None
    return np.array(t), np.array(e)

def km_at(time, event, tau=TAU):
    """KM estimate at tau with Greenwood SE; returns dict with support flag."""
    order = np.argsort(time, kind="stable")
    t, e = time[order], event[order]
    ev_times = np.unique(t[(e == 1) & (t <= tau)])
    s, var_sum = 1.0, 0.0
    sorted_t = t  # ascending
    for ti in ev_times:
        n = len(sorted_t) - np.searchsorted(sorted_t, ti, side="left")
        d = int(np.sum((t == ti) & (e == 1)))
        s *= 1.0 - d / n
        if n - d > 0:
            var_sum += d / (n * (n - d))
        else:
            var_sum = math.inf
    supported = bool(t.max() >= tau or s == 0.0)
    se = 0.0 if s == 0.0 else (math.sqrt(s * s * var_sum) if math.isfinite(var_sum) else math.nan)
    return dict(S=s, se=se, supported=supported,
                events_total=int(e.sum()), events_tau=int(((e == 1) & (t <= tau)).sum()), n=int(len(t)))

def cox_efron_binary(time, event, x, max_iter=100, tol=1e-12):
    """Cox PH with one binary covariate, Efron ties; returns beta, se, status.

    Matches survival::coxph(..., ties="efron", timefix=FALSE) for exact ties.
    """
    time = np.asarray(time, float); event = np.asarray(event, int); x = np.asarray(x, int)
    if event.sum() == 0:
        return dict(beta=math.nan, se=math.nan, status="no_events")
    if len(np.unique(x[event == 1])) < 2:
        return dict(beta=math.nan, se=math.nan, status="events_in_one_group")
    ev_times = np.unique(time[event == 1])
    t0 = np.sort(time[x == 0]); t1 = np.sort(time[x == 1])
    r0 = len(t0) - np.searchsorted(t0, ev_times, side="left")  # at risk, x=0
    r1 = len(t1) - np.searchsorted(t1, ev_times, side="left")  # at risk, x=1
    d0 = np.array([int(np.sum((time == ti) & (event == 1) & (x == 0))) for ti in ev_times])
    d1 = np.array([int(np.sum((time == ti) & (event == 1) & (x == 1))) for ti in ev_times])
    d = d0 + d1
    # Efron fractions j/d for j = 0..d-1, expanded per event time
    frac = np.concatenate([np.arange(di) / di for di in d])
    idx = np.repeat(np.arange(len(d)), d)

    def derivatives(beta):
        eb = math.exp(beta)
        s0 = r0 + r1 * eb; s1 = r1 * eb
        e0 = d0 + d1 * eb; e1 = d1 * eb
        den = s0[idx] - frac * e0[idx]
        num = s1[idx] - frac * e1[idx]
        ratio = num / den
        loglik = float(beta * d1.sum() - np.log(den).sum())
        score = float(d1.sum() - ratio.sum())
        info = float((ratio - ratio ** 2).sum())
        return loglik, score, info

    beta = 0.0
    loglik, score, info = derivatives(beta)
    for _ in range(max_iter):
        if info <= 0 or not math.isfinite(info):
            return dict(beta=math.nan, se=math.nan, status="nonfinite_information")
        step = score / info
        new = beta + step
        # step halving to guarantee ascent
        new_ll, new_score, new_info = derivatives(new)
        halvings = 0
        while new_ll < loglik - 1e-12 and halvings < 30:
            step /= 2; new = beta + step
            new_ll, new_score, new_info = derivatives(new)
            halvings += 1
        beta, loglik, score, info = new, new_ll, new_score, new_info
        if abs(step) < tol and abs(score) < 1e-9:
            break
    else:
        if abs(score) > 1e-6:
            return dict(beta=beta, se=math.nan, status="not_converged")
    if not math.isfinite(beta) or abs(beta) > 20:
        return dict(beta=beta, se=math.nan, status="nonfinite_or_unidentified")
    return dict(beta=float(beta), se=float(1.0 / math.sqrt(info)), status="available")

def analyze(rows, arms_ok):
    """Per-record-set quantities; arms_ok maps 'A'/'B' to validity flags."""
    out = dict(status="available")
    arrays = {}
    for g in ("A", "B"):
        arr = arm_arrays(rows, g) if arms_ok.get(g, True) else None
        if arr is None:
            out[f"km_{g}"] = None
            out["status"] = "missing_or_invalid_arm"
            continue
        arrays[g] = arr
        out[f"km_{g}"] = km_at(*arr)
    if len(arrays) == 2:
        time = np.concatenate([arrays["A"][0], arrays["B"][0]])
        event = np.concatenate([arrays["A"][1], arrays["B"][1]])
        x = np.concatenate([np.zeros(len(arrays["A"][0]), int), np.ones(len(arrays["B"][0]), int)])
        out["cox"] = cox_efron_binary(time, event, x)
    else:
        out["cox"] = dict(beta=math.nan, se=math.nan, status="missing_or_invalid_arm")
    return out

def errors(rec, ref):
    """Reconstructed-minus-reference error metrics; NaN where unavailable."""
    e = {m: math.nan for m in METRICS}
    e["log_hr_ref"] = ref["cox"]["beta"]; e["log_hr_rec"] = rec["cox"]["beta"]
    e["se_log_hr_ref"] = ref["cox"]["se"]; e["se_log_hr_rec"] = rec["cox"]["se"]
    if rec["cox"]["status"] == "available" and ref["cox"]["status"] == "available":
        e["abs_log_hr_error"] = abs(rec["cox"]["beta"] - ref["cox"]["beta"])
        e["abs_se_log_hr_error"] = abs(rec["cox"]["se"] - ref["cox"]["se"])
        e["rel_se_log_hr_error"] = (rec["cox"]["se"] - ref["cox"]["se"]) / ref["cox"]["se"]
        e["abs_rel_se_log_hr_error"] = abs(e["rel_se_log_hr_error"])
    for g in ("A", "B"):
        kr, kx = ref[f"km_{g}"], rec[f"km_{g}"]
        e[f"greenwood_se24_ref_{g}"] = kr["se"] if kr else math.nan
        e[f"greenwood_se24_rec_{g}"] = kx["se"] if kx else math.nan
        e[f"events_total_ref_{g}"] = kr["events_total"] if kr else math.nan
        e[f"events_total_rec_{g}"] = kx["events_total"] if kx else math.nan
        if kr and kx and kr["supported"] and kx["supported"] and math.isfinite(kx["se"]):
            e[f"abs_greenwood_se24_error_{g}"] = abs(kx["se"] - kr["se"])
        if kr and kx:
            e[f"signed_events_total_error_{g}"] = kx["events_total"] - kr["events_total"]
            e[f"abs_events_total_error_{g}"] = abs(kx["events_total"] - kr["events_total"])
            e[f"abs_events_24_error_{g}"] = abs(kx["events_tau"] - kr["events_tau"])
    return e


def main():
    with (ROOT / 'data/run_manifest.csv').open() as stream:
        runs = list(csv.DictReader(stream))
    runs = [r for r in runs if r['complete']=='True' and r['experiment']!='effort_screen']
    references, output = {}, []
    for r in runs:
        dataset = r['dataset_id']
        if dataset not in references:
            references[dataset] = analyze(read_records(ROOT/r['reference_ipd']), {'A':True,'B':True})
        rec = analyze(read_records(ROOT/r['reconstruction_output']), {'A':True,'B':True})
        values = errors(rec, references[dataset])
        if not all(math.isfinite(v) for v in values.values()):
            raise ValueError(f"Unavailable uncertainty estimate: {r['run_id']}")
        output.append({k:r[k] for k in ('run_id','dataset_id','scenario','batch','condition','experiment','method')} | values)
    (ROOT/'build').mkdir(exist_ok=True)
    with (ROOT/'build/uncertainty.csv').open('w',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=list(output[0]));writer.writeheader();writer.writerows(output)
    (ROOT/'build/reference_uncertainty.json').write_text(json.dumps(references,indent=2,allow_nan=False)+'\n')
    print(f'Calculated uncertainty and event counts for {len(output)} saved reconstructions.')


if __name__=='__main__':
    main()
