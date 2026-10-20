#!/usr/bin/env python3
"""Verify saved extraction classifications and recompute method diagnostics.

Inputs: data/run_manifest.csv, actual PDFs/reference IPDs, frozen reconstruction
IPDs and reconstruction_programs.zip (original source, decisions, replay receipts).
Outputs: build/reconstruction_audit.csv, extraction_strategy_counts.csv,
template_diagnostics.csv, template_arm_diagnostics.csv, method_audit.json.
Run: python scripts/audit_reconstruction_methods.py (repository root).
Classification labels are historical reviewed observations. This script verifies
the cited programs/anchors and aggregates those labels; it does not substitute an
automatic classifier or execute model-generated code. Geometry/KM diagnostics
are recalculated using the unchanged original diagnostic algorithms.
"""
from pathlib import Path
from collections import Counter, defaultdict
import csv
import json
import hashlib
import zipfile
import numpy as np
import pdfplumber

ROOT=Path(__file__).resolve().parents[1]
SPEC='docs/proposals/pdf_png300_extraction_strategy_audit_v1/'

def read_csv(path):
    with (ROOT/path).open(newline='') as stream: return list(csv.DictReader(stream))

def write_csv(path,rows):
    with (ROOT/path).open('w',newline='') as stream:
        w=csv.DictWriter(stream,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)

def km_steps(records: list[dict], group: str) -> dict[str, np.ndarray]:
    """KM after each unique event time; tied censors remain at risk for that event."""
    subset = [r for r in records if r["group"] == group]
    times = np.array([float(r["observed_time"]) for r in subset])
    events = np.array([int(r["event"]) for r in subset])
    require(len(times) > 0 and np.all(np.isfinite(times)), "Invalid diagnostic records")
    require(np.all(np.isin(events, [0, 1])), "Invalid event indicators")
    event_times, deaths = np.unique(times[events == 1], return_counts=True)
    risk = len(times) - np.searchsorted(np.sort(times), event_times, side="left")
    return {
        "time": event_times, "deaths": deaths, "risk": risk,
        "cumulative_events": np.cumsum(deaths),
        "height": np.cumprod(1.0 - deaths / risk),
        "individual_events": np.sort(times[events == 1]),
        "individual_censors": np.sort(times[events == 0]),
    }

def compare_arm(reference: list[dict], records: list[dict], group: str) -> dict:
    """Separate event timing from KM height; this is not a new accuracy score."""
    ref, rec = km_steps(reference, group), km_steps(records, group)
    require(len(ref["individual_events"]) == len(rec["individual_events"]),
            "Unequal event totals in template diagnostics")
    require(len(ref["individual_censors"]) == len(rec["individual_censors"]),
            "Unequal censor totals in template diagnostics")
    # Compare heights after the same cumulative number of events. Adjacent events
    # coalesced by figure precision are compared after the merged reconstructed drop.
    positions = np.searchsorted(ref["cumulative_events"], rec["cumulative_events"])
    require(np.all(ref["cumulative_events"][positions] == rec["cumulative_events"]),
            "Unmatched cumulative-event endpoints")
    within = rec["time"] <= 24
    height_delta = rec["height"][within] - ref["height"][positions][within]
    timing_delta = rec["individual_events"] - ref["individual_events"]
    censor_delta = rec["individual_censors"] - ref["individual_censors"]
    return {
        "events": len(timing_delta), "censors": len(censor_delta),
        "max_event_rank_time_difference_months": float(np.max(np.abs(timing_delta))),
        "median_signed_event_rank_time_difference_months": float(np.median(timing_delta)),
        "max_censor_rank_time_difference_months": float(np.max(np.abs(censor_delta))),
        "median_signed_censor_rank_time_difference_months": float(np.median(censor_delta)),
        "matched_endpoints_through_24_months": int(within.sum()),
        "height_mismatched_endpoints_tolerance_1e_12": int((np.abs(height_delta) > 1e-12).sum()),
        "maximum_matched_KM_height_difference_pp": float(np.max(np.abs(height_delta)) * 100),
    }

def geometry(path: Path) -> dict:
    with pdfplumber.open(path) as pdf:
        curves = [c for c in pdf.pages[0].curves if len(c["pts"]) > 10]
        require(len(curves) == 2, "Expected two native KM paths: " + str(path))
        points = np.concatenate([np.asarray(c["pts"], float) for c in curves])
    x = points[:, 0]
    # Cairo inputs show millipoint serialization. Also measure proximity to the
    # finer backend grid, without claiming which rendering operation caused it.
    return {
        "curve_vertices": len(x),
        "x_coordinates_on_0_001_pt_grid_fraction": float(np.mean(np.abs(x * 1000 - np.round(x * 1000)) < 1e-6)),
        "x_coordinates_within_0_000501_pt_of_1_over_256_grid_fraction":
            float(np.mean(np.abs(x - np.round(x * 256) / 256) <= 0.000501)),
    }

def glyph_center_bias(path: Path, reference: list[dict]) -> dict:
    """Compare raw PDF text bounding-box centers with observed censor times.

    This measures a diagnostic discrepancy; it does not modify coordinates or IPD.
    Axis ticks are fitted independently, with no import of saved reconstruction code.
    """
    colors = {"A": (0.972549, 0.462745, 0.427451), "B": (0.0, 0.74902, 0.768627)}
    results = {}
    with pdfplumber.open(path) as pdf:
        page = pdf.pages[0]
        xt = sorted(l["x0"] for l in page.lines if l["stroking_color"] == (0., 0., 0.)
                    and abs(l["x0"] - l["x1"]) < 1e-6 and 200 < l["top"] < 210)
        require(len(xt) == 7, "Diagnostic axis-tick count changed")
        ax, bx = np.polyfit(xt, np.arange(0, 37, 6), 1)
        for group, color in colors.items():
            centers = sorted((c["x0"] + c["x1"]) / 2 for c in page.chars
                             if c["text"] == "|" and c["non_stroking_color"] == color and 40 < c["top"] < 206)
            censors = np.sort([float(r["observed_time"]) for r in reference if r["group"] == group and r["event"] == "0"])
            require(len(centers) == len(censors), "Censor glyph/reference count differs")
            delta = ax * np.asarray(centers) + bx - censors
            results[group] = float(np.median(delta))
    return results


def main():
    runs=read_csv('data/run_manifest.csv');runmap={r['run_id']:r for r in runs}
    scores={r['run_id']:r for r in read_csv('build/scores.csv')}
    results=[];summaries=[];arms=[];checked_sources=set();anchors=0;replays=0
    with zipfile.ZipFile(ROOT/'frozen_results/reconstruction_programs.zip') as archive:
        def load(member): return json.loads(archive.read(member))
        decisions=load(SPEC+'CLASSIFICATIONS.json')['sessions']
        reviews={r['session_id']:r for r in load(SPEC+'template_source_review.json')['sessions']}
        def verify_source(entry):
            source=entry.get('file',entry.get('path'))
            data=archive.read(source)
            require(hashlib.sha256(data).hexdigest()==entry['sha256'],'Changed saved program: '+source)
            checked_sources.add(source)
            return data.decode().splitlines()
        for d in decisions:
            sid=d['session_id'];r=runmap[sid]
            for e in d['source_closure']: verify_source(e)
            for e in d['evidence']:
                lines=verify_source(e)
                require(lines[e['line']-1].strip()==e['excerpt'],'Changed source anchor: '+sid)
                anchors+=1
            complete=r['complete']=='True'
            results.append(dict(run_id=sid,experiment=r['experiment'],batch=r['batch'],dataset_id=r['dataset_id'],
                scenario=r['scenario'],condition=r['condition'],is_first_run=r['is_first_run'],
                complete=complete,recovery=r['assignment_status']=='recovery',
                strategy=d['curve_extraction_strategy'],taae_pp=100*float(scores[sid]['taae_mean']) if complete else '',
                classification_basis=d['classification_basis']))
        # Historical receipts prove replay at the recorded time; no fresh replay is claimed.
        for r in runs:
            if r['complete']!='True' or r['method']!='LLM' or r['experiment']=='effort_screen':continue
            receipt=load(r['replay_receipt_member'])
            digest=hashlib.sha256((ROOT/r['reconstruction_output']).read_bytes()).hexdigest()
            require(receipt['original_csv_sha256']==receipt['replay_csv_sha256']==digest,'Changed replay receipt: '+r['run_id'])
            replays+=1
        byid={r['run_id']:r for r in results}
        originals={r['dataset_id']:r for r in runs if r['primary64']=='True' and r['batch']=='first' and r['condition']=='PDF'}
        for sid,review in reviews.items():
            for e in review['source_evidence']:
                lines=verify_source(e)
                require('\n'.join(lines[e['start_line']-1:e['end_line']])==e['excerpt'],'Changed targeted source anchor')
                anchors+=1
            item=runmap[sid];ds=item['dataset_id'];original=originals[ds]
            reference=read_csv(item['reference_ipd'])
            result=dict(dataset=ds,scenario=item['scenario'],n_per_group=item['n_per_arm'],survminer_session=sid,
                        original_session=original['run_id'],original_TAAE_pp=byid[original['run_id']]['taae_pp'],
                        survminer_TAAE_pp=byid[sid]['taae_pp'],original_strategy=byid[original['run_id']]['strategy'],
                        survminer_strategy=byid[sid]['strategy'])
            result.update({k:v for k,v in review.items() if k not in ('session_id','source_evidence')})
            result.update({'raw_glyph_center_median_censor_time_bias_'+g+'_months':v for g,v in glyph_center_bias(ROOT/item['input_figure'],reference).items()})
            for label,r in [('original',original),('survminer',item)]:
                result.update({label+'_'+k:v for k,v in geometry(ROOT/r['input_figure']).items()})
                metrics=[compare_arm(reference,read_csv(r['reconstruction_output']),g) for g in 'AB']
                for g,values in zip('AB',metrics):arms.append(dict(dataset=ds,template=label,session_id=r['run_id'],group=g,**values))
                result[label+'_max_event_rank_time_difference_months']=max(m['max_event_rank_time_difference_months'] for m in metrics)
                result[label+'_median_censor_time_bias_A_months']=metrics[0]['median_signed_censor_rank_time_difference_months']
                result[label+'_median_censor_time_bias_B_months']=metrics[1]['median_signed_censor_rank_time_difference_months']
                result[label+'_max_matched_height_difference_pp']=max(m['maximum_matched_KM_height_difference_pp'] for m in metrics)
                result[label+'_height_matches_both_arms_tolerance_1e_12']=all(m['height_mismatched_endpoints_tolerance_1e_12']==0 for m in metrics)
            summaries.append(result)
        excess=set()
        for study,filename in [('pdf_png300_32_budget_v2','session_review_final.json'),('pdf_png300_raster_pdf_extension_v1','session_review.json'),('pdf_png300_survminer_template_v1','session_review.json')]:
            for row in load('outputs/evaluation/'+study+'/method_audit_v1/'+filename)['sessions']:
                if row.get('confirmed_excess_primary_corrections'):excess.add(row['session_id'])
        excess.update(r['session_id'] for r in load('outputs/evaluation/pdf_png300_replication32_v1/method_audit_v1/process_adjudications.json')['confirmed_excess_correction_sessions'])
    for r in results:r['confirmed_excess_corrections']=r['run_id'] in excess
    grouped=defaultdict(list)
    for r in results:grouped[(r['experiment'],r['condition'])].append(r)
    counts=[]
    for (experiment,condition),members in grouped.items():
        counts.append(dict(experiment=experiment,condition=condition,registered=len(members),
            complete=sum(r['complete'] for r in members),
            path_parsing=sum(r['strategy']=='path_parsing' for r in members),
            pixel_segmentation=sum(r['strategy']=='pixel_segmentation' for r in members),
            hybrid=sum(r['strategy']=='hybrid' for r in members),
            missing=sum(not r['complete'] for r in members),
            complete_first=sum(r['complete'] and r['is_first_run']=='True' for r in members),
            complete_repeat=sum(r['complete'] and r['is_first_run']=='False' for r in members)))
    pilot=[r for r in runs if r['experiment']=='effort_screen']
    summary=dict(registrations=len(results),complete=sum(r['complete'] for r in results),verified_program_files=len(checked_sources),
        verified_source_anchors=anchors,verified_historical_replay_receipts=replays,
        strategy_counts=dict(Counter(r['strategy'] for r in results)),
        excess_corrections_by_strategy={s:sum(r['run_id'] in excess for r in results if r['strategy']==s) for s in ('path_parsing','pixel_segmentation')},
        exceeded_run_ids=sorted(excess),
        effort_screen={e:dict(assigned=sum(r['effort']==e for r in pilot),complete=sum(r['effort']==e and r['complete']=='True' for r in pilot)) for e in ('medium','high')},
        replay_scope='Existing receipts checked against frozen IPD; programs not freshly executed.')
    require(len(results)==246 and summary['complete']==240 and replays==400,'Audit population changed')
    write_csv('build/reconstruction_audit.csv',results)
    write_csv('build/extraction_strategy_counts.csv',counts)
    write_csv('build/template_diagnostics.csv',summaries)
    write_csv('build/template_arm_diagnostics.csv',arms)
    (ROOT/'build/method_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('Verified source anchors and 400 historical replays; recalculated 16 matched-template diagnostics.')

if __name__=='__main__':main()
