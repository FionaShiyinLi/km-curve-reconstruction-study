#!/usr/bin/env python3
"""Export all 3 main and 11 supplementary tables from newly calculated results.

Inputs: data/run_manifest.csv and build/*.csv/json from the analysis scripts.
Outputs: build/tables/{Table1..3,TableS1..11}.{csv,tex}, tables.html, and
build/reported_statistics.json. CSVs retain unrounded numerical values; TeX/HTML
are readable views. Run: python scripts/write_manuscript_tables.py (repo root).
No published table or frozen validation output is read by this script.
"""
from pathlib import Path
from collections import Counter
import csv
import html
import json
import math
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'build/tables'
CURVE=[('taae_mean','Mean TAAE, both arms','pp'),('taae_A','TAAE, A','pp'),('taae_B','TAAE, B','pp'),
       ('maximum_A','Maximum curve error, A','pp'),('maximum_B','Maximum curve error, B','pp'),
       ('duration_above_002_A','Time above 2 pp, A','months'),('duration_above_002_B','Time above 2 pp, B','months'),
       ('absolute_rmst_error_A','Absolute RMST error, A','months'),('absolute_rmst_error_B','Absolute RMST error, B','months'),
       ('absolute_rmst_difference_error','Absolute RMST-contrast error','months'),('absolute_log_hr_error','Absolute log-HR error','log-HR')]
SE=[('abs_log_hr_error','Log HR, absolute error'),('abs_se_log_hr_error','SE(log HR), absolute error'),
    ('rel_se_log_hr_error','SE(log HR), signed relative error'),('abs_rel_se_log_hr_error','SE(log HR), absolute relative error'),
    ('abs_greenwood_se24_error_A','Greenwood SE at 24 months, A'),('abs_greenwood_se24_error_B','Greenwood SE at 24 months, B'),
    ('signed_events_total_error_A','Events, A, signed error'),('signed_events_total_error_B','Events, B, signed error'),
    ('abs_events_total_error_A','Events, A, absolute error'),('abs_events_total_error_B','Events, B, absolute error')]
TITLES={'Table1':'Analysis overview','Table2':'Native PDF versus PNG300: current 64 pairs',
        'Table3':'Automated Guyot versus LLM PNG300: 64 pairs','TableS1':'Dataset-level paired KM errors: probability scale',
        'TableS2':'Direct first–repeat discrepancy: probability scale','TableS3':'Image-only PDF contrasts',
        'TableS4':'Automated Guyot and LLM accuracy','TableS5':'Native-PDF plotting-style diagnostics',
        'TableS6':'Missing-output sensitivity analyses','TableS7':'Figure-characteristic secondary contrasts',
        'TableS8':'Exact-coordinate Guyot decomposition','TableS9':'Sources of reconstructed curve coordinates',
        'TableS10':'Standard-error and event-count accuracy','TableS11':'Accuracy by input condition'}
TABLES={}

def read(name):
    with (ROOT/name).open(newline='') as f:return list(csv.DictReader(f))

def number(value):
    if value=='':return None
    return float(value)

def mean(rows,key):return float(np.mean([float(r[key]) for r in rows]))

def display(value):
    if value is None:return '—'
    if isinstance(value,(int,np.integer)):return str(value)
    if isinstance(value,(float,np.floating)):
        if not math.isfinite(value):return 'unavailable'
        return f'{value:.6g}'
    return str(value)

def table(name,rows):
    TABLES[name]=rows
    with (OUT/(name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    def escape(v):
        s=display(v)
        for a,b in [('\\',r'\textbackslash{}'),('_',r'\_'),('%',r'\%'),('&',r'\&'),('#',r'\#')]:s=s.replace(a,b)
        return s.replace('—','--')
    lines=['% '+TITLES[name],r'\begin{tabular}{'+'l'*len(rows[0])+'}',r'\hline',
           ' & '.join(escape(k) for k in rows[0])+r' \\',r'\hline']
    lines += [' & '.join(escape(v) for v in r.values())+r' \\' for r in rows]
    lines += [r'\hline',r'\end{tabular}']
    (OUT/(name+'.tex')).write_text('\n'.join(lines)+'\n')

def ci(row,prefix='',scale=1):
    return {k:(scale*float(row[prefix+k]) if row[prefix+k] and math.isfinite(float(row[prefix+k])) else None)
            for k in ('estimate','lower','upper')}

def main():
    OUT.mkdir(parents=True,exist_ok=True)
    runs=read('data/run_manifest.csv');scores=read('build/scores.csv');pairs=read('build/primary_pairs.csv')
    current={(r['family'],r['metric'],r['contrast']):r for r in read('build/contrasts_current64.csv')}
    exploratory={(r['family'],r['metric'],r['contrast']):r for r in read('build/exploratory_contrasts.csv')}
    follow={(r['contrast'],r['metric']):r for r in read('build/followup_contrasts.csv') if r['scope']=='all'}
    primary=[r for r in scores if r['primary64']=='TRUE'];precision={r['run_id']:r for r in read('build/uncertainty.csv')}
    def n(experiment):return len({r['dataset_id'] for r in runs if r['experiment']==experiment})
    table('Table1',[dict(analysis=label,datasets=n(experiment)) for label,experiment in [
        ('File format','primary'),('Censoring marks','followup'),('Image resolution','followup'),('JPEG compression','followup'),
        ('Vector content','image_only'),('Plotting style','survminer'),('Automated Guyot workflow','guyot_png')]])
    selected={'taae_mean','absolute_rmst_difference_error','absolute_log_hr_error'}
    for name,family,contrast,left,right in [('Table2','primary','PDF_minus_PNG300','Native PDF','PNG300'),
                                           ('Table3','guyot','Guyot_minus_LLM_PNG300','Guyot PNG300','LLM PNG300')]:
        rows=[]
        for metric,label,unit in CURVE:
            if metric not in selected:continue
            row=current[(family,metric,contrast)]
            rows.append(dict(outcome=label,unit=unit,**{left:float(row['mean_left']),right:float(row['mean_right'])},**ci(row,'new_')))
        table(name,rows)
    table('TableS1',[dict(batch='first' if p['dataset'].startswith('PP32') else 'second',dataset=p['dataset'],
        PDF=float(p['taae_mean_PDF']),PNG300=float(p['taae_mean_PNG300']),difference=float(p['taae_difference']),
        later_PNG_recovery=next(r['assignment_status']=='recovery' for r in primary if r['dataset_id']==p['dataset'] and r['condition']=='PNG300')) for p in pairs])
    repeats=sorted(read('build/repeat_direct.csv'),key=lambda r:(r['dataset'],r['representation']))
    table('TableS2',[dict(dataset=r['dataset'],format=r['representation'],first_run=r['first_session'],repeat_run=r['repeat_session'],
                         direct_TAAE=float(r['taae_mean'])) for r in repeats])
    raster={(r['contrast'],r['metric']):r for r in read('build/image_only_bootstrap.csv') if r['scope']=='all'}
    rows=[]
    for metric,label,unit in CURVE:
        row=dict(outcome=label,unit=unit)
        for short,contrast in [('PDF_minus_imagePDF','PDF_minus_RASTER_PDF'),('imagePDF_minus_PNG300','RASTER_PDF_minus_PNG300')]:
            row.update({short+'_'+k:v for k,v in ci(raster[(contrast,metric)],scale=100 if unit=='pp' else 1).items()})
        rows.append(row)
    table('TableS3',rows)
    rows=[]
    for metric,label,unit in CURVE:
        g=current[('guyot',metric,'Guyot_minus_LLM_PNG300')];pdf=current[('primary',metric,'PDF_minus_PNG300')]
        rows.append(dict(outcome=label,unit=unit,Guyot=float(g['mean_left']),LLM_PNG300=float(g['mean_right']),LLM_PDF=float(pdf['mean_left']),**ci(g,'new_')))
    table('TableS4',rows)
    diagnostics=read('build/template_diagnostics.csv')
    table('TableS5',[dict(dataset=r['dataset'],n_per_arm=int(r['n_per_group']),censor_mark_adjustment=r['censor_glyph_handling'],
        original_TAAE_pp=float(r['original_TAAE_pp']),survminer_TAAE_pp=float(r['survminer_TAAE_pp']),
        matched_event_count_peak_error_pp=float(r['survminer_max_matched_height_difference_pp'])) for r in diagnostics])
    cohort={(r['analysis'],r['metric']):r for r in read('build/cohort_sensitivity.csv') if r['scope']=='all'}
    missing=json.loads((ROOT/'build/missing_bounds.json').read_text());sensitivity=[]
    for name,label,pop,pooled in [('first_batch','First batch',None,False),('historical62','Second before later recovery','historical62',False),
        ('historical62','Pooled before later recovery','historical62',True),('historical61','Second original assignments','historical61',False),
        ('historical61','Pooled original assignments','historical61',True),('second_completed','Second completed',None,False),
        ('current64','Current primary 64',None,True)]:
        if name=='current64':
            r=current[('primary','taae_mean','PDF_minus_PNG300')];v=ci(r,'new_');count=int(r['pairs']);assigned=len(pairs)
        else:
            r=next(r for r in read('build/'+name+'.csv') if r['metric']=='taae_mean') if pooled else cohort[(name,'taae_mean')]
            v=ci(r,scale=100);count=int(r['complete_pairs']);assigned=len(pairs) if pooled else len(pairs)//2
        low=high=obslo=obshi=None
        if pop:
            pp=read('build/'+pop+'_pairs.csv');pp=pp if pooled else [r for r in pp if r['dataset'].startswith('RP32')]
            low=mean(pp,'lower_bound')*100;high=mean(pp,'upper_bound')*100
            if pop=='historical62':
                block=missing['pooled62' if pooled else 'second_main']
                obslo=block['fill_observed_max'] if pooled else block['fill_observed_max']['difference']
                obshi=block['fill_observed_min'] if pooled else block['fill_observed_min']['difference']
        sensitivity.append(dict(analysis=label,complete_pairs=count,assigned_pairs=assigned,**v,worst_lower=low,worst_upper=high,observed_range_lower=obslo,observed_range_upper=obshi))
    table('TableS6',sensitivity)
    follow_metrics=CURVE[1:]+[('abs_rel_se_log_hr_error','SE(log HR), absolute relative error','fraction'),
                             ('signed_events_total_error_A','Events, A, signed error','count'),('signed_events_total_error_B','Events, B, signed error','count')]
    rows=[]
    for metric,label,unit in follow_metrics:
        row=dict(outcome=label,unit=unit)
        for contrast in ('marks_vector','marks_raster','resolution_150','resolution_100','compression_100'):
            row.update({contrast+'_'+k:v for k,v in ci(follow[(contrast,metric)],scale=100 if unit=='pp' else 1).items()})
        rows.append(row)
    table('TableS7',rows)
    rows=[]
    decomposition=[c for c in CURVE if c[0] in selected|{'maximum_A','maximum_B'}]
    decomposition += [('abs_rel_se_log_hr_error','SE(log HR), absolute relative error','fraction'),('rel_se_log_hr_error','SE(log HR), signed relative error','fraction'),
                      ('abs_events_total_error_A','Events, A, absolute error','count'),('signed_events_total_error_A','Events, A, signed error','count')]
    for metric,label,unit in decomposition:
        family='decomposition' if metric in {m for m,_,_ in CURVE} else 'precision'
        scale=100 if unit=='pp' else 1
        row=dict(outcome=label,unit=unit)
        for method in ('PDF','GuyotExact','PNG300','Guyot'):
            selected_runs=[r for r in scores if (r['primary64']=='TRUE' and r['condition']==method) or r['method']==method]
            row[method]=float(np.mean([float(r[metric]) if family=='decomposition' else float(precision[r['run_id']][metric]) for r in selected_runs]))*scale
        for contrast in ('Guyot_minus_GuyotExact','GuyotExact_minus_PDF'):
            row.update({contrast+'_'+k:v for k,v in ci(exploratory[(family,metric,contrast)],'new_').items()})
        rows.append(row)
    table('TableS8',rows)
    counts=read('build/extraction_strategy_counts.csv');rows=[]
    for exp,cond,label in [('primary','PDF','Primary native PDF'),('primary','PNG300','Primary PNG300'),('image_only','RASTER_PDF','Image-only PDF'),('survminer','PDF','Survminer native PDF'),('survminer','PNG300','Survminer PNG300')]:
        r=next(r for r in counts if r['experiment']==exp and r['condition']==cond)
        rows.append(dict(experiment=label,**{k:int(r[k]) for k in ('registered','complete','path_parsing','pixel_segmentation','hybrid')}))
    rows.append(dict(experiment='Total',**{k:sum(r[k] for r in rows) for k in rows[0] if k!='experiment'}));table('TableS9',rows)
    rows=[]
    for metric,label in SE:
        row=dict(outcome=label)
        for method in ('PDF','PNG300','Guyot'):
            rr=[r for r in scores if (r['primary64']=='TRUE' and r['condition']==method) or r['method']==method]
            row[method]=float(np.mean([float(precision[r['run_id']][metric]) for r in rr]))
        for contrast in ('PDF_minus_PNG300','Guyot_minus_PNG300'):
            row.update({contrast+'_'+k:v for k,v in ci(exploratory[('precision',metric,contrast)],'new_').items()})
        rows.append(row)
    table('TableS10',rows)
    conditions={r['condition']:r for r in read('build/followup_conditions.csv') if r['scope']=='all'}
    rows=[]
    for condition in ('PDF','NOMARK_PDF','PNG300','NOMARK_PNG300','PNG150','PNG100','JPEG100'):
        r=conditions[condition]
        rows.append(dict(condition=condition,**{k+'_pp':float(r['taae_'+k])*100 for k in ('median','q1','q3','maximum')},
            maximum_pointwise_pp=float(r['maximum_curve_error'])*100,mean_log_HR_error=float(r['mean_absolute_log_hr_error']),
            SE_within_1pct=int(r['se_log_hr_within_1pct']),exact_events=int(r['exact_event_totals_both_arms']),datasets=int(r['datasets'])))
    table('TableS11',rows)
    # Unrounded values underlying prose claims, derived from the same populations.
    prose=dict(primary_datasets=len(pairs),pdf_lower_taae_pairs=sum(float(r['taae_difference'])<0 for r in pairs),
        missing_bounds=missing,revision_sensitivity=json.loads((ROOT/'build/revision_sensitivity.json').read_text()),
        method_audit=json.loads((ROOT/'build/method_audit.json').read_text()),
        event_and_se_accuracy={})
    for method in ('PDF','PNG300','Guyot'):
        rr=[precision[r['run_id']] for r in scores if (r['primary64']=='TRUE' and r['condition']==method) or r['method']==method]
        prose['event_and_se_accuracy'][method]=dict(datasets=len(rr),
            both_event_totals_exact=sum(float(r['signed_events_total_error_A'])==float(r['signed_events_total_error_B'])==0 for r in rr),
            mean_abs_relative_cox_se=mean(rr,'abs_rel_se_log_hr_error'))
    # Make the descriptive results quoted in the prose directly discoverable.
    # The primary IQR uses the original exclusive/type-6 convention; these are
    # fresh summaries of individual scores, never copied manuscript numbers.
    prose['primary_taae_pp']={}
    for method in ('PDF','PNG300'):
        values=np.array([float(r['taae_mean_'+method])*100 for r in pairs])
        q1,median,q3=np.quantile(values,[.25,.5,.75],method='weibull')
        prose['primary_taae_pp'][method]=dict(datasets=len(values),mean=float(values.mean()),
            median=float(median),q1=float(q1),q3=float(q3),minimum=float(values.min()),maximum=float(values.max()))
    prose['repeatability']={method:dict(datasets=len([r for r in repeats if r['representation']==method]),
        mean_direct_taae_pp=100*mean([r for r in repeats if r['representation']==method],'taae_mean'))
        for method in ('PDF','PNG300')}
    image_only=read('build/image_only_pairs.csv')
    differences=[float(r['taae_mean_RASTER_PDF'])-float(r['taae_mean_PNG300']) for r in image_only]
    prose['image_only_pdf']=dict(datasets=len(image_only),
        mean_taae_pp={m:100*mean(image_only,'taae_mean_'+m) for m in ('PDF','RASTER_PDF','PNG300')},
        lower_than_png300=sum(v<0 for v in differences),higher_than_png300=sum(v>0 for v in differences))
    plotting_style=read('build/survminer_pairs.csv')
    prose['plotting_style']=dict(datasets=len(plotting_style),
        mean_taae_pp={m:100*mean(plotting_style,'taae_mean_new_'+m) for m in ('PDF','PNG300')},
        pdf_lower_pairs=sum(float(r['taae_mean_new_PDF'])<float(r['taae_mean_new_PNG300']) for r in plotting_style))
    prose['figure_characteristic_conditions']=TABLES['TableS11']
    (ROOT/'build/reported_statistics.json').write_text(json.dumps(prose,indent=2,allow_nan=False)+'\n')
    markup=['<!doctype html><meta charset="utf-8"><title>Reproduced manuscript tables</title><style>body{font:14px system-ui;margin:2rem}table{border-collapse:collapse}th,td{padding:6px;border-bottom:1px solid #ddd;text-align:right}th:first-child,td:first-child{text-align:left}.scroll{overflow:auto}h2{margin-top:3rem}</style><h1>Reproduced manuscript tables</h1><p>Calculated from frozen IPD. CSV files retain full numerical precision. pp = percentage points.</p>']
    for name,rows in TABLES.items():
        markup+=['<h2>'+html.escape(name+': '+TITLES[name])+'</h2><div class="scroll"><table><tr>']
        markup+=['<th>'+html.escape(k)+'</th>' for k in rows[0]];markup+=['</tr>']
        for r in rows:markup+=['<tr>']+['<td>'+html.escape(display(v))+'</td>' for v in r.values()]+['</tr>']
        markup+=['</table></div>']
    (OUT/'tables.html').write_text(''.join(markup))
    print(f'Exported {len(TABLES)} manuscript tables and supporting prose statistics.')

if __name__=='__main__':main()
