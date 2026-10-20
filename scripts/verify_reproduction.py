#!/usr/bin/env python3
"""Verify public inputs, study populations, and regenerated manuscript results.

Inputs: data/source_checksums.csv, run_manifest.csv, frozen_results/validation/,
numeric display references, and generated build/ files. Outputs: build/verification.json.
Run: python scripts/verify_reproduction.py --inputs
     python scripts/verify_reproduction.py --results
Validation outputs are read ONLY here, never used to compute analysis results.
Numerical tolerances allow floating-point serialization, not changed estimates:
1e-10 for R/display-scale contrasts and 1e-9 for independently fitted Python Cox.
Published rounded cells must agree at their stated display precision.
Path-redacted administrative artifacts carry separate original fingerprints;
reference IPD, reconstruction outputs, assigned figures and model programs do
not permit this transformation. Historical hashes are not rewritten.
"""
from pathlib import Path
from collections import Counter
import argparse
import csv
import hashlib
import json
import math
import re
import zipfile

ROOT=Path(__file__).resolve().parents[1]
REPORT={}

def read(path):
    with (ROOT/path).open(newline='') as f:return list(csv.DictReader(f))

def check(condition,message):
    if not condition:raise AssertionError(message)

def numeric(a,b,label,tolerance=1e-10):
    if a in ('',None) or b in ('',None):
        check(a in ('',None) and b in ('',None),label+': missing-value disagreement');return 0
    a,b=float(a),float(b)
    check(math.isfinite(a)==math.isfinite(b),label+': nonfinite disagreement')
    error=abs(a-b) if math.isfinite(a) else 0
    check(error<=tolerance,f'{label}: {a:.17g} versus {b:.17g}, difference {error:.3g}')
    return error

def inputs():
    ledger=read('data/source_checksums.csv');archives={};checked={};redacted=set()
    try:
        for r in ledger:
            key=(r['path'],r['archive_member'])
            if key not in checked:
                if r['archive_member']:
                    if r['path'] not in archives:archives[r['path']]=zipfile.ZipFile(ROOT/r['path'])
                    data=archives[r['path']].read(r['archive_member'])
                else:data=(ROOT/r['path']).read_bytes()
                checked[key]=(hashlib.sha256(data).hexdigest(),len(data))
            check(checked[key]==(r['sha256'],int(r['bytes'])),'Input changed: '+str(key))
            transformation=r.get('transformation','')
            if transformation:
                metadata=(transformation=='private-path-redaction-v1' and
                    r['path'] in ('frozen_results/reconstruction_programs.zip','data/guyot_workflow.zip') and
                    r['archive_member'].endswith('.json'))
                generator=(transformation=='portable-rscript-path-v1' and r['path']=='data/input_generation.zip' and
                    r['archive_member']=='scripts/main_study/prepare_pdf_png300_survminer_template.py')
                check(metadata or generator,'Unapproved scientific-artifact transformation: '+str(key))
                check(bool(re.fullmatch('[0-9a-f]{64}',r.get('original_sha256',''))) and
                      r['original_sha256']!=r['sha256'] and int(r['original_bytes'])>0,
                      'Missing original fingerprint for public derivative: '+str(key))
                redacted.add(key)
            else:
                check(not r.get('original_sha256') and not r.get('original_bytes'),
                      'Original fingerprint requires a documented transformation: '+str(key))
    finally:
        for z in archives.values():z.close()
    runs=read('data/run_manifest.csv')
    check(len(runs)==550 and len({r['run_id'] for r in runs})==550,'Registration inventory changed')
    complete=[r for r in runs if r['complete']=='True']
    check(len(complete)==543,'Completed output count changed')
    with zipfile.ZipFile(ROOT/'frozen_results/reconstruction_programs.zip') as archive:
        members=set(archive.namelist())
        for r in runs:
            for field in ('public_manifest_member','prompt_member','replay_receipt_member'):
                check(not r[field] or r[field] in members,'Missing archived '+field+': '+r['run_id'])
    for r in runs:
        check(bool(r['reconstruction_output'])==(r['complete']=='True'),'Missing-output flag: '+r['run_id'])
        if r['reconstruction_output']:check((ROOT/r['reconstruction_output']).is_file(),'Missing IPD: '+r['run_id'])
        for name in ('reference_ipd','input_figure'):check((ROOT/r[name]).is_file(),'Missing '+name+': '+r['run_id'])
        if r['replaces_run_id']:
            old=next(x for x in runs if x['run_id']==r['replaces_run_id'])
            check(old['complete']=='False' and old['dataset_id']==r['dataset_id'] and old['condition']==r['condition'],'Recovery mapping changed')
    for field,n in [('primary64',64),('complete_case62',62),('original61',61)]:
        population=[r for r in runs if r[field]=='True']
        check(len(population)==2*n and len({r['dataset_id'] for r in population})==n,field+' population changed')
        check(all(v==2 for v in Counter(r['dataset_id'] for r in population).values()),field+' pairing changed')
    repeated=[r for r in runs if r['repeatability']=='True']
    check(len(repeated)==48 and len({r['dataset_id'] for r in repeated})==24,'Repeatability population changed')
    check(Counter(r['scenario'] for r in runs if r['primary64']=='True')==dict.fromkeys(('S1','S2','S3','S4'),32),'Primary scenario balance changed')
    REPORT['inputs']=dict(checksum_records=len(ledger),distinct_payloads=len(checked),
        path_redacted_payloads=len(redacted),registered_runs=len(runs),
        complete_outputs=len(complete),missing_outputs=len(runs)-len(complete))


def compare_rows(actual,expected,keys,columns,tolerance=1e-10,families=None):
    aa=read(actual);ee=read(expected)
    if families is not None:ee=[r for r in ee if r['family'] in families]
    by={tuple(r[k] for k in keys):r for r in ee};worst=0;count=0
    actual_keys=[tuple(r[k] for k in keys) for r in aa]
    check(len(by)==len(ee),expected+': duplicate validation keys')
    check(len(set(actual_keys))==len(aa),actual+': duplicate result keys')
    check(set(actual_keys)==set(by),actual+': missing or extra comparison rows')
    for r in aa:
        key=tuple(r[k] for k in keys);check(key in by,'Unmatched validation key: '+str(key))
        for column in columns:
            worst=max(worst,numeric(r[column],by[key][column],actual+':'+str(key)+':'+column,tolerance));count+=1
    return dict(rows=len(aa),values=count,max_absolute_difference=worst)


def recursive_numbers(a,b,label):
    if isinstance(a,dict):
        for k in a:recursive_numbers(a[k],b[k],label+'/'+k)
    elif isinstance(a,list):
        check(len(a)==len(b),label+': list length changed')
        for i,(x,y) in enumerate(zip(a,b)):recursive_numbers(x,y,label+'/'+str(i))
    elif isinstance(a,(int,float)) and not isinstance(a,bool):numeric(a,b,label)


def published_numbers(cell):
    # Strip TeX decoration while retaining the decimal precision of every token.
    s=re.sub(r'([+-]?\d+(?:\.\d+)?)\\times10\^\{([+-]?\d+)\}',r'\1e\2',cell)
    s=s.replace('$-$','-').replace('$','').replace(r'\,','').replace(r'\;','')
    return re.findall(r'(?<![A-Za-z0-9])[-+]?\d+(?:\.\d+)?(?:e[-+]?\d+)?',s)

def check_display(value,token,label):
    mantissa,*exponent=token.lower().split('e');decimals=len(mantissa.split('.')[1]) if '.' in mantissa else 0
    quantum=10.0**((int(exponent[0]) if exponent else 0)-decimals)
    tolerance=quantum/2+1e-12*max(abs(float(token)),1e-12)
    numeric(value,float(token),label,tolerance)

def latex_rows(text):
    return [line for line in text.splitlines() if '&' in line and line.rstrip().endswith(r'\\')]

def display_tokens(rows,columns,expected,label):
    check(len(rows)==len(expected),label+': displayed row count changed')
    count=0
    for i,(r,tokens) in enumerate(zip(rows,expected)):
        values=[r[c] for c in columns]
        check(len(tokens)==len(values),f'{label} row {i+1}: unexpected display columns {tokens}')
        for value,token in zip(values,tokens):check_display(value,token,f'{label} row {i+1}');count+=1
    return count

def display_cells(rows,columns,latex,label):
    expected=[sum([published_numbers(c) for c in line.split('&')[1:]],[]) for line in latex]
    return display_tokens(rows,columns,expected,label)


def published_tables():
    validation=ROOT/'frozen_results/validation/manuscript_tables'
    # Retain only numeric tokens and their display precision, not private prose.
    inline=json.loads((ROOT/'frozen_results/validation/inline_table_cells.json').read_text())
    count=0
    for name in ('Table2','Table3'):
        rows=read('build/tables/'+name+'.csv')
        count+=display_tokens(rows,list(rows[0])[2:],inline[name],name)
    s1=read('build/tables/TableS1.csv');s2=read('build/tables/TableS2.csv')
    count+=display_cells(s1[:32],['PDF','PNG300','difference'],latex_rows((validation/'dataset_pairs.tex').read_text()),'S1 first')
    count+=display_tokens(s1[32:],['PDF','PNG300','difference'],inline['TableS1_second'],'S1 second')
    # The format column contains no numerical values.
    count+=display_cells(s2[:16],['direct_TAAE'],latex_rows((validation/'repeats.tex').read_text()),'S2 first')
    count+=display_tokens(s2[16:],['direct_TAAE'],inline['TableS2_second'],'S2 second')
    specs=[('TableS3','raster_secondary.tex',2),('TableS4','guyot_secondary.tex',2),('TableS8','decomposition.tex',2),('TableS10','se_estimands_pooled64.tex',1)]
    for table,file,skip in specs:
        rows=read('build/tables/'+table+'.csv')
        count+=display_cells(rows,list(rows[0])[skip:],latex_rows((validation/file).read_text()),table)
    # S7 has an undefined CI for a degenerate contrast: compare its single
    # estimate and preserve the blank CI, rather than fabricating endpoints.
    rows=read('build/tables/TableS7.csv');tex=latex_rows((validation/'followup_secondary.tex').read_text())
    check(len(rows)==len(tex),'S7 row count')
    for i,(r,line) in enumerate(zip(rows,tex)):
        values=[v for v in list(r.values())[2:] if v!=''];tokens=sum([published_numbers(c) for c in line.split('&')[1:]],[])
        check(len(values)==len(tokens),'S7 parse')
        for v,t in zip(values,tokens):check_display(v,t,'S7 row '+str(i));count+=1
    rows=read('build/tables/TableS9.csv')
    count+=display_cells(rows,list(rows[0])[1:],latex_rows((validation/'extraction_strategies.tex').read_text()),'S9')
    rows=read('build/tables/TableS5.csv');tex=latex_rows((validation/'extraction_template_diagnostics.tex').read_text())
    check(len(rows)==len(tex),'S5 row count')
    for r,line in zip(rows,tex):
        cells=line.split('&');values=[r['n_per_arm'],r['original_TAAE_pp'],r['survminer_TAAE_pp'],r['matched_event_count_peak_error_pp']]
        tokens=sum([published_numbers(cells[i]) for i in (1,3,4,5)],[])
        for v,t in zip(values,tokens):check_display(v,t,'S5 '+r['dataset']);count+=1
    rows=read('build/tables/TableS11.csv');tex=latex_rows((validation/'followup_conditions.tex').read_text())
    check(len(rows)==len(tex),'S11 row count')
    for r,line in zip(rows,tex):
        values=[r[k] for k in ('median_pp','q1_pp','q3_pp','maximum_pp','maximum_pointwise_pp','mean_log_HR_error','SE_within_1pct','datasets','exact_events','datasets')]
        tokens=sum([published_numbers(c) for c in line.split('&')[1:]],[])
        check(len(values)==len(tokens),'S11 parse')
        for v,t in zip(values,tokens):check_display(v,t,'S11 '+r['condition']);count+=1
    rows=read('build/tables/TableS6.csv');expected=inline['TableS6']
    check(len(expected)==len(rows),'S6 rows')
    for r,tokens in zip(rows,expected):
        values=[v for k,v in r.items() if k!='analysis' and v!='']
        check(len(values)==len(tokens),'S6 parse')
        for v,t in zip(values,tokens):check_display(v,t,'S6 '+r['analysis']);count+=1
    return count


def results():
    scored_rows=read('build/scores.csv')
    scores={r['run_id']:r for r in scored_rows}
    expected_ids={r['run_id'] for r in read('data/run_manifest.csv') if r['complete']=='True' and r['experiment']!='effort_screen'}
    check(len(scored_rows)==len(scores)==528 and set(scores)==expected_ids,'Missing, extra or duplicate scored runs')
    metrics=[k for k in next(iter(scores.values())) if k.startswith(('taae','maximum','duration_above','signed_rmst','absolute_rmst','absolute_log_hr'))]
    count=0;worst=0;seen=set()
    for p in (ROOT/'frozen_results/validation').glob('*scores.csv'):
        for r in read(p):
            sid=r.get('session_id') or ('GuyotExact_' if p.name=='exact_scores.csv' else 'Guyot_')+r['dataset']
            if sid not in scores:continue
            seen.add(sid)
            for m in metrics:
                if m in r:worst=max(worst,numeric(scores[sid][m],r[m],sid+':'+m));count+=1
    check(len(seen)==528,'Not every scored reconstruction has an independent validation authority')
    REPORT['scores']=dict(runs=len(seen),values=count,max_absolute_difference=worst)
    for actual,expected,keys in [('contrasts_current64.csv','four_scenario_contrasts.csv',['family','metric','contrast']),
        ('historical62.csv','historical62.csv',['metric']),('image_only_bootstrap.csv','image_only_bootstrap.csv',['contrast','metric','scope']),
        ('survminer_bootstrap.csv','survminer_bootstrap.csv',['contrast','metric','scope']),('followup_contrasts.csv','followup_contrasts.csv',['contrast','metric','scope'])]:
        columns=['new_'+k for k in ('estimate','lower','upper','standard_error')] if actual.startswith('contrasts') else ['estimate','lower','upper','standard_error']
        REPORT[actual]=compare_rows('build/'+actual,'frozen_results/validation/'+expected,keys,columns,
            families={'primary','guyot'} if actual.startswith('contrasts') else None)
    current={tuple(r[k] for k in ('family','metric','contrast')):r for r in read('frozen_results/validation/four_scenario_contrasts.csv')}
    corrected={(r['metric'],r['contrast']):r for r in read('frozen_results/validation/exact_peak_contrasts.csv')}
    worst=0;count=0
    for r in read('build/exploratory_contrasts.csv'):
        old=current[(r['family'],r['metric'],r['contrast'])]
        for k in ('estimate','lower','upper'):
            target=corrected[(r['metric'],r['contrast'])]['csv_'+k] if r['family']=='decomposition' and r['metric'].startswith('maximum') else old['new_'+k]
            worst=max(worst,numeric(r['new_'+k],target,'Exploratory '+r['metric']+' '+r['contrast']));count+=1
    REPORT['exploratory']=dict(values=count,max_absolute_difference=worst)
    uncertainty_rows=read('build/uncertainty.csv')
    uncertainty={r['run_id']:r for r in uncertainty_rows};count=0;worst=0
    check(len(uncertainty_rows)==len(uncertainty)==528 and set(uncertainty)==expected_ids,'Uncertainty population changed')
    for r in read('frozen_results/validation/se_estimands_runs.csv'):
        sid=r['session_id'].replace('GUYOTEXACT_','GuyotExact_').replace('GUYOT_','Guyot_')
        if sid not in uncertainty:continue
        for k,v in r.items():
            if k not in uncertainty[sid]:continue
            try:float(v)
            except ValueError:continue
            worst=max(worst,numeric(uncertainty[sid][k],v,sid+':'+k,1e-9));count+=1
    REPORT['uncertainty']=dict(values=count,max_absolute_difference=worst)
    direct={(r['dataset'],r['representation']):r for r in read('build/repeat_direct.csv')};count=0
    for path in ('original_repeat.csv','second_repeat.csv'):
        for r in read('frozen_results/validation/'+path):
            for m in metrics:numeric(direct[(r['dataset'],r['representation'])][m],r[m],'Repeat '+r['dataset']+':'+m);count+=1
    REPORT['direct_repeat_values']=count
    for name in ('missing_bounds.json','revision_sensitivity.json'):
        recursive_numbers(json.loads((ROOT/'frozen_results/validation'/name).read_text()),json.loads((ROOT/'build'/name).read_text()),name)
    for name,expected,keys in [('template_diagnostics.csv','survminer_template_comparison.csv',['dataset']),('template_arm_diagnostics.csv','survminer_arm_diagnostics.csv',['dataset','template','group'])]:
        aa=read('build/'+name);ee=read('frozen_results/validation/'+expected);by={tuple(r[k] for k in keys):r for r in ee}
        for r in aa:
            for k,v in by[tuple(r[k] for k in keys)].items():
                if k not in r:continue
                try:float(v)
                except ValueError:continue
                numeric(r[k],v,name+':'+k)
    REPORT['published_numeric_cells']=published_tables()
    figs=list((ROOT/'build/figures').glob('*.pdf'));check(len(figs)==6,'Expected six manuscript figures')
    for p in figs:check(p.read_bytes().startswith(b'%PDF-'),'Invalid PDF '+p.name)
    REPORT['figures']=[p.name for p in sorted(figs)]
    check(len(list((ROOT/'build/tables').glob('Table*.csv')))==14,'Expected fourteen manuscript tables')


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--inputs',action='store_true');parser.add_argument('--results',action='store_true')
    args=parser.parse_args()
    if not (args.inputs or args.results):parser.error('Choose --inputs or --results')
    path=ROOT/'build/verification.json'
    if path.exists():REPORT.update(json.loads(path.read_text()))
    if args.inputs:inputs()
    if args.results:results()
    (ROOT/'build').mkdir(exist_ok=True);path.write_text(json.dumps(REPORT,indent=2,allow_nan=False)+'\n')
    print('PASS:', 'input integrity' if args.inputs else '', 'numerical and manuscript validation' if args.results else '')

if __name__=='__main__':main()
