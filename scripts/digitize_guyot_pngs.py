#!/usr/bin/env python3
"""Optional replay of the original automated PNG digitizer; no model/API calls.

Inputs: the actual 64 primary PNGs, public sample sizes, data/guyot_workflow.zip,
and a separately installed frozen OCR runtime (node/bin/node, node/node_modules,
language/eng.traineddata). Outputs: build/guyot_digitization/<dataset>.json.
Run: python scripts/digitize_guyot_pngs.py --runtime /path/to/frozen/runtime
Smoke test: append --dataset PP32_S1_D01. The default manuscript reproduction
uses the saved digitization/IPD and does not require OCR. Runtime files must
match the original lock. No downloads or package installation occur here.
The archived digitizer receives ONLY the PNG path and public sample sizes.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import importlib.util
import json
import os
import subprocess
import sys
import zipfile

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--runtime',type=Path,required=True)
    p.add_argument('--dataset')
    args=p.parse_args();runtime=args.runtime.resolve()
    target=ROOT/'build/guyot_digitization';target.mkdir(parents=True,exist_ok=True)
    source=ROOT/'build/guyot_source';source.mkdir(exist_ok=True)
    members=['src/km_analysis/deterministic_extraction.py','src/runner/geometry.py','scripts/deterministic/proposal_a_ocr.cjs','scripts/main_study/extract_comparator.py']
    with zipfile.ZipFile(ROOT/'data/guyot_workflow.zip') as archive:
        lock=json.loads(archive.read('environment/comparator-ocr-lock.json'))
        for member in members:
            destination=source/member;destination.parent.mkdir(parents=True,exist_ok=True)
            destination.write_bytes(archive.read(member))
    for relative,expected in lock['files_sha256'].items():
        data=(runtime/relative).read_bytes()
        if hashlib.sha256(data).hexdigest()!=expected:raise ValueError('OCR runtime differs from study lock: '+relative)
    version=subprocess.check_output([str(runtime/'node/bin/node'),'--version'],text=True).strip()
    if version!=lock['node_version']:raise ValueError('Wrong Node version: '+version)
    os.environ['CODEX_PRIMARY_DEPS']=str(runtime)
    os.environ['KM_OCR_LANGUAGE_DIR']=str(runtime/'language')
    sys.path.insert(0,str(source/'src'))
    from km_analysis.deterministic_extraction import OCR,extract
    spec=importlib.util.spec_from_file_location('original_digitizer_validation',source/'scripts/main_study/extract_comparator.py')
    validator=importlib.util.module_from_spec(spec);spec.loader.exec_module(validator)
    with (ROOT/'data/run_manifest.csv').open() as stream:
        runs=[r for r in csv.DictReader(stream) if r['method']=='Guyot']
    if args.dataset:
        runs=[r for r in runs if r['dataset_id']==args.dataset]
        if not runs:raise ValueError('Unknown primary dataset')
    # Original workflow starts a fresh OCR worker for each image; preserve that
    # boundary even though adaptive OCR learning is also explicitly disabled.
    for r in runs:
        ocr=OCR()
        try:
            caption=dict(A=int(r['n_per_arm']),B=int(r['n_per_arm']))
            result=extract(ROOT/r['input_figure'],caption,ocr)
            validator.validate(result,caption)
        finally:ocr.close()
        expected=json.loads((ROOT/'data/guyot_digitization'/(r['dataset_id']+'.json')).read_text())
        # Isolation receipts are host-specific. Compare the actual numerical
        # coordinates/risk/sample-size payload received by the Guyot algorithm.
        if result['usable']!=expected['usable']:raise ValueError('Digitizer usability differs: '+r['dataset_id'])
        for g in ('A','B'):
            for key in ('coordinates','risk','n'):
                if result['arms'][g][key]!=expected['arms'][g][key]:raise ValueError(f'Digitization discrepancy: {r["dataset_id"]} {g} {key}')
        (target/(r['dataset_id']+'.json')).write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
        print(r['dataset_id']+': unchanged Guyot input payload',flush=True)

if __name__=='__main__':main()
