"""Experiment 2B: serial, resumable twenty-dataset screen and paired order checks."""
from __future__ import annotations
import os
for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import shutil
import time
import numpy as np
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from fit_adapter import fit_model, preprocessing_scale, native_details

ROOT = RUN
POLICY = json.loads((ROOT/'protocol.json').read_text(encoding='utf-8'))
ROSTER = POLICY['datasets']
METHODS = ('manhattan', 'euclidean')
RESULTS = ROOT/'results'
GRID = []

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def read(path, default=None):
    return json.loads(Path(path).read_text(encoding='utf-8')) if Path(path).exists() else default

def lines(path):
    return [json.loads(x) for x in Path(path).read_text(encoding='utf-8').splitlines() if x.strip()] if Path(path).exists() else []

def save(path, value):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    temporary=path.with_suffix(path.suffix+'.tmp')
    temporary.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n',encoding='utf-8')
    for attempt in range(40):
        try:
            temporary.replace(path); return
        except PermissionError:
            if attempt==39: raise
            time.sleep(.05)

def append(path,row):
    with Path(path).open('a',encoding='utf-8') as stream:
        stream.write(json.dumps(row,allow_nan=False)+'\n')

def progress(**values):
    values.update(pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                  grid_fits=len(GRID),permutation_fits=len(lines(RESULTS/'permutation.jsonl')))
    save(RESULTS/'progress.json',values); print(json.dumps(values),flush=True)
    if (ROOT/'STOP').exists(): raise SystemExit('Paused at checkpoint. Remove STOP to resume.')

def load(dataset):
    path=ROOT/'data'/f'{dataset}.npz'
    with np.load(path,allow_pickle=False) as arrays:
        X,y=arrays['X'],arrays['y']
    assert X.dtype==np.float64 and y.shape==(len(X),) and np.isfinite(X).all()
    return X,y,path

def score(y,labels):
    return dict(ari=float(adjusted_rand_score(y,labels)),ami=float(adjusted_mutual_info_score(y,labels)))

def best(rows):
    return min(rows,key=lambda row:(-row['ari'],row['min_pts'],row['radius_index']))

def plan(dataset,X,path):
    target=RESULTS/'plans'/f'{dataset}.json'
    old=read(target)
    if old:
        assert old['input_sha256']==sha(path); return old
    design=dict(dataset=dataset,kind='real',input_sha256=sha(path),radius_count=32,
                min_pts=[0,5,10,20,40,80],merge_scale=1.4,methods={})
    for method in METHODS:
        values=pdist(X,metric='cityblock' if method=='manhattan' else 'euclidean')
        dense=squareform(values)
        near=float(np.median(np.partition(dense,5,axis=1)[:,5]))
        positive=values[values>0]
        far=float(np.median(positive))
        low,high=.35*near,.9*far
        if not (0<low<high and np.isfinite(high)):
            raise ValueError(f'Invalid predeclared radius interval for {dataset}/{method}: {low}, {high}')
        scale=preprocessing_scale(method,X)
        assert scale>0
        design['methods'][method]=dict(median_5nn=near,median_pair=far,physical_low=low,physical_high=high,
                                      internal_scale=scale,radii=(np.geomspace(low,high,32)/scale).tolist())
    save(target,design); return design


def pilot():
    target=RESULTS/'readiness.json'
    if target.exists(): return
    records=[]
    for dataset in ROSTER:
        X,y,path=load(dataset)
        design=plan(dataset,X,path)
        if dataset in POLICY['anchors']: continue
        for method in METHODS:
            progress(stage='unscored_readiness',dataset=dataset,method=method)
            radius=design['methods'][method]['radii'][16]
            fitted=fit_model(method,X,radius,10,1.4)
            labels=fitted.pop('labels')
            assert labels.shape==y.shape and np.isfinite(labels).all()
            records.append(dict(dataset=dataset,method=method,radius=radius,min_pts=10,**fitted))
    save(target,dict(status='passed',records=records,note='Full-input execution/calibration pilot; no ARI/AMI used.'))

def cases(dataset,X,y,design,indices,minima):
    seen={(r['dataset'],r['method'],r['radius_index'],r['min_pts']) for r in GRID}
    for index in indices:
        for minimum in minima:
            for method in METHODS:
                if (dataset,method,index,minimum) in seen: continue
                geometry=design['methods'][method]
                radius=geometry['physical_low']*(geometry['physical_high']/geometry['physical_low'])**(index/31)/geometry['internal_scale']
                if len(GRID)%20==0: progress(stage='grid',dataset=dataset,radius_index=index,method=method)
                result=fit_model(method,X,radius,minimum,1.4)
                labels=result.pop('labels')
                row=dict(dataset=dataset,kind='real',method=method,radius=float(radius),radius_index=index,
                         min_pts=minimum,merge_scale=1.4,status='ok',origin='fresh',**result,**score(y,labels))
                append(RESULTS/'grid.jsonl',row); GRID.append(row)

def select(dataset,X,y):
    path=RESULTS/'selected.json'; selected=read(path,[])
    for method in METHODS:
        winner=best([r for r in GRID if r['dataset']==dataset and r['method']==method])
        old=next((r for r in selected if r['dataset']==dataset and r['method']==method),None)
        if old:
            assert (old['radius_index'],old['min_pts'])==(winner['radius_index'],winner['min_pts'])
            continue
        progress(stage='selected_repetitions',dataset=dataset,method=method)
        warm=fit_model(method,X,winner['radius'],winner['min_pts'],1.4)
        labels=warm['labels']; measured=score(y,labels)
        assert abs(measured['ari']-winner['ari'])<1e-12 and abs(measured['ami']-winner['ami'])<1e-12
        timings=[]
        for repeat in range(5):
            result=fit_model(method,X,winner['radius'],winner['min_pts'],1.4)
            assert adjusted_rand_score(labels,result['labels'])==1
            timings.append(result['fit_seconds'])
        target=RESULTS/'assignments'/f'{dataset}_{method}_selected.npz'
        target.parent.mkdir(parents=True,exist_ok=True); np.savez_compressed(target,labels=labels)
        selected.append(dict(winner,timing_seconds=timings,fit_seconds_median=float(np.median(timings)),
                             labels_path=str(target.relative_to(ROOT))))
        save(path,selected)

def grid(dataset):
    X,y,path=load(dataset); design=plan(dataset,X,path)
    cases(dataset,X,y,design,range(32),design['min_pts'])
    target=RESULTS/'plans'/f'{dataset}_extension.json'; extension=read(target)
    if extension is None:
        winners={method:best([r for r in GRID if r['dataset']==dataset and r['method']==method]) for method in METHODS}
        radius_edge=any(r['radius_index'] in (0,31) for r in winners.values())
        min_edge=any(r['min_pts']==80 for r in winners.values())
        extension=dict(radius_edge=radius_edge,min_edge=min_edge,
                       indices=list(range(-4,36)) if radius_edge else list(range(32)),
                       min_pts=design['min_pts']+([120,160] if min_edge else []),initial_winners=winners)
        save(target,extension)
    cases(dataset,X,y,design,extension['indices'],extension['min_pts'])
    select(dataset,X,y)
    progress(stage='dataset_complete',dataset=dataset)

def permutations():
    chosen={(r['dataset'],r['method']):r for r in read(RESULTS/'selected.json')}
    seen={(r['dataset'],r['method'],r['replicate']) for r in lines(RESULTS/'permutation.jsonl')}
    for dataset in ROSTER:
        X,y,path=load(dataset)
        for replicate in range(5):
            order=np.random.default_rng(np.random.SeedSequence([101,ROSTER.index(dataset),replicate,867])).permutation(len(X))
            for method in METHODS:
                if (dataset,method,replicate) in seen: continue
                choice=chosen[dataset,method]
                progress(stage='permutation',dataset=dataset,method=method,replicate=replicate)
                result=fit_model(method,X[order],choice['radius'],choice['min_pts'],1.4)
                labels=result.pop('labels'); aligned=np.empty_like(labels); aligned[order]=labels
                reference=np.load(ROOT/choice['labels_path'])['labels']
                target=RESULTS/'assignments'/f'{dataset}_{method}_perm{replicate}.npz'
                np.savez_compressed(target,labels=aligned)
                row=dict(dataset=dataset,method=method,replicate=replicate,status='ok',origin='fresh',
                         **{k:choice[k] for k in ('radius','radius_index','min_pts','merge_scale')},
                         **result,**score(y,aligned),stability_ari=float(adjusted_rand_score(reference,aligned)),
                         permutation_sha256=hashlib.sha256(order.tobytes()).hexdigest(),
                         labels_path=str(target.relative_to(ROOT)))
                append(RESULTS/'permutation.jsonl',row)

def main(argv=None):
    global RESULTS,GRID
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results-dir',type=Path,default=ROOT/'results')
    parser.add_argument('--pilot-only',action='store_true')
    args=parser.parse_args(argv)
    RESULTS=args.results_dir.resolve()
    assert RESULTS.is_relative_to(ROOT) and RESULTS!=ROOT
    RESULTS.mkdir(parents=True,exist_ok=True)
    manifest=dict(protocol_sha256=sha(ROOT/'protocol.json'),run_source_sha256=sha(__file__),
                  fit_adapter_sha256=sha(HELPERS/'fit_adapter.py'),native=native_details(),
                  imported_mode=False,roster=ROSTER)
    old=read(RESULTS/'run_manifest.json')
    if old: assert old==manifest,'Changed run mode/source/environment; use a fresh results directory.'
    else:
        assert not (RESULTS/'grid.jsonl').exists()
        save(RESULTS/'run_manifest.json',manifest)
    GRID=lines(RESULTS/'grid.jsonl')
    try:
        pilot()
        if args.pilot_only:
            progress(stage='pilot_complete'); return
        for dataset in ROSTER: grid(dataset)
        permutations()
        progress(stage='complete',datasets=len(ROSTER),selected_configurations=len(read(RESULTS/'selected.json')),
                 grid_origins=dict(Counter(r['origin'] for r in GRID)),
                 permutation_origins=dict(Counter(r['origin'] for r in lines(RESULTS/'permutation.jsonl'))))
    except Exception as error:
        append(RESULTS/'failures.jsonl',dict(status='failed',exception=type(error).__name__,message=str(error),
               last_progress=read(RESULTS/'progress.json'),utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())))
        raise

if __name__=='__main__': main()
