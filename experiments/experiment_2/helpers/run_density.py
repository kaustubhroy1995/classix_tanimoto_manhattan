"""Seeded, serial density baselines; search caches are never reported as fit times."""
from __future__ import annotations
import os
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[name]='1'
import argparse, hashlib, inspect, json, platform, sys, time, warnings
from pathlib import Path
from paths import REPO, HELPERS, RUN, DENSITY, MANHATTAN, PROVENANCE
import numpy as np
import scipy, sklearn
from scipy.spatial.distance import pdist, squareform
from sklearn.cluster import DBSCAN, OPTICS, cluster_optics_xi
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from threadpoolctl import threadpool_limits
from fit_adapter import fit_model, native_details

ROOT=DENSITY
RESULTS=ROOT/'results'
POLICY=json.loads((ROOT/'protocol.json').read_text())
METRICS=POLICY['metrics']
ALGORITHMS=('dbscan','optics')
GRID=[]

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def read(path,default=None):return json.loads(Path(path).read_text()) if Path(path).exists() else default
def lines(path):return [json.loads(x) for x in Path(path).read_text().splitlines() if x.strip()] if Path(path).exists() else []
def save(path,value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
    tmp.replace(path)
def append(path,row):
    with Path(path).open('a') as f:f.write(json.dumps(row,allow_nan=False)+'\n')
def progress(**values):
    values.update(pid=os.getpid(),updated_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),grid_candidates=len(GRID))
    save(RESULTS/'progress.json',values);print(json.dumps(values),flush=True)
    if (ROOT/'STOP').exists():raise SystemExit('Paused at checkpoint; remove STOP to resume.')
def canonical(labels):
    ids={};out=[]
    for v in labels:
        v=int(v)
        if v == -1:out.append(-1)
        else:
            if v not in ids:ids[v]=len(ids)
            out.append(ids[v])
    return np.asarray(out,dtype='<i8')
def partition_sha(labels):return hashlib.sha256(canonical(labels).tobytes()).hexdigest()
def score(y,labels):
    assert labels.shape==y.shape and np.all(labels>=-1)
    return dict(ari=float(adjusted_rand_score(y,labels)),ami=float(adjusted_mutual_info_score(y,labels,average_method='arithmetic')),
                n_clusters=int(np.unique(labels[labels>=0]).size),noise_points=int(np.count_nonzero(labels==-1)),
                retained_fraction=float(np.mean(labels>=0)),partition_sha256=partition_sha(labels))
def load(dataset):
    path=ROOT/'data'/f'{dataset}.npz'
    expected=read(ROOT/'input_provenance.json')['inputs'][f'data/{dataset}.npz']
    assert sha(path)==expected
    with np.load(path,allow_pickle=False) as a:X,y=a['X'].copy(),a['y'].copy()
    assert X.dtype==np.float64 and X.ndim==2 and len(X)==len(y) and np.isfinite(X).all()
    return X,y,path
def build(algorithm,metric,minimum,parameter):
    shared=dict(metric=metric,min_samples=int(minimum),algorithm='auto',n_jobs=1)
    if algorithm=='dbscan':return DBSCAN(eps=float(parameter),**shared)
    return OPTICS(xi=float(parameter),max_eps=np.inf,cluster_method='xi',min_cluster_size=None,
                  predecessor_correction=True,memory=None,**shared)
def density_fit(algorithm,metric,X,minimum,parameter):
    data=np.array(X,dtype=np.float64,order='C',copy=True)
    model=build(algorithm,metric,minimum,parameter)
    with threadpool_limits(limits=1), warnings.catch_warnings(record=True) as seen:
        warnings.simplefilter('always')
        start=time.perf_counter();model.fit(data);seconds=time.perf_counter()-start
    return model,np.asarray(model.labels_,dtype=np.int64).copy(),float(seconds),[
        {'category':w.category.__name__,'message':str(w.message)} for w in seen]
def plan(dataset,X,path):
    target=RESULTS/'plans'/f'{dataset}.json'
    found=read(target)
    if found:
        assert found['input_sha256']==sha(path);return found
    data=dict(dataset=dataset,input_sha256=sha(path),shape=list(X.shape),methods={})
    for metric in METRICS:
        values=pdist(X,metric='cityblock' if metric=='manhattan' else 'euclidean')
        dense=squareform(values)
        near=float(np.median(np.partition(dense,5,axis=1)[:,5]))
        far=float(np.median(values[values>0]))
        low,high=.35*near,.9*far
        assert 0<low<high and np.isfinite(high)
        data['methods'][metric]=dict(eps_low=low,eps_high=high,median_5nn=near,median_pair=far)
    save(target,data);return data
def parameter(algorithm,metric,design,index):
    if algorithm=='dbscan':low,high=design['methods'][metric]['eps_low'],design['methods'][metric]['eps_high']
    else:low,high=POLICY['optics']['xi_low'],POLICY['optics']['xi_high']
    return float(low*(high/low)**(index/31))
def best(rows):return min(rows,key=lambda r:(-r['ari'],r['min_samples'],r['parameter_index']))
def optics_graph(dataset,metric,minimum,X,input_sha):
    path=RESULTS/'graphs'/f'{dataset}_{metric}_{minimum}.npz'
    meta=path.with_suffix('.json')
    record=read(meta)
    if record:
        assert record['input_sha256']==input_sha and sha(path)==record['graph_sha256']
        with np.load(path,allow_pickle=False) as a:graph={k:a[k].copy() for k in a.files}
        return graph,record
    progress(stage='optics_graph',dataset=dataset,metric=metric,min_samples=minimum)
    model,labels,seconds,seen=density_fit('optics',metric,X,minimum,.05)
    graph={k:np.asarray(getattr(model,k+'_')).copy() for k in ['reachability','predecessor','ordering','core_distances']}
    check,_=cluster_optics_xi(reachability=graph['reachability'],predecessor=graph['predecessor'],ordering=graph['ordering'],min_samples=minimum,xi=.05,min_cluster_size=None,predecessor_correction=True)
    assert np.array_equal(check,labels)
    path.parent.mkdir(exist_ok=True)
    np.savez_compressed(path,**graph)
    record=dict(dataset=dataset,metric=metric,min_samples=minimum,input_sha256=input_sha,
                graph_sha256=sha(path),search_graph_fit_seconds=seconds,warnings=seen)
    save(meta,record);return graph,record
def run_candidates(dataset,X,y,design,algorithm,indices,minima):
    seen={(r['dataset'],r['algorithm'],r['metric'],r['min_samples'],r['parameter_index']) for r in GRID}
    for metric in METRICS:
        for minimum in minima:
            pending=[i for i in indices if (dataset,algorithm,metric,minimum,i) not in seen]
            if not pending:continue
            progress(stage='search',dataset=dataset,algorithm=algorithm,metric=metric,min_samples=minimum)
            graph=None
            if algorithm=='optics':graph,_=optics_graph(dataset,metric,minimum,X,design['input_sha256'])
            for index in pending:
                value=parameter(algorithm,metric,design,index)
                if graph is None:
                    _,labels,seconds,observed=density_fit(algorithm,metric,X,minimum,value)
                    timing_key='search_fit_seconds'
                else:
                    start=time.perf_counter()
                    labels,_=cluster_optics_xi(reachability=graph['reachability'],predecessor=graph['predecessor'],ordering=graph['ordering'],min_samples=minimum,xi=value,min_cluster_size=None,predecessor_correction=True)
                    seconds=time.perf_counter()-start;observed=[];timing_key='search_extraction_seconds'
                row=dict(dataset=dataset,algorithm=algorithm,metric=metric,min_samples=minimum,parameter=value,
                         parameter_name='eps' if algorithm=='dbscan' else 'xi',parameter_index=index,
                         status='ok',warnings=observed,**score(y,labels),**{timing_key:seconds})
                append(RESULTS/'grid.jsonl',row);GRID.append(row)
                if (ROOT/'STOP').exists():progress(stage='paused_after_candidate',dataset=dataset)
def sweep(dataset,X,y,design,algorithm):
    run_candidates(dataset,X,y,design,algorithm,list(range(32)),POLICY['min_samples'])
    path=RESULTS/'plans'/f'{dataset}_{algorithm}_extension.json'
    ext=read(path)
    if ext is None:
        winners={m:best([r for r in GRID if r['dataset']==dataset and r['algorithm']==algorithm and r['metric']==m and 0<=r['parameter_index']<=31 and r['min_samples'] in POLICY['min_samples']]) for m in METRICS}
        edge=any(r['parameter_index'] in (0,31) for r in winners.values())
        min_edge=any(r['min_samples']==80 for r in winners.values())
        ext=dict(parameter_edge=edge,min_samples_edge=min_edge,initial_winners=winners,
                 indices=list(range(-4,36)) if edge else list(range(32)),
                 min_samples=POLICY['min_samples']+[v for v in (120,160) if min_edge and v<=len(X)])
        save(path,ext)
    run_candidates(dataset,X,y,design,algorithm,ext['indices'],ext['min_samples'])
def select_and_time(dataset,X,y):
    target=RESULTS/'selected'/f'{dataset}.json'
    if target.exists():return
    prior=read(RESULTS/'classix_selected_prior.json')
    configs=[]
    for metric in METRICS:
        old=next(r for r in prior if r['dataset']==dataset and r['method']==metric)
        configs.append(dict(dataset=dataset,algorithm='classix',metric=metric,radius=old['radius'],min_pts=old['min_pts'],merge_scale=old['merge_scale'],ari=old['ari'],ami=old['ami'],selection_origin='existing_Experiment2_grid',prior_fit_seconds_median=old['fit_seconds_median']))
    for algorithm in ALGORITHMS:
        for metric in METRICS:
            winner=best([r for r in GRID if r['dataset']==dataset and r['algorithm']==algorithm and r['metric']==metric])
            configs.append(dict(winner,selection_origin='new_Experiment2_density_grid'))
    progress(stage='selected_fresh_timings',dataset=dataset,variants=len(configs))
    def fit(c):
        if c['algorithm']=='classix':
            r=fit_model(c['metric'],X,c['radius'],c['min_pts'],c['merge_scale'])
            return r.pop('labels'),r['fit_seconds'],r['warnings']
        _,labels,seconds,seen=density_fit(c['algorithm'],c['metric'],X,c['min_samples'],c['parameter'])
        return labels,seconds,seen
    labels_by_key={}
    for c in configs:
        labels,_,observed=fit(c);measured=score(y,labels)
        assert abs(measured['ari']-c['ari'])<1e-12 and abs(measured['ami']-c['ami'])<1e-12
        if c['algorithm']!='classix':assert c['partition_sha256']==measured['partition_sha256']
        c.update(measured,timing_seconds=[],timing_warnings=observed)
        labels_by_key[(c['algorithm'],c['metric'])]=labels
    for rep in range(5):
        order=np.random.default_rng(np.random.SeedSequence([101,POLICY['datasets'].index(dataset),rep,991])).permutation(len(configs))
        for index in order:
            c=configs[int(index)];labels,seconds,observed=fit(c)
            assert np.array_equal(canonical(labels),canonical(labels_by_key[(c['algorithm'],c['metric'])]))
            c['timing_seconds'].append(seconds);c['timing_warnings'].extend(observed)
    for c in configs:
        path=RESULTS/'assignments'/f"{dataset}_{c['algorithm']}_{c['metric']}.npz"
        path.parent.mkdir(exist_ok=True)
        np.savez_compressed(path,labels=labels_by_key[(c['algorithm'],c['metric'])])
        c.update(fit_seconds_median=float(np.median(c['timing_seconds'])),labels_path=path.relative_to(ROOT).as_posix(),labels_file_sha256=sha(path),input_sha256=sha(ROOT/'data'/f'{dataset}.npz'))
    save(target,configs)
def manifest():
    data=dict(created_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),python=sys.version,executable=sys.executable,
              platform=platform.platform(),numpy=np.__version__,scipy=scipy.__version__,sklearn=sklearn.__version__,
              native=native_details(),protocol_sha256=sha(ROOT/'protocol.json'),runner_sha256=sha(__file__),
              sklearn_source={name:{'path':inspect.getsourcefile(cls),'sha256':sha(inspect.getsourcefile(cls))} for name,cls in [('DBSCAN',DBSCAN),('OPTICS',OPTICS),('cluster_optics_xi',cluster_optics_xi)]})
    target=RESULTS/'run_manifest.json'
    previous=read(target)
    if previous:
        for key in ['protocol_sha256','runner_sha256','sklearn','numpy','scipy']:assert previous[key]==data[key],f'Run mismatch {key}'
    else:save(target,data)
def pilot():
    rows=[]
    for dataset in ['iris','segmentation','musk']:
        X,y,path=load(dataset);design=plan(dataset,X,path)
        for algorithm in ALGORITHMS:
            for metric in METRICS:
                value=parameter(algorithm,metric,design,16)
                model,labels,seconds,seen=density_fit(algorithm,metric,X,10,value)
                assert labels.shape==y.shape and np.all(labels>=-1)
                if algorithm=='optics':
                    check,_=cluster_optics_xi(reachability=model.reachability_,predecessor=model.predecessor_,ordering=model.ordering_,min_samples=10,xi=value,min_cluster_size=None,predecessor_correction=True)
                    assert np.array_equal(check,labels)
                row=dict(dataset=dataset,algorithm=algorithm,metric=metric,shape=list(X.shape),seconds=seconds,warnings=seen,unscored=True)
                rows.append(row);print(json.dumps(row),flush=True)
    save(RESULTS/'readiness.json',dict(status='passed',records=rows))
def main():
    parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');args=parser.parse_args()
    if args.pilot:pilot();return
    assert read(RESULTS/'readiness.json')['status']=='passed'
    manifest();GRID.extend(lines(RESULTS/'grid.jsonl'))
    for dataset in POLICY['datasets']:
        if (RESULTS/'selected'/f'{dataset}.json').exists():continue
        X,y,path=load(dataset);design=plan(dataset,X,path)
        for algorithm in ALGORITHMS:sweep(dataset,X,y,design,algorithm)
        select_and_time(dataset,X,y)
        progress(stage='dataset_complete',dataset=dataset,complete_datasets=len(list((RESULTS/'selected').glob('*.json'))),total_datasets=len(POLICY['datasets']))
    combined=[r for d in POLICY['datasets'] for r in read(RESULTS/'selected'/f'{d}.json')]
    save(RESULTS/'selected.json',combined)
    progress(stage='complete',complete_datasets=len(POLICY['datasets']),selected_variants=len(combined))
if __name__=='__main__':main()
