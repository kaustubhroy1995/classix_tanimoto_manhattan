"""AMI and sampled-query/full-reference Tanimoto silhouettes of SAVED fits.

No estimator is imported or refitted. All writes stay in --output-dir.
Repeated partitions must agree; differing repetitions require an explicit analysis.
"""
from __future__ import annotations
from pathlib import Path
import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "common"))
from experiment_paths import REPO, COMMON, VENDOR, work_root, read_configuration, read_repetitions
SOURCE_ROOT = Path(__file__).resolve().parents[1]
CONFIG = SOURCE_ROOT / "config.py"
import os
for key in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMBA_NUM_THREADS'):
    os.environ[key] = '1'
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import subprocess
import time
import numpy as np
from sklearn.metrics import adjusted_mutual_info_score, adjusted_rand_score

CODE_ROOT = SOURCE_ROOT
ROOT = CODE_ROOT / 'derived/scores'
RESULTS = ROOT / 'results'
ARGS = None
METHODS = ['classix_packed','dbscan','bitbirch','bblean','butina']
SEED, SAMPLE_SIZE, BLOCK_SIZE = 101, 2000, 64


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda:stream.read(1024**2),b''): h.update(chunk)
    return h.hexdigest()


def array_sha(array):
    return hashlib.sha256(memoryview(np.ascontiguousarray(array)).cast('B')).hexdigest()


def save(path, value):
    path = Path(path)
    assert path.resolve().is_relative_to(ROOT.resolve()), path
    path.parent.mkdir(parents=True, exist_ok=True)
    temp = path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False)+'\n',encoding='utf-8')
    for attempt in range(60):
        try:
            temp.replace(path)
            break
        except PermissionError:
            if attempt == 59: raise
            time.sleep(.05)


def utc():
    return datetime.now(timezone.utc).isoformat()


def progress(stage, **fields):
    value = dict(stage=stage,pid=os.getpid(),updated_utc=utc(),**fields)
    save(ROOT/'progress.json', value)
    print(json.dumps(value), flush=True)
    if (ROOT/'STOP').exists(): raise SystemExit('Paused at a safe checkpoint; remove STOP explicitly to resume.')


def canonical(labels):
    _,first,inverse = np.unique(labels,return_index=True,return_inverse=True)
    order = np.argsort(first)
    remap = np.empty(len(order), dtype=np.int64)
    remap[order] = np.arange(len(order), dtype=np.int64)
    return remap[inverse]


def policy(scope, method):
    return 'nonnoise_cluster_size_strictly_greater_than_10' if scope=='landrum' and method in ('bitbirch','bblean') else 'all_nonnoise_clusters'


def retained(labels, rule):
    keep = labels >= 0
    names, counts = np.unique(labels[keep],return_counts=True)
    if rule=='nonnoise_cluster_size_strictly_greater_than_10':
        keep &= np.isin(labels,names[counts>10])
    names,inverse,counts = np.unique(labels[keep], return_inverse=True,return_counts=True)
    compact = np.full(len(labels),-1,dtype=np.int32)
    compact[keep] = inverse.astype(np.int32)
    return keep,compact,counts.astype(np.int64)


def prepare_inventory():
    """Check raw evidence and aligned arrays; derive exact AMI without clustering."""
    RESULTS.mkdir(parents=True,exist_ok=True)
    sources, datasets, records, states = {}, {}, [], {}
    ami_cache = {}
    def source(path):
        relative = str(Path(path).resolve())
        if relative not in sources: sources[relative] = sha(path)
        return relative
    source(Path(__file__))
    protected = {}
    def dataset(scope,case):
        key = scope+'_'+case['data_id']
        if key in datasets: return datasets[key]
        directory = ARGS.data_dir/case['data_id']
        manifest_path = directory/'manifest.json'; manifest=read(manifest_path)
        assert manifest['shape']==[case['rows'],case['features']]
        if scope=='synthetic': assert manifest['parameters']==case['generator']
        arrays = {}
        for kind in ('packed_little','truth'):
            path = directory/manifest['files'][kind]['file']
            value = np.load(path,allow_pickle=False,mmap_mode='r')
            assert array_sha(value)==manifest['files'][kind]['sha256'], (key,kind)
            arrays[kind]=value; source(path)
        packed,truth = arrays['packed_little'],arrays['truth']
        assert packed.dtype==np.uint8 and packed.shape==(case['rows'],(case['features']+7)//8)
        assert truth.shape==(case['rows'],) and truth.dtype==np.int32
        if case['features']%8:
            assert not np.any(packed[:,-1] & np.uint8(255 ^ ((1<<(case['features']%8))-1)))
        rows = min(ARGS.pilot_rows, case['rows']) if ARGS.pilot_rows else case['rows']
        packed, truth = packed[:rows], truth[:rows]
        logical_hash = manifest['files']['int32']['sha256']
        if rows != manifest['shape'][0]:
            dense = np.unpackbits(packed, axis=1, bitorder='little', count=case['features']).astype(np.int32)
            logical_hash = array_sha(dense)
        item = dict(key=key,scope=scope,data_id=case['data_id'],rows=rows,features=case['features'],
                    directory=str(directory.resolve()),
                    packed_file=source(directory/'packed_little.npy'), truth_file=source(directory/'truth.npy'),
                    manifest_file=source(manifest_path),manifest_sha256=sha(manifest_path),
                    packed_array_sha256=array_sha(packed),truth_array_sha256=array_sha(truth),
                    logical_int32_array_sha256=logical_hash,
                    reference_classes=int(np.unique(truth).size))
        datasets[key]=item
        return item
    plan = read_configuration(ARGS.config); source(ARGS.config)
    cases = plan['cases']
    if ARGS.cases:
        assert ARGS.pilot_rows, 'A reduced case grid must be explicitly marked as pilot.'
        requested = set(ARGS.cases)
        assert len(requested) == len(ARGS.cases)
        assert requested.issubset({c['case_id'] for c in cases})
        cases = [c for c in cases if c['case_id'] in requested]
    assert bool(plan.get('diagnostic_pilot', False)) == bool(ARGS.pilot_rows) or ARGS.scope == 'synthetic'
    methods = METHODS if ARGS.scope == 'synthetic' else ['classix_packed','bitbirch','bblean']
    contract = dict(scope=ARGS.scope, is_pilot=bool(ARGS.pilot_rows), pilot_rows=ARGS.pilot_rows,
                    expected_cases=[c['case_id'] for c in cases], expected_methods=methods,
                    complete_case_results=True)
    from render_scores import summary_contract
    summary_contract(contract)
    amendment = read_repetitions(ARGS.config)
    for scope in [ARGS.scope]:
        for i,case in enumerate(cases):
            progress('inventory_and_exact_ami',scope=scope,case=i+1,total_cases=len(cases),case_id=case['case_id'])
            input_info = dataset(scope,case)
            truth=np.load(input_info['truth_file'],allow_pickle=False)[:input_info['rows']]
            raw_path = ARGS.results_dir/(case['case_id']+'.json')
            raw = read(raw_path); source(raw_path)
            assert raw['case']==case and raw.get('completed') is True
            if raw.get('data_manifest_sha256'):
                assert raw['data_manifest_sha256']==input_info['manifest_sha256']
            if raw.get('reproduction'):
                assert all(run.get('protocol')==raw['protocol'] for run in raw['runs'])
            if scope == 'synthetic': assert raw.get('pilot_rows') == ARGS.pilot_rows
            elif not ARGS.pilot_rows: assert raw.get('pilot_rows') is None

            for method in METHODS if scope=='synthetic' else ['classix_packed','bitbirch','bblean']:
                runs = [r for r in raw['runs'] if r['variant']==method]
                schedules = ({m:[0] for m in methods} if ARGS.pilot_rows else
                             case.get('method_repeat_ids') if scope=='landrum' else
                             amendment.get('schedules', {}).get(case['case_id']))
                assert schedules is not None, 'Missing authoritative repeat schedule.'
                if raw.get('method_repeat_ids') is not None:
                    assert raw['method_repeat_ids']==schedules, 'Saved schedule differs from the intended schedule.'
                expected = schedules[method]
                assert [r['repeat'] for r in runs]==expected[:len(runs)]
                failed = [r for r in runs if r['status']!='ok']
                base=dict(scope=scope,case_id=case['case_id'],data_key=input_info['key'],method=method,
                          study=case['study'],x=case['x'] if scope=='synthetic' else case['grid']['similarity'],
                          n_total=input_info['rows'],n_features=case['features'],source_record=source(raw_path),
                          rule=policy(scope,method),case_parameters=case)
                if failed:
                    assert len(failed)==1 and failed[0] is runs[-1]
                    assert failed[0]['status'] not in ('pending','running','ok'), 'Fit is not terminal.'
                    records.append(dict(base,status=failed[0]['status'],ami_all=None,ami_retained=None,
                                        silhouette_mean=None,silhouette_se=None,silhouette_status='not_available',
                                        silhouette_reason='No complete measured repetition schedule: '+str(failed[0]['status'])+'.',
                                        n_retained=None,coverage=None,n_clusters=None,n_queries=0,exact_silhouette=False))
                    continue
                assert [r['repeat'] for r in runs]==expected
                measured=[r for r in runs if not r['warmup']]
                assert measured and all(r['repeat']>=0 for r in measured)
                reference=None; evidence=[]; partitions=set()
                for run in measured:
                    checked=run['validated']
                    relative = checked.get('output_relative_assignments_file')
                    labels_path = ARGS.results_dir.parent/relative if relative else Path(checked['assignments_file'])
                    assert labels_path.is_file() and labels_path.resolve().is_relative_to(raw_path.parent.resolve())
                    labels=np.load(labels_path,allow_pickle=False)
                    assert labels.dtype==np.int64 and labels.shape==truth.shape and np.all(labels>=-1)
                    canonical_hash=array_sha(canonical(labels))
                    assert canonical_hash==checked['canonical_partition_sha256']
                    assert checked['logical_input_sha256']==input_info['logical_int32_array_sha256'] and (checked['input_unchanged'] is True or (ARGS.pilot_rows and checked['input_unchanged'] is None))
                    assert int((labels<0).sum())==checked['noise_points']
                    assert int(np.unique(labels[labels>=0]).size)==checked['clusters']
                    ari=float(adjusted_rand_score(truth,labels))
                    assert abs(ari-checked['ari'])<1e-12
                    noise_hash=array_sha(labels<0)
                    partitions.add((canonical_hash,noise_hash))
                    evidence.append(dict(repeat=run['repeat'],assignment_file=source(labels_path),
                                         assignment_file_sha256=sha(labels_path),assignment_array_sha256=array_sha(labels),
                                         canonical_partition_sha256=canonical_hash,noise_mask_sha256=noise_hash,ari=ari))
                    if reference is None: reference=labels
                assert len(partitions)==1, 'Repeated partitions differ; score every partition rather than choosing one.'
                canon,noise=next(iter(partitions)); rule=base['rule']
                state_key=hashlib.sha256((input_info['key']+'|'+input_info['packed_array_sha256']+'|'+input_info['truth_array_sha256']+'|'+canon+'|'+noise+'|'+rule).encode()).hexdigest()[:24]
                if state_key not in states:
                    keep,compact,counts=retained(reference,rule)
                    retained_n=int(keep.sum()); cluster_n=len(counts)
                    key=(input_info['key'],canon,noise)
                    if key not in ami_cache:
                        ami_cache[key]=float(adjusted_mutual_info_score(truth,reference,average_method='arithmetic'))
                    valid=2<=cluster_n<=retained_n-1
                    reason=None if valid else ('Fewer than two retained clusters.' if cluster_n<2 else 'Every retained point is a singleton cluster.')
                    order=np.random.default_rng(SEED).permutation(len(reference))
                    query=order[keep[order]][:min(SAMPLE_SIZE,retained_n)] if valid else np.empty(0,dtype=np.int64)
                    query=np.sort(query).astype(np.int64)
                    target=RESULTS/'states'/f'{state_key}.npz'
                    target.parent.mkdir(parents=True,exist_ok=True)
                    np.savez_compressed(target,labels=compact,counts=counts,queries=query)
                    states[state_key]=dict(state_key=state_key,data_key=input_info['key'],rule=rule,
                        canonical_partition_sha256=canon,noise_mask_sha256=noise,
                        source_assignment=evidence[0],state_file=target.relative_to(ROOT).as_posix(),
                        compact_array_sha256=array_sha(compact),retained_mask_sha256=array_sha(keep),
                        counts_array_sha256=array_sha(counts),query_array_sha256=array_sha(query),
                        n_total=len(reference),n_retained=retained_n,coverage=retained_n/len(reference),
                        n_clusters=cluster_n,n_queries=len(query),exact_silhouette=len(query)==retained_n,
                        ami_all=ami_cache[key],
                        ami_retained=float(adjusted_mutual_info_score(truth[keep],reference[keep],average_method='arithmetic')) if retained_n else None,
                        silhouette_status='pending' if valid else 'undefined',silhouette_reason=reason,
                        silhouette_mean=None,silhouette_se=None)
                state=states[state_key]
                fields={k:state[k] for k in ('ami_all','ami_retained','n_retained','coverage','n_clusters','n_queries',
                        'exact_silhouette','silhouette_status','silhouette_reason','silhouette_mean','silhouette_se')}
                records.append(dict(base,status='ok',state_key=state_key,measured_repeats=len(measured),
                                    repeated_partitions_identical=True,assignments=evidence,**fields))
    assert len(records)==len(cases)*len(methods)
    protocol=dict(seed=SEED,query_limit=SAMPLE_SIZE,block_size=BLOCK_SIZE,
        ami='sklearn adjusted_mutual_info_score, arithmetic normalizer, full reference labels/assignments including noise as the recorded -1 category; no sample and no tuning.',
        silhouette='Tanimoto/Jaccard distance on binary fingerprints; exact distances to every retained cluster member for each uniformly selected query. Self excluded, singleton silhouette zero, a=b=0 gives zero; valid only for 2 <= K <= N_retained-1.',
        query_sampling='One shared random row priority per dataset from default_rng(101). For each retained population take its first min(2000,N_retained) rows; sort query IDs for computation. Equal populations use identical queries. Every population is a uniform sample without replacement. Query unions allow distance reuse without changing any mean.',
        estimate='Mean of sampled individual silhouettes estimates the full retained-population mean, not a silhouette computed only within a sampled submatrix.',
        uncertainty='Design-based Monte Carlo standard error = sample SD/sqrt(m)*sqrt(1-m/N_retained). Exact when m=N_retained; pointwise approximate normal bands describe query sampling only, not clustering stability or statistical significance.',
        populations={'synthetic':'Exclude recorded noise (-1); keep all nonnoise clusters, including singletons.',
                     'landrum_classix':'All final nonnoise clusters, matching the current minPts=10 study.',
                     'landrum_bitbirch_lean':'Keep final clusters strictly larger than 10 (>10), matching the existing minimum-size quality evaluation. Do not reinterpret as >=10.'},
        aggregation='All measured assignments are checked. Repeated partitions must agree; otherwise scoring stops explicitly for a repetition-aware analysis, rather than choosing a fit. Score each distinct partition/population once. Repeated deterministic fits are not independent observations.',
        scope=f'{ARGS.scope}: {len(cases)} settings, {len(methods)} methods; pilot={bool(ARGS.pilot_rows)}. No estimator fits or tuning.',
        definitions={'silhouette':'https://scikit-learn.org/1.6/modules/generated/sklearn.metrics.silhouette_samples.html',
                     'ami':'https://scikit-learn.org/1.6/modules/generated/sklearn.metrics.adjusted_mutual_info_score.html'})
    inventory=dict(schema=2,**contract,created_utc=utc(),protocol=protocol,datasets=datasets,states=states,records=records,
                   input_file_sha256=sources,protected_file_sha256=protected,
                   environment=dict(python=sys.version,platform=platform.platform(),
                                    packages={p:importlib.metadata.version(p) for p in ('numpy','scipy','scikit-learn','matplotlib','threadpoolctl')}))
    save(ROOT/'protocol.json',protocol)
    save(RESULTS/'inventory.json',inventory)
    refresh_summary(inventory)
    progress('exact_ami_complete',records=len(records),complete_settings=sum(r['status']=='ok' for r in records),failed_settings=sum(r['status']!='ok' for r in records),distinct_states=len(states))
    return inventory


def scoring_identity(inventory):
    identity={name:sha((SOURCE_ROOT/'code'/name) if (SOURCE_ROOT/'code'/name).exists() else (COMMON/name)) for name in ('saved_scores.py','native_silhouette.py','silhouette_kernel.c')}
    identity['environment']=inventory['environment']
    identity['silhouette_kernel.dll']=sha(Path(os.environ['CLASSIX_SILHOUETTE_DIR'])/'silhouette_kernel.dll')
    return identity


def refresh_summary(inventory):
    records=[]
    producer=scoring_identity(inventory)
    for record in inventory['records']:
        record=dict(record)
        if record['status']=='ok':
            score_path=RESULTS/'scores'/(record['state_key']+'.json')
            if score_path.exists():
                score=read(score_path)
                state=inventory['states'][record['state_key']]
                valid=(score.get('scoring_source')==producer and all(score.get(k)==state[k] for k in
                       ('query_array_sha256','compact_array_sha256','counts_array_sha256')))
                if valid:
                    for key in ('silhouette_mean','silhouette_se','silhouette_status','silhouette_reason'):
                        record[key]=score[key]
        records.append(record)
    summary=dict(schema=2,**{k:inventory[k] for k in ('scope','is_pilot','pilot_rows','expected_cases','expected_methods','complete_case_results')},updated_utc=utc(),protocol=inventory['protocol'],records=records,
                 counts={'settings':len(records),'successful_fits':sum(r['status']=='ok' for r in records),
                         'silhouette_status':dict(Counter(r['silhouette_status'] for r in records))},
                 limitations=['AMI uses synthetic truth or Landrum target categories; Landrum target categories need not uniquely define molecular clusters. The same SMILES may occur under multiple targets.',
                              'Silhouette uses method-specific retained populations; read coverage alongside scores. All-point AMI has the same full population for every method in each setting.',
                              'Sample-query silhouette is explicitly an estimate unless all retained rows were evaluated. No clustering reruns or timing claims are made.'])
    save(RESULTS/'summary.json',summary)


def evaluate(inventory, only_datasets=None):
    from native_silhouette import distance_block, prepare_packed, prepare_labels, reduce_block
    producer={name:sha((SOURCE_ROOT/'code'/name) if (SOURCE_ROOT/'code'/name).exists() else (COMMON/name)) for name in ('saved_scores.py','native_silhouette.py','silhouette_kernel.c')}
    producer = scoring_identity(inventory)
    save(RESULTS/'native_build.json', {p.name:sha(p) for p in Path(os.environ['CLASSIX_SILHOUETTE_DIR']).glob('*') if p.is_file()})
    save(RESULTS/'scoring_source.json',producer)
    grouping=defaultdict(list)
    for key,state in inventory['states'].items(): grouping[state['data_key']].append(state)
    for dataset_index,(data_key,states) in enumerate(grouping.items()):
        if only_datasets and data_key not in only_datasets: continue
        info=inventory['datasets'][data_key]
        packed=np.load(info['packed_file'],allow_pickle=False)[:info['rows']]
        assert array_sha(packed)==info['packed_array_sha256']
        reference=prepare_packed(packed,n_features=info['features'])
        active=[]
        for state in states:
            target=RESULTS/'scores'/(state['state_key']+'.json')
            if target.exists():
                old=read(target)
                assert old['scoring_source']==producer and all(old[k]==state[k] for k in ('query_array_sha256','compact_array_sha256','counts_array_sha256')), 'Stale score cache; use a fresh output directory.'
                continue
            if state['silhouette_status']=='undefined':
                save(target,dict(state,scoring_source=producer,computed_utc=utc()))
                continue
            with np.load(ROOT/state['state_file'],allow_pickle=False) as values:
                labels,counts,queries=values['labels'],values['counts'],values['queries']
            assert array_sha(labels)==state['compact_array_sha256'] and array_sha(counts)==state['counts_array_sha256'] and array_sha(queries)==state['query_array_sha256']
            active.append(dict(state=state,labels=labels,counts=counts,native_labels=prepare_labels(labels,counts),queries=queries,
                               silhouette=np.full(len(queries),np.nan),a=np.full(len(queries),np.nan),b=np.full(len(queries),np.nan)))
        if not active:
            refresh_summary(inventory)
            continue
        union=np.unique(np.concatenate([item['queries'] for item in active]))
        started=time.perf_counter()
        progress('silhouette_dataset',dataset=data_key,dataset_index=dataset_index+1,datasets=len(grouping),
                 rows=info['rows'],features=info['features'],partitions=len(active),union_queries=len(union))
        for start in range(0,len(union),BLOCK_SIZE):
            query=union[start:start+BLOCK_SIZE]
            distances=distance_block(reference,query)
            for item in active:
                where=np.searchsorted(item['queries'],query)
                good=where<len(item['queries'])
                good[good] &= item['queries'][where[good]]==query[good]
                if not good.any(): continue
                positions=where[good]
                block=distances if good.all() else np.ascontiguousarray(distances[good])
                q=query if good.all() else query[good]
                values,a,b=reduce_block(block,item['native_labels'],q,check_distances=False)
                item['silhouette'][positions]=values; item['a'][positions]=a; item['b'][positions]=b
            if start==0 or (start//BLOCK_SIZE+1)%8==0:
                progress('silhouette_block',dataset=data_key,queries_done=min(start+BLOCK_SIZE,len(union)),
                         total_queries=len(union),partitions=len(active),elapsed_seconds=time.perf_counter()-started)
        for item in active:
            state=item['state']; values=item['silhouette']; n=state['n_retained']; m=len(values)
            assert np.isfinite(values).all() and np.all((values>=-1-1e-12)&(values<=1+1e-12))
            mean=float(values.mean())
            se=0.0 if m==n else float(values.std(ddof=1)/np.sqrt(m)*np.sqrt(1-m/n))
            arrays_path=RESULTS/'samples'/(state['state_key']+'.npz')
            arrays_path.parent.mkdir(parents=True,exist_ok=True)
            np.savez_compressed(arrays_path,queries=item['queries'],silhouette=values,a=item['a'],b=item['b'])
            save(RESULTS/'scores'/(state['state_key']+'.json'),dict(state,silhouette_status='ok',silhouette_reason=None,
                 silhouette_mean=mean,silhouette_se=se,scoring_source=producer,computed_utc=utc(),
                 samples_file=arrays_path.relative_to(ROOT).as_posix(),samples_file_sha256=sha(arrays_path),
                 sample_silhouette_array_sha256=array_sha(values),dataset_elapsed_seconds=time.perf_counter()-started))
        refresh_summary(inventory)
        progress('silhouette_dataset_complete',dataset=data_key,elapsed_seconds=time.perf_counter()-started)
    refresh_summary(inventory)
    remaining=sum(row['silhouette_status']=='pending' for row in read(RESULTS/'summary.json')['records'])
    progress('complete' if remaining==0 else 'partial',remaining_silhouette_settings=remaining)


def main():
    global ARGS, ROOT, RESULTS
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--scope', choices=['synthetic','landrum'], required=True)
    for name in ('config','results-dir','data-dir','output-dir'):
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--pilot-rows',type=int)
    parser.add_argument('--cases',nargs='+')
    ARGS=parser.parse_args()
    if ARGS.pilot_rows is not None: assert ARGS.pilot_rows>=2
    for name in ('config','results_dir','data_dir','output_dir'):
        setattr(ARGS,name,getattr(ARGS,name).resolve())
    ROOT=ARGS.output_dir; RESULTS=ROOT/'results'
    assert ROOT != CODE_ROOT
    assert not ROOT.is_relative_to(ARGS.data_dir) and not ROOT.is_relative_to(ARGS.results_dir)
    assert not list(ARGS.results_dir.parent.rglob('run.lock')), 'Timed fits must finish before scoring.'
    inventory=prepare_inventory()  # Recheck input and assignment integrity on every invocation.
    evaluate(inventory)
    for path, expected in inventory['input_file_sha256'].items():
        assert sha(path)==expected, 'Input changed during scoring: '+path
    summary=read(RESULTS/'summary.json')
    assert not any(r['silhouette_status']=='pending' for r in summary['records'])
    subprocess.run([sys.executable,'-B',str(CODE_ROOT/'code/render_scores.py'),'--root',str(ROOT)],check=True)
    save(ROOT/'verification.json',dict(status='passed',checked_inputs=len(inventory['input_file_sha256']),
         settings=len(summary['records']),summary_sha256=sha(RESULTS/'summary.json'),
         definitions='Arithmetic AMI and sampled-query/full-reference Tanimoto silhouette.',
         timed_fits_rerun=False,checked_utc=utc()))


if __name__=='__main__': main()
