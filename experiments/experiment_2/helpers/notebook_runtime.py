"""Environment, worker execution and shared inputs for the numeric notebook."""
from pathlib import Path
import hashlib
import json
import os
import subprocess
import sys

from config import CLASSIX_CONFIG, DENSITY_CONFIG
from paths import REPO, HELPERS, RUN, DENSITY, ensure_source_provenance

NOTEBOOK_DIR = HELPERS.parent
PYTHON = NOTEBOOK_DIR / '_generated/environment' / ('Scripts/python.exe' if os.name == 'nt' else 'bin/python')

def checked_write(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and path.read_bytes() != data:
        raise RuntimeError(f'Inputs changed: {path}. Choose a new RUN_NAME.')
    if not path.exists():
        path.write_bytes(data)

def setup():
    RUN.mkdir(parents=True, exist_ok=True)
    DENSITY.mkdir(parents=True, exist_ok=True)
    requirements = NOTEBOOK_DIR / 'requirements.txt'
    identity = PYTHON.parents[1] / 'experiment_requirements.sha256'
    requirement_sha = hashlib.sha256(requirements.read_bytes()).hexdigest()
    if not PYTHON.is_file():
        candidates = [[sys.executable]] + ([['py', '-3.12']] if os.name == 'nt' else [['python3.12']])
        base = None
        for command in candidates:
            try:
                version = subprocess.check_output(command + ['-c', 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")'], text=True).strip()
                if version == '3.12':
                    base = command
                    break
            except (OSError, subprocess.CalledProcessError):
                pass
        if base is None:
            raise RuntimeError('Python 3.12 is required.')
        subprocess.run(base + ['-m', 'venv', str(PYTHON.parents[1])], check=True)
    if not identity.exists():
        subprocess.run([str(PYTHON), '-m', 'pip', 'install', '-r', str(requirements)], check=True)
        identity.write_text(requirement_sha, encoding='utf-8')
    elif identity.read_text(encoding='utf-8') != requirement_sha:
        raise RuntimeError('Dependency pins differ from this environment.')
    version_check = 'import importlib.metadata as m; import sys; assert sys.version_info[:2] == (3,12); ' + '; '.join(
        f"assert m.version({line.split('==')[0]!r}) == {line.split('==')[1]!r}"
        for line in requirements.read_text().splitlines() if line)
    subprocess.run([str(PYTHON), '-c', version_check], check=True)
    for target, policy in ((RUN, CLASSIX_CONFIG), (DENSITY, DENSITY_CONFIG)):
        checked_write(target/'protocol.json', (json.dumps(policy, indent=2) + '\n').encode())
    ensure_source_provenance()
    with (RUN/'environment.txt').open('w', encoding='utf-8') as stream:
        subprocess.run([str(PYTHON), '-m', 'pip', 'freeze'], stdout=stream, check=True)
    return PYTHON

def run_script(name, *arguments, density=False):
    work = DENSITY if density else RUN
    log_dir = work/'logs'
    log_dir.mkdir(parents=True, exist_ok=True)
    command = [str(PYTHON), '-u', str(HELPERS/name), *map(str, arguments)]
    environment = dict(os.environ, CLASSIX_NUMERIC_RUN=str(RUN), PYTHONIOENCODING='utf-8', MPLBACKEND='Agg')
    for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMBA_NUM_THREADS'):
        environment[key] = '1'
    with (log_dir/(Path(name).stem+'.log')).open('a', encoding='utf-8') as log:
        process = subprocess.Popen(command, cwd=work, env=environment,
                                   stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                   text=True, encoding='utf-8', errors='replace')
        try:
            for line in process.stdout:
                log.write(line)
                log.flush()
                try:
                    row = json.loads(line)
                except json.JSONDecodeError:
                    row = {}
                if density or row.get('stage') != 'grid' or row.get('grid_fits', 0) % 200 == 0:
                    print(line, end='', flush=True)
            if process.wait():
                raise RuntimeError(f'{name} failed; see {log.name}')
        except KeyboardInterrupt:
            (work/'STOP').touch()
            process.terminate()
            process.wait()
            raise

def prepare_density():
    def density_bytes(relative, payload):
        target = (DENSITY/relative).resolve()
        if not target.is_relative_to(DENSITY.resolve()):
            raise ValueError('Generated files must remain within the run directory.')
        checked_write(target, payload)
    def write_density_source(relative, content):
        density_bytes(relative, content.encode('utf-8'))
    density_policy = json.loads((DENSITY/'protocol.json').read_text(encoding='utf-8'))
    classix_policy = json.loads((RUN/'protocol.json').read_text(encoding='utf-8'))
    if density_policy['datasets'] != classix_policy['datasets']:
        raise RuntimeError('The density and CLASSIX dataset rosters must agree.')
    prior_path = RUN/'results/selected.json'
    prior = json.loads(prior_path.read_text(encoding='utf-8'))
    expected = {(d, m) for d in density_policy['datasets'] for m in ('euclidean', 'manhattan')}
    observed = {(r['dataset'], r['method']) for r in prior}
    if observed != expected or len(prior) != len(expected):
        raise RuntimeError('Complete the CLASSIX screen before preparing density baselines.')

    inputs = {}
    for dataset in density_policy['datasets']:
        for suffix in ('.npz', '.manifest.json'):
            relative = f'data/{dataset}{suffix}'
            payload = (RUN/relative).read_bytes()
            density_bytes(relative, payload)
            if suffix == '.npz':
                inputs[relative] = hashlib.sha256(payload).hexdigest()
    density_bytes('results/classix_selected_prior.json', prior_path.read_bytes())
    grid_path = RUN/'results/grid.jsonl'
    density_bytes('results/classix_grid_prior.jsonl', grid_path.read_bytes())
    prior_assignments = {}
    for row in prior:
        relative_labels = str(row['labels_path']).replace('\\', '/')
        source_labels = RUN/relative_labels
        if not source_labels.resolve().is_relative_to(RUN.resolve()):
            raise ValueError('Prior CLASSIX assignment path escaped this run.')
        target_labels = 'results/classix_assignments_prior/'+Path(relative_labels).name
        if target_labels in prior_assignments:
            raise RuntimeError('Prior assignment filenames must be unique.')
        payload = source_labels.read_bytes()
        density_bytes(target_labels, payload)
        prior_assignments[target_labels] = hashlib.sha256(payload).hexdigest()
    provenance = {
        'prior_root': '..',
        'origin': 'Prepared original features and selected CLASSIX configurations from this notebook run',
        'inputs': inputs,
        'prior_classix_selected_sha256': hashlib.sha256(prior_path.read_bytes()).hexdigest(),
        'prior_classix_grid_sha256': hashlib.sha256(grid_path.read_bytes()).hexdigest(),
        'prior_classix_assignments': prior_assignments,
        'prior_protocol_sha256': hashlib.sha256((RUN/'protocol.json').read_bytes()).hexdigest(),
    }
    write_density_source('input_provenance.json', json.dumps(provenance, indent=2)+'\n')
    for name in ('logs', 'results/plans', 'results/selected', 'results/assignments',
                 'results/graphs', 'plots', 'tables'):
        (DENSITY/name).mkdir(parents=True, exist_ok=True)
    print(f'Prepared identical inputs for {len(inputs)} datasets and {len(prior)} CLASSIX configurations.')

def render_paper_tables():
    subprocess.run([str(PYTHON), str(REPO/'scripts/render_numeric_tables.py'),
                    '--run-dir', str(DENSITY)], check=True)
