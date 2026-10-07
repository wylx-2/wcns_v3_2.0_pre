"""Audit this delivery's bounded test outputs and source/configuration/grid snapshot."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import tarfile

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def config(path):
    values = {}
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        line = line.split('#', 1)[0].strip()
        if '=' in line:
            key, value = line.split('=', 1)
            values[key.strip()] = value.strip()
    return values


def rows(path):
    with path.open(encoding='utf-8') as stream:
        return list(csv.DictReader(stream))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('package', type=Path)
    args = parser.parse_args()
    package = args.package.resolve()
    work = ROOT / 'tmp'
    checked_source = 0
    for directory in ['include', 'src', 'third_party/cgns', 'third_party/fftw']:
        originals = {p.relative_to(ROOT) for p in (ROOT / directory).rglob('*') if p.is_file()}
        packaged = {p.relative_to(package / 'source') for p in (package / 'source' / directory).rglob('*') if p.is_file()}
        assert originals == packaged
        for path in originals:
            assert digest(ROOT / path) == digest(package / 'source' / path), path
            checked_source += 1
    assert digest(ROOT / 'CMakeLists.txt') == digest(package / 'source/CMakeLists.txt')
    matrix = json.loads((package / 'campaign.json').read_text())['cases']
    assert [(c['n'], c['riemann']) for c in matrix] == [(n, r) for n in [32, 64, 128] for r in ['roe', 'rusanov', 'roe_all_speed']]
    for item in matrix:
        v = config(package / item['config'])
        expected = {'algorithm.profile': 'scmm6_wcns', 'algorithm.reconstruction': 'mdcd_linear',
                    'algorithm.mdcd.diss': '0.001', 'time.integrator': 'ssprk3', 'run.cfl': '0.3',
                    'turbulence.model': 'none', 'hit.initialization': 'shell_spectrum',
                    'hit.preparation.time': '0.10381382891686522', 'hit.type': 'decay',
                    'run.t_end': '0.31885676024465742', 'hit.seed': '20261003',
                    'hit.n': str(item['n']), 'algorithm.riemann': item['riemann']}
        assert all(v[k] == value for k, value in expected.items()), item['id']
        assert (package / item['config']).parent.joinpath(v['mesh.path']).resolve().is_file()
        assert (package / item['config']).parent.joinpath(v['hit.spectrum_file']).resolve().is_file()
    grids = []
    for n in [16, 32, 64, 128]:
        log = (work / f'cbc9-grid-{n}.log').read_text()
        extents = re.findall(r' cells=(\d+):(\d+):(\d+)', log)
        assert extents == [(str(n // 2),) * 3] * 8
        assert len(re.findall(r'^connectivity ', log, re.M)) == 48
        path = package / f'inputs/grids/hit{n}.cgns'
        assert digest(path) == digest(ROOT / f'cases/manual/case11_hit_decay/grids/hit{n}.cgns')
        grids.append({'n': n, 'blocks': 8, 'cells': n**3, 'sha256': digest(path)})
    smoke = json.loads((work / 'cbc9-smoke/campaign_state.json').read_text())
    assert smoke['mode'] == 'smoke' and len(smoke['cases']) == 9
    max_parseval = max_rematch = 0.0
    for name, item in smoke['cases'].items():
        assert item['status'] == 'complete'
        output = work / 'cbc9-smoke' / item['attempts'][-1]['segment'] / 'output'
        meta = config(output / f'{name}_hit_metadata.txt')
        assert meta['phase'] == 'production' and meta['preparation_complete'] == '1'
        assert meta['FFT'] == 'distributed-slab-MPI+FFTW3'
        for row in rows(output / f'{name}_hit_history.csv'):
            max_parseval = max(max_parseval, abs(float(row['parseval_error'])))
        for row in rows(output / 'preparation' / f'{name}_rematch.csv'):
            max_rematch = max(max_rematch, abs(float(row['E_after']) - float(row['E_target_discrete'])))
    assert max_parseval < 1e-12 and max_rematch < 1e-12
    resumed = json.loads((work / 'cbc9-resume-smoke/campaign_state.json').read_text())['cases']['cbc_n032_roe']
    assert resumed['status'] == 'complete' and len(resumed['attempts']) == 2
    assert [a['mpi_ranks'] for a in resumed['attempts']] == [2, 4]
    assert resumed['attempts'][0]['reason'] == 'wall_time_checkpoint'
    last = rows(work / 'cbc9-resume-smoke' / resumed['attempts'][1]['segment'] / 'output/cbc_n032_roe_hit_history.csv')[-1]
    baseline = rows(work / 'cbc9-smoke/cbc_n032_roe/segment-0001/output/cbc_n032_roe_hit_history.csv')[-1]
    restart_error = max(abs(float(last[k]) - float(baseline[k])) / max(1.0, abs(float(baseline[k])))
                        for k in ['K', 'rho_mean', 'total_energy', 'epsilon_viscous', 'enstrophy'])
    assert restart_error < 1e-10
    preflight = json.loads((work / 'cbc9-preflight/campaign_state.json').read_text())
    actual_dry_runs = [name for name, item in preflight['cases'].items() if item['status'] == 'validated']
    archive = work / 'cbc9-smoke-statistics.tar.gz'
    with tarfile.open(archive) as tar:
        names = tar.getnames()
        assert not any(name.endswith('.cgns') for name in names)
        assert sum(name.endswith('_hit_history.csv') for name in names) == 18
        for line in tar.extractfile('RESULT_CONTENTS.sha256').read().decode().splitlines():
            sha, name = line.split('  ', 1)
            assert hashlib.sha256(tar.extractfile(name).read()).hexdigest() == sha
    report = {'source_files_checked': checked_source + 1, 'production_configs_checked': 9,
              'grids': grids, 'linux_smoke_grid_n': 16, 'linux_smoke_cases_completed': 9,
              'linux_smoke_mpi_ranks': 2, 'max_parseval_error': max_parseval,
              'max_rematch_spectrum_error': max_rematch, 'restart_mpi_ranks': [2, 4],
              'restart_relative_statistics_error': restart_error, 'actual_grid_dry_runs': actual_dry_runs,
              'statistics_archive_checksums_verified': True, 'production_time_evolution_performed': False}
    (package / 'validation.json').write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
