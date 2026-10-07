"""Stage the CURRENT source worktree and prepared CBC matrix; finalize after validation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
ASSETS = ROOT / 'campaigns/cbc_9cases'
NAME = 'WCNS_v2.6_CBC_9cases_20261007'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def copy_tree(src, dst):
    for path in sorted(src.rglob('*')):
        if path.is_file() and '__pycache__' not in path.parts and path.suffix not in ['.pyc', '.exe']:
            target = dst / path.relative_to(src)
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, target)


def manifest(root):
    files = sorted(p for p in root.rglob('*') if p.is_file() and
                   p.name != 'PACKAGE_CONTENTS.sha256' and '__pycache__' not in p.parts)
    (root / 'PACKAGE_CONTENTS.sha256').write_text(''.join(digest(p) + '  ' + p.relative_to(root).as_posix() + '\n'
                                                         for p in files), encoding='utf-8')
    return len(files)


def stage(root):
    if root.exists():
        raise SystemExit('Refusing existing directory: ' + str(root))
    root.mkdir(parents=True)
    copy_tree(ASSETS, root)
    source = root / 'source'
    for directory in ['include', 'src', 'third_party/cgns', 'third_party/fftw']:
        copy_tree(ROOT / directory, source / directory)
    for path in [ROOT / n for n in ['CMakeLists.txt', 'LICENSE.md', 'THIRD_PARTY_NOTICES.md', '算法补充.md']] + list((ROOT / 'tools').glob('*.cpp')):
        target = source / path.relative_to(ROOT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
    (source / 'README.md').write_text('WCNS v2.6 current source snapshot, 2026-10-07.\n\n'
                                    'Build through ../build_linux.sh. See ../README.md for the CBC campaign.\n'
                                    'MPI, FFTW and CGNS must remain enabled. No production run was performed locally.\n', encoding='utf-8')
    (source / 'docs').mkdir()
    shutil.copyfile(ROOT / 'docs/all-speed-roe.md', source / 'docs/all-speed-roe.md')
    for name in ['channel_retau180_from_t250_all_speed_roe.wcns', 'full_case_template_v2.3.wcns']:
        target = source / 'examples' / name
        target.parent.mkdir(exist_ok=True)
        shutil.copyfile(ROOT / 'examples' / name, target)
    source_sha = hashlib.sha256(''.join(p.relative_to(source).as_posix() + ':' + digest(p) + '\n'
                                       for p in sorted(source.rglob('*')) if p.is_file()).encode()).hexdigest()
    head = subprocess.check_output(['git', 'rev-parse', '--short=12', 'HEAD'], cwd=ROOT, text=True).strip()
    revision = head + '-worktree-cbc9-' + source_sha[:16]
    (source / 'WCNS_SOURCE_REVISION').write_text(revision + '\n', encoding='utf-8')
    for n in [16, 32, 64, 128]:
        dest = root / 'inputs/grids' / f'hit{n}.cgns'
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / 'cases/manual/case11_hit_decay/grids' / dest.name, dest)
    copy_tree(ROOT / 'cases/manual/case11_hit_decay/reference', root / 'inputs/reference')
    values = {}
    for line in (ROOT / 'cases/manual/case11_hit_decay/prepared_production.wcns').read_text(encoding='utf-8-sig').splitlines():
        line = line.split('#', 1)[0].strip()
        if line and '=' in line:
            key, value = line.split('=', 1)
            values[key.strip()] = value.strip()
    for key in ['hit.forcing', 'hit.forcing_power', 'hit.forcing_kmax', 'hit.remove_mean_acceleration']:
        values.pop(key, None)
    values.update({'algorithm.reconstruction': 'mdcd_linear', 'run.cfl': '0.3', 'run.max_wall_time': '0',
                   'hit.spectrum_file': '../inputs/reference/cbc_spectrum_42_nondimensional.dat',
                   'output.field.every_time': '0',
                   'output.field.explicit_times': values['hit.sample_times'],
                   'output.checkpoint.every_time': '0.05',
                   'output.checkpoint.explicit_times': values['hit.sample_times']})
    cases = []
    (root / 'cases').mkdir()
    for n in [32, 64, 128]:
        for solver in ['roe', 'rusanov', 'roe_all_speed']:
            name = f'cbc_n{n:03d}_{solver}'
            config = dict(values)
            config.update({'case.name': name, 'mesh.path': f'../inputs/grids/hit{n}.cgns',
                           'hit.n': str(n), 'algorithm.riemann': solver,
                           'output.directory': 'manual_output/' + name})
            if solver != 'roe_all_speed':
                config = {k: v for k, v in config.items() if not k.startswith('algorithm.roe_all_speed.')}
            (root / 'cases' / (name + '.wcns')).write_text('# Server production configuration; use run_all.sh for sequential/resumable execution.\n' +
                                                        ''.join(f'{k} = {v}\n' for k, v in config.items()), encoding='utf-8')
            cases.append({'id': name, 'n': n, 'riemann': solver, 'config': f'cases/{name}.wcns'})
    (root / 'campaign.json').write_text(json.dumps({'version': '2.6', 'package_date': '2026-10-07',
                                                    'source_revision': revision, 'source_snapshot_sha256': source_sha,
                                                    'cases': cases}, indent=2) + '\n', encoding='utf-8')
    manifest(root)
    print(json.dumps({'staged': str(root), 'source_revision': revision}, indent=2))


def finalize(root):
    if not (root / 'validation.txt').is_file():
        raise SystemExit('Write validation.txt after bounded validation, then finalize.')
    if any((root / p).exists() for p in ['build', 'results', 'results_smoke', 'results_dry-run']):
        raise SystemExit('Build/results belong outside the upload tree during validation.')
    archive = root.with_name(root.name + '.tar.gz')
    if archive.exists():
        raise SystemExit('Refusing existing archive: ' + str(archive))
    count = manifest(root)
    files = [(root / line.split('  ', 1)[1]) for line in (root / 'PACKAGE_CONTENTS.sha256').read_text(encoding='utf-8').splitlines()]
    files.append(root / 'PACKAGE_CONTENTS.sha256')
    with archive.open('xb') as raw, tarfile.open(fileobj=raw, mode='w:gz') as tar:
        for path in files:
            info = tar.gettarinfo(str(path), arcname=root.name + '/' + path.relative_to(root).as_posix())
            info.mode = 0o755 if path.suffix == '.sh' else 0o644
            info.uid = info.gid = 0
            info.uname = info.gname = ''
            with path.open('rb') as stream:
                tar.addfile(info, stream)
    Path(str(archive) + '.sha256').write_text(digest(archive) + '  ' + archive.name + '\n', encoding='utf-8')
    print(json.dumps({'archive': str(archive), 'bytes': archive.stat().st_size, 'sha256': digest(archive), 'files': count}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--stage', type=Path)
    group.add_argument('--finalize', type=Path)
    args = parser.parse_args()
    if args.stage:
        stage(args.stage.resolve())
    else:
        finalize(args.finalize.resolve())
