#!/usr/bin/env python3
"""Build a solver-only v2.6 source package from a Git commit (or an explicit working-tree rehearsal)."""
import argparse
import gzip
import hashlib
import io
import json
from pathlib import Path
import subprocess
import tarfile

ROOT = Path(__file__).resolve().parents[1]
PREFIXES = ('include/', 'src/', 'third_party/cgns/', 'third_party/fftw/')
EXPLICIT = {
    'CMakeLists.txt': 'CMakeLists.txt', 'LICENSE.md': 'LICENSE.md',
    'THIRD_PARTY_NOTICES.md': 'THIRD_PARTY_NOTICES.md', '算法补充.md': '算法补充.md',
    'docs/minimal-build-readme.md': 'README.md', 'docs/minimal-runtime.md': 'docs/runtime.md',
    'examples/full_case_template_v2.3.wcns': 'examples/full_case_template_v2.3.wcns',
    'examples/channel_retau180_from_t250_all_speed_roe.wcns': 'examples/channel_retau180_from_t250_all_speed_roe.wcns',
}
BUILD = '''#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cmake -S "$ROOT" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release \\
  -DWCNS_ENABLE_MPI=ON -DWCNS_ENABLE_CGNS=ON -DWCNS_ENABLE_FFTW=ON \\
  -DWCNS_BUILD_TESTS=OFF -DWCNS_BUILD_TOOLS=OFF "$@"
cmake --build "$ROOT/build" --target wcns_run --parallel "${JOBS:-2}"
'''


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def create(args):
    destination = args.output.resolve()
    archive = destination.with_name(destination.name + '.tar.gz')
    if destination.exists() or archive.exists():
        raise SystemExit('Refusing an existing destination or archive.')
    revision = git('rev-parse', args.ref + '^{commit}').decode().strip()
    if args.working_tree:
        names = [p.relative_to(ROOT).as_posix() for prefix in PREFIXES
                 for p in (ROOT / prefix).rglob('*') if p.is_file() and '__pycache__' not in p.parts]
    else:
        names = [n for n in git('ls-tree', '-r', '--name-only', revision).decode('utf-8').splitlines()
                 if n.startswith(PREFIXES)]
    mapping = {name: name for name in names}
    mapping.update(EXPLICIT)
    payload = {}
    for original, target in sorted(mapping.items()):
        data = (ROOT / original).read_bytes() if args.working_tree else git('show', revision + ':' + original)
        if Path(target).suffix in ['.cpp', '.hpp', '.md', '.txt', '.wcns'] or target == 'CMakeLists.txt':
            data = data.replace(b'\r\n', b'\n')
        payload[target] = data
    snapshot = hashlib.sha256(''.join(name + ':' + hashlib.sha256(data).hexdigest() + '\n'
                                     for name, data in sorted(payload.items())).encode()).hexdigest()
    identity = revision + ('-worktree-' + snapshot[:16] if args.working_tree else '')
    payload['WCNS_SOURCE_REVISION'] = (identity + '\n').encode()
    payload['build_linux.sh'] = BUILD.encode()
    payload['package-info.json'] = (json.dumps({'program_version': '2.6', 'git_commit': revision,
                                               'working_tree': args.working_tree,
                                               'source_snapshot_sha256': snapshot,
                                               'build_tools': False, 'build_tests': False}, indent=2) + '\n').encode()
    payload['PACKAGE_CONTENTS.sha256'] = ''.join(hashlib.sha256(data).hexdigest() + '  ' + name + '\n'
                                                for name, data in sorted(payload.items())).encode()
    destination.mkdir(parents=True)
    for name, data in payload.items():
        target = destination / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
    with archive.open('xb') as raw, gzip.GzipFile(fileobj=raw, mode='wb', filename='', mtime=0) as gz, \
            tarfile.open(fileobj=gz, mode='w') as tar:
        for name, data in sorted(payload.items()):
            info = tarfile.TarInfo(destination.name + '/' + name)
            info.size = len(data)
            info.mode = 0o755 if name.endswith('.sh') else 0o644
            tar.addfile(info, io.BytesIO(data))
    sha = hashlib.sha256(archive.read_bytes()).hexdigest()
    Path(str(archive) + '.sha256').write_text(sha + '  ' + archive.name + '\n', encoding='utf-8')
    print(json.dumps({'directory': str(destination), 'archive': str(archive), 'files': len(payload),
                      'bytes': archive.stat().st_size, 'sha256': sha, 'source_revision': identity}, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--ref', default='HEAD')
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--working-tree', action='store_true', help='Explicit rehearsal; final delivery uses a committed ref')
    create(parser.parse_args())
