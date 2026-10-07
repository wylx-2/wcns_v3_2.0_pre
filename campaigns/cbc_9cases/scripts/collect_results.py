"""Collect every continuation segment; default excludes large fields/checkpoints."""
import argparse
import datetime
import fcntl
import io
import json
from pathlib import Path
import tarfile
from verify_package import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT / 'results')
    parser.add_argument('--full', action='store_true', help='Include CGNS fields/checkpoints; potentially many GB')
    parser.add_argument('--archive', type=Path)
    args = parser.parse_args()
    results = args.results.resolve()
    if not (results / 'campaign_state.json').is_file():
        raise SystemExit('No campaign_state.json in ' + str(results))
    lock = (results / '.campaign.lock').open('a')
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        raise SystemExit('Campaign is running; collect after it stops to get a consistent archive.')
    stamp = datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    archive = (args.archive or ROOT / f'CBC9_results_{stamp}{"_full" if args.full else "_statistics"}.tar.gz').resolve()
    if archive.is_relative_to(results):
        raise SystemExit('Place the archive outside the results directory.')
    files = {}
    for path in sorted(results.rglob('*')):
        if path.is_file() and path.name != '.campaign.lock' and (args.full or path.suffix.lower() != '.cgns'):
            files['results/' + path.relative_to(results).as_posix()] = path
    for directory in ['inputs/reference', 'build-logs', 'cases']:
        for path in sorted((ROOT / directory).rglob('*')):
            if path.is_file():
                files['provenance/' + path.relative_to(ROOT).as_posix()] = path
    for name in ['campaign.json', 'README.md', 'PACKAGE_CONTENTS.sha256', 'source/WCNS_SOURCE_REVISION', 'validation.txt']:
        if (ROOT / name).exists():
            files['provenance/' + name] = ROOT / name
    digest_lines = ''.join(sha256(path) + '  ' + name + '\n' for name, path in files.items())
    with archive.open('xb') as raw, tarfile.open(fileobj=raw, mode='w:gz') as tar:
        for name, path in files.items():
            tar.add(path, arcname=name, recursive=False)
        content = digest_lines.encode()
        info = tarfile.TarInfo('RESULT_CONTENTS.sha256')
        info.size = len(content)
        tar.addfile(info, io.BytesIO(content))
    sidecar = Path(str(archive) + '.sha256')
    sidecar.write_text(sha256(archive) + '  ' + archive.name + '\n', encoding='utf-8')
    state = json.loads((results / 'campaign_state.json').read_text(encoding='utf-8'))
    print(json.dumps({'archive': str(archive), 'sha256': sha256(archive), 'full': args.full,
                      'mode': state['mode'], 'status': {k: v['status'] for k, v in state['cases'].items()}}, indent=2))


if __name__ == '__main__':
    main()
