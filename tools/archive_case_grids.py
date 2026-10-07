#!/usr/bin/env python3
"""Losslessly archive untracked input grids for Git; original local grids remain unchanged."""
import lzma
import hashlib
import json
from pathlib import Path
import subprocess

from restore_case_grids import sha256

ROOT = Path(__file__).resolve().parents[1]


def main():
    tracked = set(subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode('utf-8').split('\0'))
    directory = ROOT / 'cases/grid-archives'
    directory.mkdir(exist_ok=True)
    records = []
    candidates = sorted(p for p in (ROOT / 'cases/manual').rglob('*.cgns')
                        if {'grids', 'grid'} & set(p.relative_to(ROOT).parts)
                        and p.relative_to(ROOT).as_posix() not in tracked)
    for path in candidates:
        digest = sha256(path)
        target = directory / (digest + '.cgns.xz')
        if not target.exists():
            with target.open('xb') as stream, lzma.LZMAFile(stream, mode='wb', preset=6) as out:
                with path.open('rb') as src:
                    for chunk in iter(lambda: src.read(1024 * 1024), b''):
                        out.write(chunk)
        check = hashlib.sha256()
        with lzma.open(target, 'rb') as stream:
            for chunk in iter(lambda: stream.read(1024 * 1024), b''):
                check.update(chunk)
        if check.hexdigest() != digest:
            raise ValueError('Archive content mismatch: ' + str(path))
        records.append({'path': path.relative_to(ROOT).as_posix(), 'bytes': path.stat().st_size,
                        'sha256': digest, 'archive': target.relative_to(ROOT).as_posix(),
                        'archive_sha256': sha256(target)})
    (directory / 'manifest.json').write_text(json.dumps({'schema_version': 1, 'grids': records}, indent=2) + '\n', encoding='utf-8')
    print(f'Archived {len(records)} grid paths using {len({r["archive"] for r in records})} unique XZ files.')


if __name__ == '__main__':
    main()
