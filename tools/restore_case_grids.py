#!/usr/bin/env python3
"""Restore large case grids from content-addressed XZ files; never replace differing data."""
import argparse
import lzma
import hashlib
import json
from pathlib import Path
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def inside(root, relative):
    path = (root / relative).resolve()
    path.relative_to(root)
    return path


def restore(root, verify_only=False):
    manifest = json.loads((root / 'cases/grid-archives/manifest.json').read_text(encoding='utf-8'))
    for item in manifest['grids']:
        target = inside(root, item['path'])
        compressed = inside(root, item['archive'])
        if sha256(compressed) != item['archive_sha256']:
            raise ValueError('Compressed file checksum mismatch: ' + str(compressed))
        if target.exists():
            if target.stat().st_size != item['bytes'] or sha256(target) != item['sha256']:
                raise ValueError('Existing grid differs; preserve it and resolve manually: ' + str(target))
            print('VERIFIED ' + item['path'], flush=True)
            continue
        if verify_only:
            raise FileNotFoundError(target)
        target.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target.parent, prefix='.restore-grid-', delete=False) as stream:
            temporary = Path(stream.name)
            try:
                with lzma.open(compressed, 'rb') as src:
                    for chunk in iter(lambda: src.read(1024 * 1024), b''):
                        stream.write(chunk)
                stream.close()
                if temporary.stat().st_size != item['bytes'] or sha256(temporary) != item['sha256']:
                    raise ValueError('Restored checksum mismatch: ' + item['path'])
                # Recheck to avoid overwriting a file created while decompression ran.
                if target.exists():
                    raise FileExistsError(target)
                temporary.rename(target)
            finally:
                temporary.unlink(missing_ok=True)
        print('RESTORED ' + item['path'], flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--verify-only', action='store_true')
    args = parser.parse_args()
    restore(args.root.resolve(), args.verify_only)
