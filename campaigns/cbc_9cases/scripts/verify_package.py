"""Check immutable release files; build/results created afterwards are allowed."""
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    count = 0
    for line in (ROOT / 'PACKAGE_CONTENTS.sha256').read_text(encoding='utf-8').splitlines():
        digest, name = line.split('  ', 1)
        path = (ROOT / name).resolve()
        path.relative_to(ROOT)
        if sha256(path) != digest:
            raise SystemExit('CHECKSUM MISMATCH: ' + name)
        count += 1
    print(f'Verified {count} packaged files (SHA256).', flush=True)


if __name__ == '__main__':
    main()
