#!/usr/bin/env python3
"""Verify immutable NASA TMR inputs used by the stage-X SA-neg gate."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--cache", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    failures: list[str] = []
    for asset in manifest["assets"]:
        path = args.cache / asset["path"]
        if not path.is_file():
            failures.append(f"missing: {path}")
            continue
        actual = digest(path)
        expected = asset["sha256"].lower()
        if actual != expected:
            failures.append(f"sha256 mismatch: {path} expected={expected} actual={actual}")
        else:
            print(f"verified {asset['path']} {actual}")

    if failures:
        for failure in failures:
            print(f"ERROR {failure}")
        return 1
    print(f"verified {len(manifest['assets'])} immutable TMR assets")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
