#!/usr/bin/env python3
"""Check audited floor-repair columns for a two-equation RANS run."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("history", type=Path)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("fields", nargs=2)
    args = parser.parse_args()

    lines = [line.strip() for line in args.history.read_text(encoding="utf-8").splitlines()]
    if not lines or not lines[0].startswith("# "):
        raise SystemExit("two-equation history is missing its header")
    header = lines[0][2:].split()
    repair_columns = [f"{field}_floor_repairs" for field in args.fields]
    missing = sorted(set(repair_columns) - set(header))
    if missing:
        raise SystemExit(f"two-equation history is missing columns: {missing}")

    rows = [line.split() for line in lines[1:] if line and not line.startswith("#")]
    if not rows:
        raise SystemExit("two-equation history contains no data rows")
    for row_number, row in enumerate(rows, start=1):
        if len(row) != len(header):
            raise SystemExit(
                f"history row {row_number} has {len(row)} columns; expected {len(header)}"
            )
        for column in repair_columns:
            value = int(row[header.index(column)])
            if value < 0:
                raise SystemExit(f"history column {column} must be non-negative")

    manifest_lines = set(args.manifest.read_text(encoding="utf-8").splitlines())
    for column in repair_columns:
        prefix = f"{column}="
        if not any(line.startswith(prefix) and int(line.removeprefix(prefix)) >= 0
                   for line in manifest_lines):
            raise SystemExit(f"manifest is missing a non-negative {column}")

    print(
        "two-equation floor diagnostics passed: "
        f"rows={len(rows)} fields={','.join(args.fields)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
