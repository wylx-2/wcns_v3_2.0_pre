#!/usr/bin/env python3
"""Check that unsteady SA history rows preserve the fixed model-column schema."""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("history", type=Path)
    args = parser.parse_args()

    lines = [line.strip() for line in args.history.read_text(encoding="utf-8").splitlines()]
    if not lines or not lines[0].startswith("# "):
        raise SystemExit("SA history is missing its header")
    header = lines[0][2:].split()
    required = {
        "nu_tilde_l2",
        "nu_tilde_reference_l2",
        "nu_tilde_normalized_l2",
        "nu_tilde_linf",
        "nu_tilde_reference_linf",
        "nu_tilde_normalized_linf",
    }
    if not required.issubset(header):
        raise SystemExit(f"SA history is missing columns: {sorted(required - set(header))}")

    rows = [line.split() for line in lines[1:] if line and not line.startswith("#")]
    if not rows:
        raise SystemExit("SA history contains no data rows")
    for row_number, row in enumerate(rows, start=1):
        if len(row) != len(header):
            raise SystemExit(
                f"SA history row {row_number} has {len(row)} columns; expected {len(header)}"
            )
        for name in ("nu_tilde_reference_l2", "nu_tilde_reference_linf"):
            if not math.isnan(float(row[header.index(name)])):
                raise SystemExit(f"unsteady SA history column {name} must be NaN")

    print(f"SA unsteady history schema passed: rows={len(rows)} columns={len(header)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
