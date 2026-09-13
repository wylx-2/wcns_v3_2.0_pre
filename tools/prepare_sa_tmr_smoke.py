#!/usr/bin/env python3
"""Materialize a stage-X TMR smoke config from verified local assets."""

from __future__ import annotations

import argparse
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--output-directory", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    if not args.grid.is_file():
        raise SystemExit(f"verified TMR grid is missing: {args.grid}")
    text = args.template.read_text(encoding="utf-8")
    replacements = {
        "@SA_TMR_GRID@": args.grid.resolve().as_posix(),
        "@SA_TMR_OUTPUT_DIRECTORY@": args.output_directory.resolve().as_posix(),
    }
    for marker, value in replacements.items():
        if text.count(marker) != 1:
            raise SystemExit(f"template marker must occur exactly once: {marker}")
        text = text.replace(marker, value)
    args.config.parent.mkdir(parents=True, exist_ok=True)
    args.config.write_text(text, encoding="utf-8", newline="\n")
    print(args.config.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
