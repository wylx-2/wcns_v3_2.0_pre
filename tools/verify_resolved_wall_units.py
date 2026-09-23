#!/usr/bin/env python3
"""Verify resolved-wall unit columns against the exported viscous traction."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def read_table(path: Path) -> tuple[list[str], list[list[float]]]:
    header: list[str] | None = None
    rows: list[list[float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if line.startswith("# patch_index "):
            header = line[2:].split()
        elif line and not line.startswith("#"):
            rows.append([float(value) for value in line.split()])
    if header is None or not rows or any(len(row) != len(header) for row in rows):
        raise ValueError(f"invalid or empty boundary table: {path}")
    return header, rows


def expected_class(y_plus: float) -> int:
    if y_plus <= 5.0:
        return 0
    if y_plus < 30.0:
        return 1
    if y_plus <= 300.0:
        return 2
    return 3


def relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--boundary", type=Path, required=True)
    parser.add_argument("--reynolds", type=float, required=True)
    parser.add_argument("--relative-tolerance", type=float, default=2.0e-12)
    args = parser.parse_args()
    if not math.isfinite(args.reynolds) or args.reynolds <= 0.0:
        raise ValueError("--reynolds must be finite and positive")

    header, values = read_table(args.boundary)
    required = {
        "nx", "ny", "nz", "mu_w", "wall_distance", "friction_velocity",
        "wall_y_plus", "wall_y_plus_class", "viscous_traction_x",
        "viscous_traction_y", "viscous_traction_z",
    }
    missing = sorted(required.difference(header))
    if missing:
        raise ValueError("boundary table is missing columns: " + ",".join(missing))

    maximum_density_identity_error = 0.0
    maximum_normal_error = 0.0
    minimum_inferred_density = math.inf
    maximum_inferred_density = 0.0
    classification_mismatches = 0
    for values_row in values:
        row = dict(zip(header, values_row, strict=True))
        normal = [row["nx"], row["ny"], row["nz"]]
        traction = [
            row["viscous_traction_x"],
            row["viscous_traction_y"],
            row["viscous_traction_z"],
        ]
        normal_norm = sum(value * value for value in normal)
        maximum_normal_error = max(maximum_normal_error, abs(normal_norm - 1.0))
        normal_traction = sum(left * right for left, right in zip(traction, normal))
        tangential = [
            value - normal_traction * normal_component
            for value, normal_component in zip(traction, normal)
        ]
        shear = math.sqrt(sum(value * value for value in tangential))
        mu_w = row["mu_w"]
        distance = row["wall_distance"]
        friction_velocity = row["friction_velocity"]
        y_plus = row["wall_y_plus"]
        if not all(math.isfinite(value) for value in values_row):
            raise ValueError("boundary wall-unit row contains a non-finite value")
        if mu_w <= 0.0 or distance <= 0.0 or friction_velocity <= 0.0 or y_plus <= 0.0:
            raise ValueError("resolved-wall unit row must be strictly positive")
        density_from_stress = shear / (friction_velocity * friction_velocity)
        density_from_y_plus = (
            y_plus * mu_w / (friction_velocity * distance * args.reynolds)
        )
        maximum_density_identity_error = max(
            maximum_density_identity_error,
            relative_error(density_from_stress, density_from_y_plus),
        )
        minimum_inferred_density = min(minimum_inferred_density, density_from_stress)
        maximum_inferred_density = max(maximum_inferred_density, density_from_stress)
        if row["wall_y_plus_class"] != float(expected_class(y_plus)):
            classification_mismatches += 1

    checks = {
        "unit_normals": maximum_normal_error <= args.relative_tolerance,
        "stress_y_plus_identity": (
            maximum_density_identity_error <= args.relative_tolerance
        ),
        "classification": classification_mismatches == 0,
    }
    report = {
        "boundary": str(args.boundary),
        "face_count": len(values),
        "maximum_normal_error": maximum_normal_error,
        "maximum_density_identity_relative_error": maximum_density_identity_error,
        "minimum_inferred_density": minimum_inferred_density,
        "maximum_inferred_density": maximum_inferred_density,
        "classification_mismatches": classification_mismatches,
        "checks": checks,
        "passed": all(checks.values()),
    }
    print(json.dumps(report, indent=2, sort_keys=True))
    return 0 if report["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
