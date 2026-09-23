#!/usr/bin/env python3
"""Compare continuous and restarted lightweight LES time statistics."""

import math
import pathlib
import sys


def read_statistics(path: pathlib.Path):
    metadata = {}
    rows = {}
    covariances = {}
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith("# accepted_events "):
            metadata["accepted_events"] = int(line.split()[2])
        elif line.startswith("# first_step "):
            parts = line.split()
            metadata["first_step"] = int(parts[2])
            metadata["first_time"] = float(parts[4])
        elif line.startswith("# last_step "):
            parts = line.split()
            metadata["last_step"] = int(parts[2])
            metadata["last_time"] = float(parts[4])
        elif line.startswith("covariance "):
            parts = line.split()
            covariances[(parts[1], parts[2])] = [float(value) for value in parts[3:]]
        elif not line.startswith("#"):
            parts = line.split()
            rows[parts[0]] = [float(value) for value in parts[1:]]
    return metadata, rows, covariances


def close(lhs: float, rhs: float) -> bool:
    # Cross-rank restart intentionally mixes deterministic serial checkpoint
    # state with a differently ordered MPI spatial reduction.  The integral
    # means retain the tight relative tolerance; near-zero RMS values need an
    # absolute floor for that last-bit reduction difference.
    return math.isclose(lhs, rhs, rel_tol=2.0e-13, abs_tol=1.0e-12)


def main() -> int:
    if len(sys.argv) < 3:
        raise SystemExit("usage: check_les_time_statistics.py CONTINUOUS RESTART [REFERENCE ...]")
    paths = [pathlib.Path(argument) for argument in sys.argv[1:]]
    datasets = [read_statistics(path) for path in paths]
    reference_metadata, reference_rows, reference_covariances = datasets[0]
    if reference_metadata.get("accepted_events") != 2:
        raise RuntimeError("continuous run did not record exactly two accepted physical steps")
    if reference_metadata.get("first_step") != 1 or reference_metadata.get("last_step") != 2:
        raise RuntimeError("continuous run recorded the wrong accepted-step interval")
    if set(reference_rows) != {"total_mass", "total_energy"}:
        raise RuntimeError("time-statistics quantity identity differs")
    if set(reference_covariances) != {("total_mass", "total_energy")}:
        raise RuntimeError("time-statistics covariance identity differs")
    for name, values in reference_rows.items():
        if int(values[0]) != 2 or not close(values[1], 0.002):
            raise RuntimeError(f"{name} has the wrong count or accepted-dt weight")

    for path, (metadata, rows, covariances) in zip(paths[1:], datasets[1:]):
        if (
            metadata != reference_metadata
            or set(rows) != set(reference_rows)
            or set(covariances) != set(reference_covariances)
        ):
            raise RuntimeError(f"metadata or quantity identity differs: {path}")
        for name, expected in reference_rows.items():
            actual = rows[name]
            if len(actual) != len(expected) or any(
                not close(lhs, rhs) for lhs, rhs in zip(actual, expected)
            ):
                raise RuntimeError(f"time statistics differ for {name}: {path}")
        for pair, expected in reference_covariances.items():
            actual = covariances[pair]
            if len(actual) != len(expected) or any(
                not close(lhs, rhs) for lhs, rhs in zip(actual, expected)
            ):
                raise RuntimeError(f"time-statistics covariance differs for {pair}: {path}")
    print(f"validated LES time-statistics restart continuity across {len(paths)} runs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
