import glob
import math
import os
import sys


def rows(path):
    with open(path, encoding="utf-8") as stream:
        return [
            [float(value) for value in line.split()]
            for line in stream
            if line.strip() and not line.startswith("#")
        ]


def require_close(reference, candidate, label, tolerance=1.0e-12):
    if len(reference) != len(candidate):
        raise SystemExit(f"{label}: row count differs")
    for row_index, (left, right) in enumerate(zip(reference, candidate)):
        if len(left) != len(right):
            raise SystemExit(f"{label}: column count differs at row {row_index}")
        for column, (a, b) in enumerate(zip(left, right)):
            scale = max(1.0, abs(a), abs(b))
            if not math.isfinite(a) or not math.isfinite(b) or abs(a - b) > tolerance * scale:
                raise SystemExit(f"{label}: mismatch at ({row_index},{column}): {a} vs {b}")


output = sys.argv[1]
reference_faces = None
reference_loads = None
reference_spanwise = None
for rank in (1, 2, 4):
    face_paths = sorted(glob.glob(os.path.join(output, f"*.boundary.r{rank}.step*.txt")))
    load_paths = glob.glob(os.path.join(output, f"*.loads.r{rank}.txt"))
    span_paths = sorted(glob.glob(os.path.join(output, f"*.spanwise_loads.r{rank}.step*.txt")))
    if len(face_paths) != 2 or len(load_paths) != 1 or len(span_paths) != 2:
        raise SystemExit(
            f"rank {rank}: expected two boundary/spanwise snapshots and one load history, "
            f"got {face_paths}, {span_paths}, {load_paths}"
        )
    face_rows = [rows(path) for path in face_paths]
    span_rows = [rows(path) for path in span_paths]
    load_rows = rows(load_paths[0])
    if reference_faces is None:
        reference_faces = face_rows
        reference_loads = load_rows
        reference_spanwise = span_rows
    else:
        for event, (left, right) in enumerate(zip(reference_faces, face_rows)):
            require_close(left, right, f"rank {rank} boundary event {event}")
        require_close(reference_loads, load_rows, f"rank {rank} loads")
        for event, (left, right) in enumerate(zip(reference_spanwise, span_rows)):
            require_close(left, right, f"rank {rank} spanwise event {event}")

print("boundary 1/2/4-rank equivalence passed")
