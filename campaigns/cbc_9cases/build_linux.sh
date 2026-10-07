#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
JOBS=${JOBS:-2}
mkdir -p "$ROOT/build-logs"
python3 "$ROOT/scripts/verify_package.py"
cmake -S "$ROOT/source" -B "$ROOT/build" -DCMAKE_BUILD_TYPE=Release \
  -DWCNS_ENABLE_MPI=ON -DWCNS_ENABLE_FFTW=ON -DWCNS_ENABLE_CGNS=ON \
  -DWCNS_BUILD_TESTS=OFF "$@" 2>&1 | tee "$ROOT/build-logs/configure.log"
cmake --build "$ROOT/build" --target wcns_run --parallel "$JOBS" 2>&1 | tee "$ROOT/build-logs/build.log"
