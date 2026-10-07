#!/usr/bin/env bash
set -euo pipefail
ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
exec python3 "$ROOT/scripts/run_campaign.py" "$@"
