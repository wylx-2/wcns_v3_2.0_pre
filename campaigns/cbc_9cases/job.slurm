#!/usr/bin/env bash
# Edit resources/partition/account to match your cluster before sbatch.
#SBATCH --job-name=cbc9
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00
#SBATCH --output=cbc9-slurm-%j.log
set -euo pipefail
# Load the SAME compiler and MPI modules used by build_linux.sh here.
cd "${SLURM_SUBMIT_DIR:?}"
# Solver gets the remaining job budget minus a 300 s checkpoint margin.
# Adjust 85800 when changing --time. Re-submit the same job to resume.
exec bash run_all.sh --np "${SLURM_NTASKS}" --job-seconds 85800
