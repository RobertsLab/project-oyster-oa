#!/bin/bash

## Submit the variance partitioning (04.7-variance-partitioning-plan.md). One job.
## Needs the 04.4 DSS noSNP prep, the 04.4 BS-SNPer genotypes, and the 04.6 gene region counts.
##
## Usage:
##   bash 04.7-variance-partitioning-submit.sh
## Settings (defaults in brackets; see 04.7-variance-partitioning.R):
##   SAMPLE_SETS ["All DropGeno4"]  GENO_K [3]  NPERM [9999]  NNULL [100]
## Resources:
##   OUT_DIR [../../analyses/Haws_04.7-variance-partitioning]  ACCOUNT [coenv]  PARTITION [ckpt-all]  MEM [48G]  TIME [6:00:00]

set -euo pipefail

ACCOUNT="${ACCOUNT:-coenv}"
PARTITION="${PARTITION:-ckpt-all}" #Checkpoint nodes, so the shared coenv nodes stay free. --requeue restarts preempted jobs
MEM="${MEM:-48G}"
TIME="${TIME:-6:00:00}"

export CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outDir="$(realpath -m "${OUT_DIR:-${CODE_DIR}/../../analyses/Haws_04.7-variance-partitioning}")" #OUT_DIR: run next to another checkout's 04.4/04.6 outputs (e.g. from a git worktree)
mkdir -p "${outDir}/logs"
[ -f "${outDir}/.gitignore" ] || printf "logs/\n" > "${outDir}/.gitignore"

# Run from a copy of the R script, so editing it cannot change a queued or running job
export R_SCRIPT="${outDir}/logs/04.7-variance-partitioning-$(date +%Y%m%d-%H%M%S).R"
cp "${CODE_DIR}/04.7-variance-partitioning.R" "$R_SCRIPT"

id=$(sbatch --parsable --requeue --account="$ACCOUNT" --partition="$PARTITION" --chdir="$outDir" \
  --job-name=04.7-varpart --cpus-per-task=2 --mem="$MEM" --time="$TIME" \
  --output=logs/%x_%j.out "${CODE_DIR}/04.7-variance-partitioning.job")
echo "04.7-varpart: ${id}"
echo "Outputs: ${outDir}"
