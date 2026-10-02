#!/bin/bash

## Submit the cis-genotype test (04.8-cis-genotype-plan.md). One job.
## Needs the 04.4 DSS noSNP prep, the 04.4 BS-SNPer genotypes, and the 04.6 gene region counts (same inputs as 04.7).
##
## Usage:
##   bash 04.8-cis-genotype-submit.sh
## Settings (defaults in brackets; see 04.8-cis-genotype.R):
##   SAMPLE_SETS ["All DropGeno4"]  GENO_K [3]  MAX_DIST [50000]  MIN_CARRIERS [3]  NNULL [20]  FAR_DIST [1000000]
## Resources:
##   OUT_DIR [../../analyses/Haws_04.8-cis-genotype]  ACCOUNT [coenv]  PARTITION [ckpt-all]  MEM [48G]  TIME [6:00:00]

set -euo pipefail

ACCOUNT="${ACCOUNT:-coenv}"
PARTITION="${PARTITION:-ckpt-all}" #Checkpoint nodes, so the shared coenv nodes stay free. --requeue restarts preempted jobs
MEM="${MEM:-48G}"
TIME="${TIME:-6:00:00}"

export CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outDir="$(realpath -m "${OUT_DIR:-${CODE_DIR}/../../analyses/Haws_04.8-cis-genotype}")" #OUT_DIR: run next to another checkout's 04.4/04.6 outputs (e.g. from a git worktree)
mkdir -p "${outDir}/logs"
[ -f "${outDir}/.gitignore" ] || printf "logs/\nrds/\n" > "${outDir}/.gitignore"

# Run from a copy of the R script, so editing it cannot change a queued or running job
export R_SCRIPT="${outDir}/logs/04.8-cis-genotype-$(date +%Y%m%d-%H%M%S).R"
cp "${CODE_DIR}/04.8-cis-genotype.R" "$R_SCRIPT"

id=$(sbatch --parsable --requeue --account="$ACCOUNT" --partition="$PARTITION" --chdir="$outDir" \
  --job-name=04.8-cis --cpus-per-task=2 --mem="$MEM" --time="$TIME" \
  --output=logs/%x_%j.out "${CODE_DIR}/04.8-cis-genotype.job")
echo "04.8-cis: ${id}"
echo "Outputs: ${outDir}"
