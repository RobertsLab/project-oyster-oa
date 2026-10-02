#!/bin/bash

## Submit the Slurm version of 04.2 methylKit parameter testing.
## prep (once) -> dml array (12 tasks per overdispersion setting) -> summary
##
## Usage:
##   bash 04.3-methylKit-slurm-submit.sh                #overdispersion = "none" (matches 04.2)
##   bash 04.3-methylKit-slurm-submit.sh none MN        #Run both side by side
##
## Optional environment variables (defaults in brackets):
##   ACCOUNT [coenv]  PARTITION [cpu-g2]  DML_CPUS [32]  DML_MEM [150G]  DML_TIME [1-00:00:00]
##   HI_PERC [99.9]   Upper coverage percentile filter; "none" reproduces 04.2 as written (no upper filter)
## Example: HI_PERC=none PARTITION=ckpt bash 04.3-methylKit-slurm-submit.sh none MN

set -euo pipefail

ACCOUNT="${ACCOUNT:-coenv}"
PARTITION="${PARTITION:-cpu-g2}"
DML_CPUS="${DML_CPUS:-32}"
DML_MEM="${DML_MEM:-150G}"
DML_TIME="${DML_TIME:-1-00:00:00}"
export HI_PERC="${HI_PERC:-99.9}"

overdispersionSettings=("$@")
[ ${#overdispersionSettings[@]} -eq 0 ] && overdispersionSettings=(none)
for od in "${overdispersionSettings[@]}"; do
  case "$od" in none|MN|shrinkMN) ;; *) echo "Unknown overdispersion setting: $od (use none, MN, or shrinkMN)" >&2; exit 1 ;; esac
done

export CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outDir="$(realpath -m "${CODE_DIR}/../../analyses/Haws_04.3-methylKit-slurm")"
mkdir -p "${outDir}/logs"
[ -f "${outDir}/.gitignore" ] || printf "prep-*/\nDML-*/rds/\nlogs/\n" > "${outDir}/.gitignore" #Keep large intermediate files out of git

common=(--parsable --requeue --account="$ACCOUNT" --partition="$PARTITION" --chdir="$outDir")

# Stage 1: prep, skipped if it already finished for this HI_PERC setting
dependency=()
if [ -f "${outDir}/prep-hiperc-${HI_PERC}/prep-complete" ]; then
  echo "Prep already complete for HI_PERC=${HI_PERC}, skipping"
else
  prepID=$(sbatch "${common[@]}" --job-name=04.3-prep --cpus-per-task=4 --mem=200G --time=12:00:00 \
    --output=logs/%x_%j.out "${CODE_DIR}/04.3-methylKit-slurm.job" prep)
  dependency=(--dependency=afterok:${prepID})
  echo "prep: ${prepID}"
fi

# Stages 2 and 3 for each overdispersion setting
for od in "${overdispersionSettings[@]}"; do
  dmlID=$(sbatch "${common[@]}" "${dependency[@]}" --job-name=04.3-dml-${od} --array=1-12 \
    --cpus-per-task="$DML_CPUS" --mem="$DML_MEM" --time="$DML_TIME" \
    --output=logs/%x_%A_%a.out "${CODE_DIR}/04.3-methylKit-slurm.job" dml "$od")
  summaryID=$(sbatch "${common[@]}" --dependency=afterok:${dmlID} --job-name=04.3-summary-${od} \
    --cpus-per-task=1 --mem=8G --time=1:00:00 \
    --output=logs/%x_%j.out "${CODE_DIR}/04.3-methylKit-slurm.job" summary "$od")
  echo "overdispersion=${od}: dml array ${dmlID}, summary ${summaryID}"
done

echo "Outputs: ${outDir}"
