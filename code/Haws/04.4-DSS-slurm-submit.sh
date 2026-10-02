#!/bin/bash

## Submit the DSS analysis (04.4-DSS-plan.md).
## prep-raw (once) -> qc (once) and prep (once per filter setting) -> fit array (NCHUNK tasks) -> summary
##
## Usage:
##   bash 04.4-DSS-slurm-submit.sh qc      #prep-raw if needed, then the sample QC table only
##   bash 04.4-DSS-slurm-submit.sh         #DSS run for the current settings (runs prep-raw and prep if needed)
##   bash 04.4-DSS-slurm-submit.sh all     #qc and the DSS run, sharing one prep-raw job
##   bash 04.4-DSS-slurm-submit.sh compare #Sample-set sensitivity table (needs All, Drop3H2, DropPC2, OutRM summaries)
##   bash 04.4-DSS-slurm-submit.sh permsummary #Observed vs permuted DML counts (needs PERM = 0 and PERM > 0 runs)
##   bash 04.4-DSS-slurm-submit.sh snp-prep    #Mark CpGs with BS-SNPer SNPs (VCFs in analyses/Haws_04.4-DSS/prep-bssnper/)
##   bash 04.4-DSS-slurm-submit.sh snp-summary #SNP enrichment among DML; genotype PCA
##
## Settings (defaults in brackets; see 04.4-DSS-slurm.R):
##   LO_COV [5]  HI_PERC [99.9]  PRESENCE [cell5]  MIN_METH [10]  SNP_FILTER [none]  SAMPLES [All]  MODEL [interaction]  PERM [0]  NCHUNK [1]
## Resources:
##   ACCOUNT [coenv]  PARTITION [cpu-g2]  FIT_MEM [64G]  FIT_TIME [4:00:00]
## Examples:
##   for s in Drop3H2 DropPC2 OutRM; do SAMPLES=$s bash 04.4-DSS-slurm-submit.sh; done   #Sample-set sensitivity checks
##   MODEL=additive bash 04.4-DSS-slurm-submit.sh
##   for s in $(seq 1 20); do PERM=$s bash 04.4-DSS-slurm-submit.sh; done; bash 04.4-DSS-slurm-submit.sh permsummary

set -euo pipefail

target="${1:-dss}"
case "$target" in qc|dss|all|compare|permsummary|snp-prep|snp-summary) ;; *) echo "Unknown target: $target (use qc, dss, all, compare, permsummary, snp-prep, or snp-summary)" >&2; exit 1 ;; esac

ACCOUNT="${ACCOUNT:-coenv}"
PARTITION="${PARTITION:-cpu-g2}"
FIT_MEM="${FIT_MEM:-64G}"
FIT_TIME="${FIT_TIME:-4:00:00}"
export LO_COV="${LO_COV:-5}" HI_PERC="${HI_PERC:-99.9}" PRESENCE="${PRESENCE:-cell5}" MIN_METH="${MIN_METH:-10}" SNP_FILTER="${SNP_FILTER:-none}" SAMPLES="${SAMPLES:-All}"
export MODEL="${MODEL:-interaction}" PERM="${PERM:-0}" NCHUNK="${NCHUNK:-1}"

export CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outDir="$(realpath -m "${CODE_DIR}/../../analyses/Haws_04.4-DSS")"
mkdir -p "${outDir}/logs"
[ -f "${outDir}/.gitignore" ] || printf "prep-*/\nDSS-*/chunks/\nDSS-*/rds/\nDSS-*-perm*/\nlogs/\n" > "${outDir}/.gitignore" #Keep large intermediate files out of git

# Run jobs from a snapshot of the R script, so editing 04.4-DSS-slurm.R cannot change jobs already queued or running
# (Rscript reads the file as it goes)
export R_SCRIPT="${outDir}/logs/04.4-DSS-slurm-$(date +%Y%m%d-%H%M%S).R"
cp "${CODE_DIR}/04.4-DSS-slurm.R" "$R_SCRIPT"

common=(--parsable --requeue --account="$ACCOUNT" --partition="$PARTITION" --chdir="$outDir")
job="${CODE_DIR}/04.4-DSS-slurm.job"

if [ "$target" = "compare" ] || [ "$target" = "permsummary" ] || [ "$target" = "snp-prep" ] || [ "$target" = "snp-summary" ]; then
  id=$(sbatch "${common[@]}" --job-name=04.4-${target} --cpus-per-task=2 --mem=150G --time=3:00:00 \
    --output=logs/%x_%j.out "$job" "$target")
  echo "${target}: ${id}"
  exit 0
fi

# prep-raw: skipped if already finished
rawDependency=()
if [ -f "${outDir}/prep-raw/prep-complete" ]; then
  echo "prep-raw already complete, skipping"
else
  rawID=$(sbatch "${common[@]}" --job-name=04.4-prep-raw --cpus-per-task=4 --mem=120G --time=6:00:00 \
    --output=logs/%x_%j.out "$job" prep-raw)
  rawDependency=(--dependency=afterok:${rawID})
  echo "prep-raw: ${rawID}"
fi

if [ "$target" != "dss" ]; then
  qcID=$(sbatch "${common[@]}" "${rawDependency[@]}" --job-name=04.4-qc --cpus-per-task=2 --mem=120G --time=4:00:00 \
    --output=logs/%x_%j.out "$job" qc)
  echo "qc: ${qcID} (outputs: ${outDir}/sample-QC)"
  [ "$target" = "qc" ] && exit 0
fi

# prep for this filter setting: skipped if already finished
setting="cov${LO_COV}-hiperc${HI_PERC}-${PRESENCE}$([ "$MIN_METH" = "none" ] || echo "-meth${MIN_METH}")$([ "$SNP_FILTER" = "none" ] || echo "-noSNP")-${SAMPLES}"
prepDependency=("${rawDependency[@]}")
if [ -f "${outDir}/prep-${setting}/prep-complete" ]; then
  echo "prep already complete for ${setting}, skipping"
else
  prepID=$(sbatch "${common[@]}" "${rawDependency[@]}" --job-name=04.4-prep --cpus-per-task=2 --mem=120G --time=4:00:00 \
    --output=logs/%x_%j.out "$job" prep)
  prepDependency=(--dependency=afterok:${prepID})
  echo "prep ${setting}: ${prepID}"
fi

run="${setting}-${MODEL}$([ "$PERM" = "0" ] || echo "-perm${PERM}")"
fitID=$(sbatch "${common[@]}" "${prepDependency[@]}" --job-name=04.4-fit --array=1-${NCHUNK} \
  --cpus-per-task=1 --mem="$FIT_MEM" --time="$FIT_TIME" \
  --output=logs/%x_%A_%a.out "$job" fit)
summaryID=$(sbatch "${common[@]}" --dependency=afterok:${fitID} --job-name=04.4-summary \
  --cpus-per-task=1 --mem=64G --time=2:00:00 \
  --output=logs/%x_%j.out "$job" summary)
echo "DSS-${run}: fit array ${fitID}, summary ${summaryID}"
echo "Outputs: ${outDir}/DSS-${run}"
