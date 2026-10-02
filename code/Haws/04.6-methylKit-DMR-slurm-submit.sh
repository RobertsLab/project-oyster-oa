#!/bin/bash

## Submit the methylKit DMR analysis (04.6-methylKit-DMR-plan.md).
## prep (once) -> regions or segments (once) -> count (once per region set and COV_BASES) -> dmr array (6 tasks per run) -> summary
## Steps that already finished are skipped, so later runs only submit what is new.
##
## Usage:
##   bash 04.6-methylKit-DMR-slurm-submit.sh              #Core grid: REGIONS="tile250 tile1000", observed labels
##   REGIONS=gene bash 04.6-methylKit-DMR-slurm-submit.sh
##   PERM="$(seq 1 20)" bash 04.6-methylKit-DMR-slurm-submit.sh            #20 label permutations of the core grid
##   OVERDISPERSION="MN none" bash 04.6-methylKit-DMR-slurm-submit.sh      #MN plus the negative control
##
## Settings (defaults in brackets; see 04.6-methylKit-DMR-slurm.R). REGIONS, OVERDISPERSION, and PERM take space-separated lists:
##   REGIONS ["tile250 tile1000"]  COV_BASES [3]  LO_COUNT [10]  HI_PERC [99.9]  SAMPLES [All]  OVERDISPERSION [MN]  PERM [0]
## Resources:
##   ACCOUNT [coenv]  PARTITION [ckpt-all]  COUNT_MEM [200G]  DMR_CPUS [16]  DMR_MEM [64G]  DMR_TIME [12:00:00]

set -euo pipefail

ACCOUNT="${ACCOUNT:-coenv}"
PARTITION="${PARTITION:-ckpt-all}" #Checkpoint nodes, so the shared coenv nodes stay free. Jobs use --requeue, so preempted ones restart
COUNT_MEM="${COUNT_MEM:-200G}"
DMR_CPUS="${DMR_CPUS:-16}"
DMR_MEM="${DMR_MEM:-64G}"
DMR_TIME="${DMR_TIME:-12:00:00}"
regionSets=(${REGIONS:-tile250 tile1000})
overdispersionSettings=(${OVERDISPERSION:-MN})
permSettings=(${PERM:-0})
export COV_BASES="${COV_BASES:-3}" LO_COUNT="${LO_COUNT:-10}" HI_PERC="${HI_PERC:-99.9}" SAMPLES="${SAMPLES:-All}"
unset REGIONS OVERDISPERSION PERM #Set for each job below

export CODE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
outDir="$(realpath -m "${CODE_DIR}/../../analyses/Haws_04.6-methylKit-DMR")"
mkdir -p "${outDir}/logs" "${outDir}/downloads"
[ -f "${outDir}/.gitignore" ] || printf "prep-*/\ncount-*/\nsegments-*/*.rds\nDMR-*/rds/\nregions/*.rds\ndownloads/\nlogs/\n" > "${outDir}/.gitignore" #Keep large intermediate files out of git

# RepeatMasker output for the TE region sets. Downloaded here because compute nodes may not have internet access
rmOut="${outDir}/downloads/GCF_902806645.1_cgigas_uk_roslin_v1_rm.out.gz"
[ -s "$rmOut" ] || curl -sSf -o "$rmOut" \
  https://ftp.ncbi.nlm.nih.gov/genomes/all/annotation_releases/29159/102/GCF_902806645.1_cgigas_uk_roslin_v1/GCF_902806645.1_cgigas_uk_roslin_v1_rm.out.gz

common=(--parsable --requeue --account="$ACCOUNT" --partition="$PARTITION" --chdir="$outDir")
job="${CODE_DIR}/04.6-methylKit-DMR-slurm.job"
afterok() { [ $# -gt 0 ] && echo "--dependency=afterok:$(IFS=:; echo "$*")"; return 0; }
# ID of a queued or running job with this name. Lets a second submission wait for prep, regions, segments, or count jobs from
# an earlier submission instead of running them again. Job names include the settings that change their output
queued() { squeue -h -u "$USER" -n "$1" -t PENDING,CONFIGURING,RUNNING,REQUEUED -o %i | head -1; }

# prep: skipped if already finished for this HI_PERC
prepIDs=()
prepName="04.6-prep-hiperc${HI_PERC}"
if [ -f "${outDir}/prep-hiperc${HI_PERC}/prep-complete" ]; then
  echo "prep already complete for HI_PERC=${HI_PERC}, skipping"
elif [ -n "$(queued "$prepName")" ]; then
  prepIDs=($(queued "$prepName"))
  echo "prep: waiting for queued job ${prepIDs[0]}"
else
  prepIDs=($(sbatch "${common[@]}" --job-name="$prepName" --cpus-per-task=4 --mem=120G --time=6:00:00 \
    --output=logs/%x_%j.out "$job" prep))
  echo "prep: ${prepIDs[0]}"
fi

# regions and segments: only submitted if a requested region set needs them
regionIDs=($(queued 04.6-regions))
segmentIDs=($(queued "04.6-segments-hiperc${HI_PERC}"))
for r in "${regionSets[@]}"; do
  if [[ "$r" == segments ]]; then
    if [ ! -f "${outDir}/segments-hiperc${HI_PERC}/segments-complete" ] && [ ${#segmentIDs[@]} -eq 0 ]; then
      segmentIDs=($(sbatch "${common[@]}" $(afterok "${prepIDs[@]}") --job-name="04.6-segments-hiperc${HI_PERC}" --cpus-per-task=4 --mem=200G \
        --time=1-00:00:00 --output=logs/%x_%j.out "$job" segments))
      echo "segments: ${segmentIDs[0]}"
    fi
  elif [[ "$r" != tile* ]]; then
    if [ ! -f "${outDir}/regions/regions-complete" ] && [ ${#regionIDs[@]} -eq 0 ]; then
      regionIDs=($(sbatch "${common[@]}" --job-name=04.6-regions --cpus-per-task=1 --mem=32G --time=2:00:00 \
        --output=logs/%x_%j.out "$job" regions))
      echo "regions: ${regionIDs[0]}"
    fi
  fi
done

# count, then one dmr array for each overdispersion setting and permutation
dmrIDs=()
for r in "${regionSets[@]}"; do
  countIDs=()
  countName="04.6-count-hiperc${HI_PERC}-${r}-cb${COV_BASES}"
  if [ -f "${outDir}/count-hiperc${HI_PERC}-${r}-cb${COV_BASES}/count-complete" ]; then
    echo "count already complete for ${r}, COV_BASES=${COV_BASES}, skipping"
  elif [ -n "$(queued "$countName")" ]; then
    countIDs=($(queued "$countName"))
    echo "count ${r}: waiting for queued job ${countIDs[0]}"
  else
    if [[ "$r" == segments ]]; then upstream=("${prepIDs[@]}" "${segmentIDs[@]}")
    elif [[ "$r" == tile* ]]; then upstream=("${prepIDs[@]}")
    else upstream=("${prepIDs[@]}" "${regionIDs[@]}"); fi
    countIDs=($(REGIONS="$r" sbatch "${common[@]}" $(afterok "${upstream[@]}") --job-name="$countName" --cpus-per-task=4 \
      --mem="$COUNT_MEM" --time=12:00:00 --output=logs/%x_%j.out "$job" count))
    echo "count ${r}: ${countIDs[0]}"
  fi
  for od in "${overdispersionSettings[@]}"; do
    for p in "${permSettings[@]}"; do
      id=$(REGIONS="$r" OVERDISPERSION="$od" PERM="$p" sbatch "${common[@]}" $(afterok "${countIDs[@]}") \
        --job-name=04.6-dmr-${r}-${od}-perm${p} --array=1-6 --cpus-per-task="$DMR_CPUS" --mem="$DMR_MEM" --time="$DMR_TIME" \
        --output=logs/%x_%A_%a.out "$job" dmr)
      dmrIDs+=("$id")
      echo "dmr ${r} ${od} perm ${p}: ${id}"
    done
  done
done

# summary: combines every run in the output directory, so it can also be rerun on its own later
summaryID=$(sbatch "${common[@]}" --dependency=afterany:$(IFS=:; echo "${dmrIDs[*]}") --job-name=04.6-summary \
  --cpus-per-task=1 --mem=8G --time=1:00:00 --output=logs/%x_%j.out "$job" summary)
echo "summary: ${summaryID}"
echo "Outputs: ${outDir}"
