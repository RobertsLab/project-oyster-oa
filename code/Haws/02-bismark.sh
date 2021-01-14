#!/bin/bash
## Job Name
#SBATCH --job-name=hawes-bismark
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=20-00:00:00
## Memory per node
#SBATCH --mem=500G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Hawes/analyses/bismark

set -e #Stop script if any command fails

# Directories and programs
bismark_dir="/gscratch/srlab/programs/Bismark-0.22.3/"
bowtie2_dir="/gscratch/srlab/programs/bowtie2-2.4.1-linux-x86_64/"
samtools="/gscratch/srlab/programs/samtools-1.10/samtools"

source /gscratch/srlab/programs/scripts/paths.sh

#Genome Preparation
#Sam White prepared the bisulfite genome: https://robertslab.github.io/sams-notebook/2019/02/21/Data-Management-Create-C.gigas-Bisulfite-Genome-with-Bismark-on-Mox.html

genome_folder="/gscratch/scrubbed/yaaminiv/Hawes/data/Cg-genome/Crassostrea_gigas.oyster_v9.dna_sm.toplevel/"

#Align to Bisulfite Genome

reads_dir="/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/"

find ${reads_dir}*_R1_val_1.fq.gz \
| xargs basename -s _R1_val_1.fq.gz | xargs -I{} ${bismark_dir}/bismark \
--path_to_bowtie ${bowtie2_dir} \
-genome ${genome_folder} \
-p 4 \
-score_min L,0,-0.6 \
--non_directional \
-1 ${reads_dir}{}_R1_val_1.fq.gz \
-2 ${reads_dir}{}_R2_val_2.fq.gz

#Deduplication

find *.bam | \
xargs basename -s .bam | \
xargs -I{} ${bismark_dir}/deduplicate_bismark \
--bam \
--paired \
{}.bam

#Methylation Extraction

${bismark_dir}/bismark_methylation_extractor \
--bedGraph --counts --scaffolds \
--multicore 14 \
--buffer_size 75% \
*deduplicated.bam

#Sort for Downstream Applications

find *deduplicated.bam | \
xargs basename -s .bam | \
xargs -I{} ${samtools} \
sort --threads 28 {}.bam \
-o {}.sorted.bam

#Index for Downstream Applications

find *.sorted.bam | \
xargs basename -s .sorted.bam | \
xargs -I{} ${samtools} \
index -@ 28 {}.sorted.bam

#HTML Processing Report

${bismark_dir}/bismark2report

#Summary Report

${bismark_dir}/bismark2summary

#Merge CpG Coverage Information

find *deduplicated.bismark.cov.gz \
| xargs basename -s deduplicated.bismark.cov.gz \
| xargs -I{} ${bismark_dir}/coverage2cytosine \
--genome_folder ${genome_folder} \
-o {} \
--merge_CpG \
--zero_based \
{}deduplicated.bismark.cov.gz

#Create Merged bedgraphs

for f in *merged_CpG_evidence.cov
do
  STEM=$(basename "${f}" .CpG_report.merged_CpG_evidence.cov)
  cat "${f}" | awk -F $'\t' 'BEGIN {OFS = FS} {if ($5+$6 >= 5) {print $1, $2, $3, $4}}' \
  > "${STEM}"_5x.bedgraph
done
