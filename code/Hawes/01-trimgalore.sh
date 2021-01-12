#!/bin/bash
## Job Name
#SBATCH --job-name=hawes-fastqc
## Allocation Definition
#SBATCH --account=srlab
#SBATCH --partition=srlab
## Resources
## Nodes
#SBATCH --nodes=1
## Walltime (days-hours:minutes:seconds format)
#SBATCH --time=3-00:00:00
## Memory per node
#SBATCH --mem=200G
##turn on e-mail notification
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yaaminiv@uw.edu
## Specify the working directory for this job
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore

#TrimGalore
/gscratch/srlab/programs/TrimGalore-0.6.6/trim_galore \
--output_dir /gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore \
--threads 28" \
--illumina \
--clip_R1 10 \
--clip_R2 10 \
--three_prime_clip_R1 10 \
--three_prime_clip_R2 10 \
--path_to_cutadapt /usr/lusers/yaaminiv/.local/bin/cutadapt \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_1_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_1_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_2_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_2_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_3_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_3_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_4_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_4_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_5_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_5_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_6_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_6_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_7_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_7_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_8_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_8_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_9_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_9_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_10_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_10_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_11_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_11_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_12_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_12_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_13_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_13_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_14_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_14_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_15_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_15_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_16_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_16_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_17_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_17_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_18_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_18_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_19_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_19_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_20_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_20_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_21_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_21_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_22_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_22_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_23_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_23_R2.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_24_R1.fq.gz \
/gscratch/scrubbed/yaaminiv/Hawes/data/zr3644_24_R2.fq.gz

#MultiQC
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/.
