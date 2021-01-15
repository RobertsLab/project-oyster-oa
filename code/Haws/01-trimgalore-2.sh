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
#SBATCH --chdir=/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2

# I will take the files I already trimmed and rerun TrimGalore and cutadapt. This is because cutadapt removes the most prominent adapter sequence from the file, not all the adapter sequences. There are still adapter sequences in the files that need to be removed. I also need to remove poly-G tails from the samples.

#TrimGalore: Remove adapter sequences
/gscratch/srlab/programs/TrimGalore-0.6.6/trim_galore \
--output_dir /gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2 \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2 \
--threads 28" \
--illumina \
--path_to_cutadapt /usr/lusers/yaaminiv/.local/bin/cutadapt \
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_1_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_1_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_2_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_2_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_3_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_3_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_4_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_4_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_5_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_5_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_6_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_6_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_7_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_7_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_8_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_8_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_9_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_9_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_10_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_10_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_11_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_11_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_12_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_12_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_13_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_13_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_14_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_14_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_15_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_15_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_16_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_16_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_17_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_17_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_18_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_18_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_19_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_19_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_20_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_20_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_21_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_21_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_22_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_22_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_23_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_23_R2_val_2.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_24_R1_val_1.fq.gz	\
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore/zr3644_24_R2_val_2.fq.gz

#MultiQC: Files after second round of adapter trimming
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2/.

#TrimGalore: Remove poly-G tails. Sequence obtained from FastQC reports
/gscratch/srlab/programs/TrimGalore-0.6.6/trim_galore \
--output_dir /gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2/poly-G \
--paired \
--fastqc_args \
"--outdir /gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2/poly-G \
--threads 28" \
--adapter GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG \
--adapter2 GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG \
--path_to_cutadapt /usr/lusers/yaaminiv/.local/bin/cutadapt \
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2/zr3644_1_R1_val_1.fq.gz	\

#MultiQC: Files after poly-G tail trimming
/gscratch/srlab/programs/anaconda3/bin/multiqc \
/gscratch/scrubbed/yaaminiv/Hawes/analyses/trimgalore-2/poly-G/.
