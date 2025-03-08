#!/bin/bash
# Set directories and files
reads_dir="/gscratch/srlab/sr320/github/project-oyster-oa/data/Haws-11/"
genome_folder="/gscratch/srlab/sr320/github/project-oyster-oa/data/Haws-11/"
output_dir="."
checkpoint_file="completed_samples.log"
read_suffix="_R1_val_1_val_1_val_1.fq.gz"
read_suffix2="_R2_val_2_val_2_val_2.fq.gz"

# Create the checkpoint file if it doesn't exist
touch ${checkpoint_file}

# Use SLURM_ARRAY_TASK_ID to pick one sample
files=(${reads_dir}*${read_suffix})
file="${files[$SLURM_ARRAY_TASK_ID]}"
sample_name=$(basename "$file" "${read_suffix}")

# Skip if sample already processed
if grep -q "^${sample_name}$" ${checkpoint_file}; then
    echo "Sample ${sample_name} already processed. Skipping..."
    exit 0
fi

# Create a directory to hold chunks
chunk_dir="${sample_name}_chunks"
mkdir -p ${chunk_dir}

# Decompress and split the FASTQ files into smaller parts,
# appending a .fq suffix so that each chunk is valid FASTQ.
zcat ${reads_dir}${sample_name}${read_suffix} | split -l 4000000 --additional-suffix=.fq - ${chunk_dir}/${sample_name}_R1_chunk_
zcat ${reads_dir}${sample_name}${read_suffix2} | split -l 4000000 --additional-suffix=.fq - ${chunk_dir}/${sample_name}_R2_chunk_

# Process each chunk in a loop
for chunk in ${chunk_dir}/${sample_name}_R1_chunk_*.fq; do
    # Remove the .fq extension and extract the chunk identifier (e.g., aa, ab, etc.)
    base=$(basename "$chunk" .fq)
    chunk_id=${base##*_}
    echo "Processing chunk ${chunk_id} for sample ${sample_name}..."
    
    # Define corresponding R2 chunk file (it will also have the .fq suffix)
    r2_chunk="${chunk_dir}/${sample_name}_R2_chunk_${chunk_id}.fq"
    
    bismark \
        -genome ${genome_folder} \
        -p 12 \
        -score_min L,0,-0.8 \
        -1 ${chunk} \
        -2 ${r2_chunk} \
        --non_directional \
        -o ${output_dir} \
        --basename ${sample_name}_chunk${chunk_id} \
        2> "${sample_name}_chunk${chunk_id}-bismark_summary.txt"
done

# Merge the individual chunk BAM files into one file
merged_bam="${sample_name}_merged.bam"
sorted_bam="${sample_name}_merged.sorted.bam"
chunk_bams=(${output_dir}/${sample_name}_chunk*_bismark_bt2*.bam)

echo "Looking for chunk BAM files with pattern: ${output_dir}/${sample_name}_chunk*_bismark_bt2*.bam"
echo "Found chunk BAM files: ${chunk_bams[@]}"

if [ ${#chunk_bams[@]} -gt 0 ] && [ -e "${chunk_bams[0]}" ]; then
    echo "Merging chunk BAM files for sample ${sample_name}..."
    samtools merge ${merged_bam} "${chunk_bams[@]}"
    echo "Sorting merged BAM file..."
    samtools sort -o ${sorted_bam} ${merged_bam}
    echo "Indexing sorted BAM file..."
    samtools index ${sorted_bam}
    echo "Merged, sorted, and indexed BAM saved as ${sorted_bam}."
else
    echo "No chunk BAM files found to merge. Please check Bismark outputs for errors."
fi

# Mark the sample as processed
echo ${sample_name} >> ${checkpoint_file}
echo "Sample ${sample_name} processing complete."