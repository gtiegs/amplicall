#!/bin/bash
#SBATCH --job-name=align
#SBATCH --output=pipe_logs/align_%j.out
#SBATCH --error=pipe_logs/align_%j.err
#SBATCH --mem=10g
#SBATCH --time=01:00:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=georgiatiegs@uvic.ca

# Set output directory
OUTPUT_DIR="align_output_${SLURM_JOBID}"

# Create necessary directories
mkdir -p "$OUTPUT_DIR"/{fastq_filtered,fastq_plots,bam_files,bam_plots}

# Step 1: Filter FASTQ files
echo "Starting to filter FASTQ files"

# Load modules
module load StdEnv/2023 gcc/12.3 python/3.10

# Activate nanofilt virtual environment
source /home/gtiegs/projects/def-ld172/gtiegs/nanofilt/nanofilt_env/bin/activate

# Run nanofilt on FASTQ files
for FASTQ in *.fastq; do
    BASENAME=$(basename "$FASTQ" .fastq)
    NanoFilt -l 900 --maxlength 1300 "$FASTQ" -q 20 > "${OUTPUT_DIR}/fastq_filtered/${BASENAME}_filtered.fastq"
done
deactivate

echo "FASTQ files filtered successfully. Filtered files are located in ${OUTPUT_DIR}/fastq_filtered/."

# Step 2: Generate read quality plots
echo "Starting to generate nanoplot report for filtered FASTQ files."

# Load modules
module load StdEnv/2023 gcc/12.3 python/3.10
module load gcc arrow/17.0.0
module load scipy-stack/2024b numpy/2.1.1

# Activate nanoplot virtual environment
source /home/gtiegs/projects/def-ld172/gtiegs/nanoplot/nanoplot_env/bin/activate

# Run nanoplot on FASTQ files
for FASTQ in ${OUTPUT_DIR}/fastq_filtered/*_filtered.fastq; do
    BASENAME=$(basename "$FASTQ" _filtered.fastq)
    NanoPlot --store --fastq "$FASTQ" \
             -o "${OUTPUT_DIR}/fastq_plots/${BASENAME}" \
             -c blueviolet -cm Earth -f png
    echo "Nanoplot for $FASTQ complete."
done
deactivate

echo "Nanoplot report generated successfully. Plots located in ${OUTPUT_DIR}/fastq_plots/."

# Step 3: Run minimap2 for alignment
echo "Starting to run minimap2 for alignment"

# Load modules
module load StdEnv/2023 gcc/12.3 minimap2

# Run minimap2 on filtered FASTQ files
for FASTQ in ${OUTPUT_DIR}/fastq_filtered/*_filtered.fastq; do
    BASENAME=$(basename "$FASTQ" _filtered.fastq)
    minimap2 -ax map-ont WTamuc1412.fa "$FASTQ" \
             > "${OUTPUT_DIR}/bam_files/${BASENAME}.sam"
    echo "$FASTQ aligned successfully"
done

echo "FASTQ files aligned successfully. Sam files are located in ${OUTPUT_DIR}/bam_files/."

# Step 4: Convert SAM files to BAM files
echo "Starting to convert SAM files to BAM files"

# Load modules
module load StdEnv/2023 gcc/12.3 samtools/1.20

# Convert SAM files to BAM files
for SAM in ${OUTPUT_DIR}/bam_files/*.sam; do
    BASENAME=$(basename "$SAM" .sam)
    samtools view -bS "$SAM" | \
    samtools sort -o "${OUTPUT_DIR}/bam_files/${BASENAME}.sorted.bam"
    samtools index "${OUTPUT_DIR}/bam_files/${BASENAME}.sorted.bam"
    rm "$SAM"  # Remove SAM file after conversion
done

echo "SAM files converted to BAM files successfully. BAM and BAI files are located in ${OUTPUT_DIR}/bam_files/."

# Step 5: Generate alignment quality plots
echo "Starting to generate alignment quality plots"

# Load modules
module load StdEnv/2023 gcc/12.3 python/3.10
module load arrow/17.0.0
module load scipy-stack/2024b numpy/2.1.1

# Activate nanoplot virtual environment
source /home/gtiegs/projects/def-ld172/gtiegs/nanoplot/nanoplot_env/bin/activate

# Run nanoplot on BAM files
for BAM in ${OUTPUT_DIR}/bam_files/*.sorted.bam; do
    BASENAME=$(basename "$BAM" .sorted.bam)
    NanoPlot --store --bam "$BAM" \
             -o "${OUTPUT_DIR}/bam_plots/${BASENAME}" \
             -c forestgreen -cm Earth -f png
done
deactivate

echo "Alignment quality plots generated successfully. Plots are located in ${OUTPUT_DIR}/bam_plots/."
