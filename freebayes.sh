#!/bin/bash
#SBATCH --job-name=freebayes
#SBATCH --output=pipe_logs/freebayes_%j.out
#SBATCH --error=pipe_logs/freebayes_%j.err
#SBATCH --mem=5g
#SBATCH --time=00:40:00
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=georgiatiegs@uvic.ca

# Set output directory
OUTPUT_DIR="freebayes_output"

# Create output directory
mkdir -p $OUTPUT_DIR

# Run freebayes for variant calling

# Load modules
module load StdEnv/2023 gcc/12.3 freebayes/1.3.7

# Run freebayes on BAM files
for BAM in align_output_64269857/bam_files/*.sorted.bam; do
    BASENAME=$(basename "$BAM" .sorted.bam)
    freebayes -f WTamuc1412.fa \
              -F 0.03 -p 1 \
              -C 1 \
              --min-base-quality 20 \
              --min-mapping-quality 20 \
              --pooled-continuous \
              "$BAM" > "${OUTPUT_DIR}/${BASENAME}_fb.vcf"
    echo "Variant calling for $BAM complete."
done

echo "Variant calling complete. VCF files are located in ${OUTPUT_DIR}/vcf_files/."
