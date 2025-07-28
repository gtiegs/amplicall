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

# Filter VCF files

# Load modules
module load StdEnv/2023 gcc/12.3 bcftools/1.19

# Normalize indels
for VCF in "${OUTPUT_DIR}/"*.vcf; do
    bcftools norm -f WTamuc1412.fa "$VCF" > "${VCF%.vcf}_norm.vcf"
done

# Filter for Q20, Q15, Q10, and Q5
for VCF in "${OUTPUT_DIR}/"*_norm.vcf; do
    bcftools filter -i "QUAL>=20" "$VCF" > "${VCF%.vcf}_Q20.vcf"
    bcftools filter -i "QUAL>=15" "$VCF" > "${VCF%.vcf}_Q15.vcf"
    bcftools filter -i "QUAL>=10" "$VCF" > "${VCF%.vcf}_Q10.vcf"
    bcftools filter -i "QUAL>=5" "$VCF" > "${VCF%.vcf}_Q5.vcf"
done

echo "Filtering complete. Filtered VCF files are located in ${OUTPUT_DIR}/vcf_files/."
