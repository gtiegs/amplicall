#!/bin/bash
#SBATCH --job-name=vcffilter
#SBATCH --output=vcffilter_%j.out
#SBATCH --error=vcffilter_%j.err
#SBATCH --mem=3g
#SBATCH --time=00:05:00

# Filter VCF files

# Load modules
module load StdEnv/2023 gcc/12.3 bcftools/1.19

# Normalize indels
for VCF in *.vcf; do
    bcftools norm -f WTamuc1412.fa "$VCF" > "${VCF%.vcf}_norm.vcf"
done

# Filter for Q20, Q15, Q10, and Q5
for VCF in *_norm.vcf; do
    bcftools filter -i "QUAL>=20" "$VCF" > "${VCF%.vcf}_Q20.vcf"
    bcftools filter -i "QUAL>=15" "$VCF" > "${VCF%.vcf}_Q15.vcf"
    bcftools filter -i "QUAL>=10" "$VCF" > "${VCF%.vcf}_Q10.vcf"
    bcftools filter -i "QUAL>=5" "$VCF" > "${VCF%.vcf}_Q5.vcf"
done

echo "Filtering complete. Filtered VCF files are located in ${OUTPUT_DIR}/vcf_files/."
