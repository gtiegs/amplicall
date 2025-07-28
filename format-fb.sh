#!/bin/bash
# Change GT from 0 to 1
for VCF in *fb*norm*.vcf; do
    bcftools +setGT "$VCF" \
        -Ov -o "${VCF%.vcf}_GT.vcf" \
        -- -t a -n c:1
done

# Compress vcfs
for VCF in *fb*norm*GT.vcf; do
    # Compress
    bcftools view "$VCF" -Oz -o "${VCF}.gz"

    # Index
    tabix -p vcf "$VCF.gz"
done