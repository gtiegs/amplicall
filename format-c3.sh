#!/bin/bash
# Compress vcfs
for VCF in *c3*norm*.vcf; do
    # Compress
    bcftools view "$VCF" -Oz -o "${VCF}.gz"

    # Index
    tabix -p vcf "$VCF.gz"
done
