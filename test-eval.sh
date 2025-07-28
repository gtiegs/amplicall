#!/bin/bash
# Configuration
BASE_OUT_DIR="vcfeval_output"
mkdir -p "$BASE_OUT_DIR"

# Freebayes loop
for VCF in all_vcfs/*fb*norm*.vcf.gz; do
    filename=$(basename "$VCF" .vcf.gz)
    rtg vcfeval \
    -t WTamuc1412.sdf \
    --baseline truth.vcf.gz \
    --calls "$VCF" \
    --all-records \
    --squash-ploidy \
    --threads 4 \
    --output "${BASE_OUT_DIR}/${filename}"
done