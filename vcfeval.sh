#!/bin/bash

# Compress truth VCF
bcftools view truth.vcf -Oz -o truth.vcf.gz
tabix -p vcf truth.vcf.gz

# Create reference
rtg format --output WTamuc1412.sdf WTamuc1412.fa

# Run vcfeval on 6_pooled_20_c3_norm.vcf
rtg vcfeval \
  -t WTamuc1412.sdf \
  --baseline truth.vcf.gz \
  --calls all_vcfs/6_pooled_20_c3_norm.vcf.gz \
  --all-records \
  --squash-ploidy \
  --threads 4 \
  --output c3_eval_output/


# Run vcfeval on all sample/caller/filter combinations

# Configuration
BASE_OUT_DIR="vcfeval_output"
mkdir -p "$BASE_OUT_DIR"

# lofreq loop
for VCF in all_vcfs/*lf*norm*GT.vcf.gz; do
  filename=${VCF##*/indelqual_}
  filename=${filename%.vcf.gz}
  rtg vcfeval \
    -t WTamuc1412.sdf \
    --baseline truth.vcf.gz \
    --calls "$VCF" \
    --all-records \
    --squash-ploidy \
    --threads 4 \
    --output "${BASE_OUT_DIR}/${filename}"
done

# Clair3 loop
for VCF in all_vcfs/*c3*norm*.vcf.gz; do
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
