#!/bin/bash
for VCF in *lf*norm*.vcf.gz; do
  base="${VCF%.vcf.gz}"
  tmp="tmp_${base}.vcf.gz"
  out="${base}_GT.vcf.gz"

  # Create a header file with the FORMAT line
  echo '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">' > header.txt
  
  # Add the header line properly using annotate
  bcftools annotate -h header.txt "$VCF" -Oz -o "$tmp"

  # Index temp file
  tabix -p vcf "$tmp"
 
  # This adds a sample called "SAMPLE1" with GT=1 
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t./.\n' "$tmp" | \
  bcftools view -h "$tmp" | \
  sed 's/#CHROM.*/#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE/' > "${tmp%.gz}"
  
  # Add the body
  bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO\tGT\t./.\n' "$VCF" >> "${tmp%.gz}"
  
  # Compress and index
  bgzip "${tmp%.gz}"
  tabix -p vcf "$tmp"

  # Now use setGT to set actual genotypes based on your criteria
  bcftools +setGT "$tmp" -Oz -o "$out" -- -t a -n c:1

  # Index the output
  tabix -p vcf "$out"

  rm "$tmp"
done