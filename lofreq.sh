#!/bin/bash

# Add indel quality scores
for BAM in bam_files/*.sorted.bam; do
    BASENAME=$(basename "$BAM") 
    OUTBAM="bam_files/indelqual_${BASENAME}"

    lofreq indelqual --dindel -f WTamuc1412.fa -o "$OUTFILE" "$BAM"
    samtools index "$OUTFILE"
done

# Run lofreq
lofreq faidx WTamuc1412.fa

# Run single‐threaded lofreq
for BAM in bam_files/indelqual_*.sorted.bam; do
    BASENAME=$(basename "$BAM" .sorted.bam)
    OUTVCF="lofreq_output/${BASENAME}_lf.vcf"

    lofreq call \
      --call-indels \
      -f WTamuc1412.fa \
      -o "$OUTVCF" \
      "$BAM"
done

#Filter lofreq vcf based on AF

# Pure cultures
for VCF in lofreq_output/*FG*.vcf; do
    BASENAME=$(basename "$VCF" .vcf)
    OUTFILE="lofreq_output/${BASENAME}_AF0.9.vcf"

    lofreq filter \
        --af-min 0.90 \
        -i "$VCF" \
        -o "$OUTFILE"
done





# Pooled cultures

for VCF in lofreq_output1/*.vcf; do
    BASENAME=$(basename "$VCF" .vcf)
    OUTFILE="lofreq_output1/${BASENAME}_AF0.03.vcf"

    lofreq filter \
        --af-min 0.03 \
        -i "$VCF" \
        -o "$OUTFILE"
done

lofreq filter \
    --af-min 0.03 \
    -i lofreq_output/indelqual_3P6VV8_5_5_pooled_lf.vcf \
    -o lofreq_output/indelqual_3P6VV8_5_5_pooled_lf_AF0.03.vcf

lofreq filter \
    --af-min 0.07 \
    -i lofreq_output/indelqual_3P6VV8_4_10_pooled_lf.vcf \
    -o lofreq_output/indelqual_3P6VV8_4_10_pooled_lf_AF0.07.vcf

lofreq filter \
    --af-min 0.12 \
    -i lofreq_output/indelqual_3P6VV8_3_15_pooled_lf.vcf\
    -o lofreq_output/indelqual_3P6VV8_3_15_pooled_lf_AF0.12.vcf

lofreq filter \
    --af-min 0.17 \
    -i lofreq_output/indelqual_3P6VV8_2_20_pooled_lf.vcf  \
    -o lofreq_output/indelqual_3P6VV8_2_20_pooled_lf_AF0.17.vcf

echo "LoFreq variant calling complete."
