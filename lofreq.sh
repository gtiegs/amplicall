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

#Filter lofreq vcf to match threshold of other callers

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
    OUTFILE="lofreq_output/${BASENAME}_AF0.03.vcf"

    lofreq filter \
        --af-min 0.03 \
        -i "$VCF" \
        -o "$OUTFILE"
done

echo "LoFreq variant calling complete."
