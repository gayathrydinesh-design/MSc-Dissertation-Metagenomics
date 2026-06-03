#!/bin/bash
# Step 7: MAG binning using MetaBAT2
# Usage: bash MAG_binning.sh SAMPLE_NAME

SAMPLE=$1
echo "Mapping reads to contigs for $SAMPLE..."

# Map reads to assembled contigs
bwa mem -t 24 ${SAMPLE}_metaspades/contigs.fasta \
    ${SAMPLE}_1.fastq.gz \
    ${SAMPLE}_2.fastq.gz | \
    samtools view -@ 4 -bS - | \
    samtools sort -@ 8 -o ${SAMPLE}_sorted.bam

# Index sorted BAM file
samtools index ${SAMPLE}_sorted.bam

# Generate contig depth/coverage
conda activate Basic_protocol_MAG1
jgi_summarize_bam_contig_depths \
    --outputDepth depth_${SAMPLE}.txt \
    ${SAMPLE}_sorted.bam

# Bin contigs with MetaBAT2
mkdir -p ${SAMPLE}_bins
metabat2 \
    -i ${SAMPLE}_metaspades/contigs.fasta \
    -a depth_${SAMPLE}.txt \
    -o ${SAMPLE}_bins/bin \
    -m 2000 \
    --minContigDepth 2

echo "MAG binning complete: $SAMPLE"

