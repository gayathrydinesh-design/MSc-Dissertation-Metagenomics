#!/bin/bash
# Step 6: AMR gene annotation using RGI against CARD database
# Usage: bash rgi_AMR.sh SAMPLE_NAME

conda activate rgi

SAMPLE=$1
echo "Running RGI on $SAMPLE..."

rgi main \
    -i ${SAMPLE}_metaspades/contigs.fasta \
    --output_file rgi_${SAMPLE}

echo "AMR annotation complete: $SAMPLE"
echo "Output: rgi_${SAMPLE}.txt"

