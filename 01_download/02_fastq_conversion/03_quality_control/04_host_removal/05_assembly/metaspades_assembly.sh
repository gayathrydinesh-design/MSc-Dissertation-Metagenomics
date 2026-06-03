#!/bin/bash
# Step 5: Metagenomic assembly using MetaSPAdes

conda activate metagenomics

# Select best 5 samples by longest contig
echo "Selecting best 5 samples by longest contig..."
for d in *_metaspades; do
    max=$(awk 'BEGIN{m=0;l=0} /^>/ {if (l>m) m=l; l=0; next} {l+=length($0)} END {if (l>m) m=l; print m}' "$d/contigs.fasta")
    echo "$d $max"
done | sort -k2 -nr | head -5

# Run assembly — pass sample name as argument: bash metaspades_assembly.sh SAMPLE_NAME
SAMPLE=$1
echo "Assembling $SAMPLE..."
metaspades.py \
    -1 ${SAMPLE}_1.fastq.gz \
    -2 ${SAMPLE}_2.fastq.gz \
    -o ${SAMPLE}_metaspades \
    --threads 16 \
    --meta
echo "Assembly complete: $SAMPLE"

