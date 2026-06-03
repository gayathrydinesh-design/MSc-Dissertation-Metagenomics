#!/bin/bash
# Step 11: Genome coverage estimation using CoverM
# Installation: conda create -n coverm_env -c bioconda -c conda-forge coverm -y

conda activate coverm_env
mkdir -p coverage_output

for r1 in host_removed/*/*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    base=$(basename $r1 _1.fastq.gz)
    echo "Running CoverM on $base..."
    coverm genome \
        --coupled ${r1} ${r2} \
        --genome-fasta-directory HQ_MAGs \
        --output-file coverage_output/${base}_coverage.tsv \
        --threads 16
    echo "Done: $base"
done

echo "Coverage estimation complete. Files in coverage_output/"

