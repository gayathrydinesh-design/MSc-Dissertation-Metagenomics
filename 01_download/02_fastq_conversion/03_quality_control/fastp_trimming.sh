#!/bin/bash
# Step 3: Quality control and adapter trimming using fastp

conda activate metagenomics
mkdir -p trimmed qc

for r1 in raw_fastq/*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    base=$(basename $r1 _1.fastq.gz)
    echo "Running fastp on $base..."
    fastp \
        -i $r1 \
        -I $r2 \
        -o trimmed/${base}_1.trim.fastq.gz \
        -O trimmed/${base}_2.trim.fastq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --unqualified_percent_limit 40 \
        --length_required 50 \
        --low_complexity_filter \
        --complexity_threshold 30 \
        -h qc/${base}_fastp.html \
        -j qc/${base}_fastp.json \
        -w 8
    echo "Done: $base"
done
echo "Trimming complete."

