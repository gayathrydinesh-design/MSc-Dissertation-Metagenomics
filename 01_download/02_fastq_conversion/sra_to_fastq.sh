#!/bin/bash
# Step 2: Convert .sra files to FASTQ format

mkdir -p fastq_files

while read srr; do
    echo "Converting $srr..."
    fasterq-dump "$srr" -O fastq_files
done < selected_files.txt

echo "Compressing FASTQ files..."
gzip fastq_files/*.fastq
echo "FASTQ conversion complete."

