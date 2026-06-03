#!/bin/bash
# Step 4: Host read removal using KneadData

conda activate kneaddata_env
mkdir -p host_removed logs

# Download human genome database (run once only)
# kneaddata_database --download human_genome bowtie2 kneaddata_db

for r1 in trimmed/*_1.trim.fastq.gz; do
    r2=${r1/_1.trim.fastq.gz/_2.trim.fastq.gz}
    base=$(basename $r1 _1.trim.fastq.gz)
    echo "Running KneadData on $base..."
    mkdir -p host_removed/${base}
    kneaddata \
        -i1 $r1 \
        -i2 $r2 \
        -db kneaddata_db \
        -o host_removed/${base} \
        --output-prefix ${base} \
        -t 8 \
        --remove-intermediate-output \
        --log logs/${base}_kneaddata.log
    echo "Done: $base"
done
echo "Host removal complete."

