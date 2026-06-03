#!/bin/bash
# Step 10: Microbial community profiling using MetaPhlAn
# Installation: conda create -n metaphlan_env -c conda-forge -c bioconda metaphlan -y

conda activate metaphlan_env
mkdir -p metaphlan_output

for r1 in host_removed/*/*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    base=$(basename $r1 _1.fastq.gz)
    echo "Running MetaPhlAn on $base..."
    metaphlan \
        ${r1},${r2} \
        --input_type fastq \
        --nproc 16 \
        --mapout metaphlan_output/${base}_map.txt \
        -o metaphlan_output/${base}_profile.txt
    echo "Done: $base"
done

echo "Community profiling complete. Profiles in metaphlan_output/"

