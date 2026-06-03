#!/bin/bash
# Step 8: MAG quality assessment using CheckM
# Usage: bash checkm_quality.sh SAMPLE_NAME

SAMPLE=$1
conda activate checkm

echo "Running CheckM on ${SAMPLE} bins..."
checkm lineage_wf ${SAMPLE}_bins checkm_out_${SAMPLE} -t 8 -x fa

checkm qa checkm_out_${SAMPLE}/lineage.ms checkm_out_${SAMPLE} \
    -f checkm_quality_summary_${SAMPLE}.tsv \
    --tab_table

# Filter high quality MAGs: completeness >= 80%, contamination <= 10%
awk -F'\t' '{print $1"\t"$12"\t"$13"\t"($12-5*$13)}' \
    checkm_quality_summary_${SAMPLE}.tsv | \
    awk -F'\t' '$2 >= 80 && $3 <= 10 && $4 >= 50' \
    > HQ_MAGs_${SAMPLE}.tsv

# Copy HQ MAGs
mkdir -p HQ_MAGs
for k in $(awk '{print $1}' HQ_MAGs_${SAMPLE}.tsv); do
    cp ${SAMPLE}_bins/${k}.fa HQ_MAGs/
done

echo "CheckM complete. HQ MAGs saved to HQ_MAGs/"

