#!/bin/bash
# Step 9: Taxonomic classification of MAGs using GTDB-Tk

# Download reference data (run once only)
# wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxiliary_files/gtdbtk_r214_data.tar.gz
# tar -xvzf gtdbtk_r214_data.tar.gz

# Convert .fa to .fna for GTDB-Tk compatibility
echo "Converting .fa to .fna..."
cd HQ_MAGs
for f in *.fa; do
    mv "$f" "${f%.fa}.fna"
done
cd ..

echo "Running GTDB-Tk classification..."
gtdbtk classify_wf \
    --genome_dir HQ_MAGs \
    --out_dir gtdbtk_out \
    --cpus 8 \
    --extension fna

echo "Taxonomic classification complete. Results in gtdbtk_out/"

