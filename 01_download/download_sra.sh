#!/bin/bash
# Step 1: Download SRA files
# Create selected_files.txt with SRR accession numbers first
# Example contents of selected_files.txt:
# SRR8114384
# SRR8114392

prefetch --option-file selected_files.txt
echo "Download complete."

