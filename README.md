# MSc Dissertation — Metagenomic Pipeline for AMR Gene Reservoir Analysis

**Exploring Antimicrobial Resistance Gene Reservoirs in the Gut Microbiome of Type 1 Diabetes Patients through Metagenomics**

*Gayathry Dinesh — M.Sc. Genomic Science, Central University of Kerala (2026)*

---

## Overview

This repository contains the complete bioinformatics pipeline developed for my M.Sc. dissertation. The project investigated the distribution and diversity of antimicrobial resistance (AMR) gene-harboring microbial taxa in the gut microbiome of Type 1 Diabetes (T1D) patients using shotgun metagenomics.

The pipeline covers the full workflow from raw sequencing data download to taxonomic profiling, AMR gene annotation, metagenome-assembled genome (MAG) reconstruction, and coverage analysis.

---

This repository contains all scripts used for downloading, preprocessing, assembly, AMR annotation, MAG reconstruction, taxonomic classification, coverage estimation, and downstream visualization of shotgun metagenomic datasets.



## Pipeline Overview

```
Raw SRA Data
    ↓
FASTQ Conversion (fasterq-dump)
    ↓
Quality Control and Adapter Trimming (fastp)
    ↓
Host Read Removal (KneadData + Bowtie2)
    ↓
Metagenomic Assembly (MetaSPAdes)
    ↓
AMR Gene Annotation (RGI against CARD database)
    ↓
MAG Binning (MetaBAT2 + MaxBin2)
    ↓
MAG Quality Assessment (CheckM)
    ↓
Taxonomic Classification (GTDB-Tk)
    ↓
Community Profiling (MetaPhlAn)
    ↓
Coverage Estimation (CoverM)
```

---

## Tools Used

| Step | Tool | Version |
|------|------|---------|
| Data download | SRA Toolkit (prefetch, fasterq-dump) | - |
| Quality control | fastp | - |
| Host removal | KneadData | - |
| Assembly | MetaSPAdes | - |
| AMR annotation | RGI (CARD database) | - |
| MAG binning | MetaBAT2 | - |
| MAG quality | CheckM | - |
| Taxonomy | GTDB-Tk r214 | - |
| Community profiling | MetaPhlAn | 4.x |
| Coverage | CoverM | 0.7.0 |

---

---

## Python Analysis Scripts

The `12_analysis/` folder contains a Python script for visualizing and interpreting RGI output across all samples.

| Script | Description |
|--------|-------------|
| `rgi_analysis_plots.py` | Plots 1–11: AMR hit counts, drug class distribution, resistance mechanisms, heatmaps, genus-level analysis, and Shannon diversity index |

**Dependencies:**
```bash
pip install pandas numpy matplotlib seaborn scipy
```

**Usage:**
```bash
cd 12_analysis
python rgi_analysis_plots.py
```

Place your `*_combined_rgi_bins.tsv` files in the same directory before running.

---

## Repository Structure

```
├── 01_download/           # SRA download scripts
├── 02_fastq_conversion/   # FASTQ conversion
├── 03_quality_control/    # fastp trimming
├── 04_host_removal/       # KneadData host removal
├── 05_assembly/           # MetaSPAdes assembly
├── 06_AMR_annotation/     # RGI AMR gene annotation
├── 07_MAG_binning/        # MetaBAT2 binning + bwa mapping
├── 08_MAG_quality/        # CheckM quality assessment
├── 09_taxonomy/           # GTDB-Tk classification
├── 10_community_profiling/ # MetaPhlAn profiling
├── 11_coverage/           # CoverM coverage estimation
├── 12_analysis/
│   └── rgi_analysis_plots.py  # All visualization scripts (Plots 1-10)
└── README.md

```

---

## Step 1 — Data Download

```bash
# Create file with SRR accession numbers
# prefetch downloads .sra files
prefetch --option-file selected_files.txt
```

---

## Step 2 — FASTQ Conversion

```bash
mkdir -p fastq_files
while read srr; do
    fasterq-dump "$srr" -O fastq_files
done < selected_files.txt

# Compress FASTQ files
gzip fastq_files/*.fastq
```

---

## Step 3 — Quality Control (fastp)

```bash
conda activate metagenomics
mkdir -p trimmed qc

for r1 in raw_fastq/*_1.fastq.gz; do
    r2=${r1/_1.fastq.gz/_2.fastq.gz}
    base=$(basename $r1 _1.fastq.gz)
    fastp \
        -i $r1 -I $r2 \
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
done
```

---

## Step 4 — Host Read Removal (KneadData)

```bash
conda activate kneaddata_env
mkdir -p host_removed logs

# Download human genome reference
kneaddata_database --download human_genome bowtie2 kneaddata_db

# Run KneadData for each sample
kneaddata \
    -i1 trimmed/${SAMPLE}_1.trim.fastq.gz \
    -i2 trimmed/${SAMPLE}_2.trim.fastq.gz \
    -db kneaddata_db \
    -o host_removed/${SAMPLE} \
    --output-prefix ${SAMPLE} \
    -t 8 \
    --remove-intermediate-output \
    --log logs/${SAMPLE}_kneaddata.log
```

---

## Step 5 — Metagenomic Assembly (MetaSPAdes)

```bash
conda activate metagenomics

metaspades.py \
    -1 ${SAMPLE}_1.fastq.gz \
    -2 ${SAMPLE}_2.fastq.gz \
    -o ${SAMPLE}_metaspades \
    --threads 16 \
    --meta
```

---

## Step 6 — AMR Gene Annotation (RGI)

```bash
conda activate rgi

rgi main \
    -i contigs.fasta \
    --output_file rgi_output.txt
```

AMR genes annotated against the **CARD (Comprehensive Antibiotic Resistance Database)**.

---

## Step 7 — MAG Binning

```bash
# Map reads to contigs
bwa mem -t 24 contigs.fasta \
    ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz | \
    samtools view -@ 4 -bS - | \
    samtools sort -@ 8 -o ${SAMPLE}_sorted.bam

# Index BAM
samtools index ${SAMPLE}_sorted.bam

# Generate coverage depth
jgi_summarize_bam_contig_depths \
    --outputDepth depth.txt \
    ${SAMPLE}_sorted.bam

# Bin with MetaBAT2
conda activate Basic_protocol_MAG1
metabat2 \
    -i contigs.fasta \
    -a depth.txt \
    -o bins/bin \
    -m 2000 \
    --minContigDepth 2
```

---

## Step 8 — MAG Quality Assessment (CheckM)

```bash
conda activate checkm

checkm lineage_wf sample_bins checkm_out -t 8 -x fa
checkm qa checkm_out/lineage.ms checkm_out \
    -f checkm_quality_summary.tsv --tab_table

# Filter: completeness >= 80%, contamination <= 10%
awk -F'\t' '{print $1"\t"$12"\t"$13"\t"($12-5*$13)}' \
    checkm_quality_summary.tsv | \
    awk -F'\t' '$2 >= 80 && $3 <= 10 && $4 >= 50' \
    > HQ_MAGs_filtered.tsv

# Copy HQ MAGs
mkdir HQ_MAGs
for k in $(awk '{print $1}' HQ_MAGs_filtered.tsv); do
    cp sample_bins/${k}.fa HQ_MAGs/
done
```

---

## Step 9 — Taxonomic Classification (GTDB-Tk)

```bash
# Download GTDB-Tk reference data
wget https://data.gtdb.ecogenomic.org/releases/release214/214.1/auxiliary_files/gtdbtk_r214_data.tar.gz
tar -xvzf gtdbtk_r214_data.tar.gz

# Classify MAGs
gtdbtk classify_wf \
    --genome_dir HQ_MAGs \
    --out_dir gtdbtk_out \
    --cpus 8 \
    --extension fa
```

---

## Step 10 — Community Profiling (MetaPhlAn)

```bash
# Installation
conda create -n metaphlan_env -c conda-forge -c bioconda metaphlan -y
conda activate metaphlan_env

# Run MetaPhlAn
metaphlan \
    ${SAMPLE}_1.fastq.gz,${SAMPLE}_2.fastq.gz \
    --input_type fastq \
    --nproc 16 \
    --mapout ${SAMPLE}_map.txt \
    -o ${SAMPLE}_profile.txt
```

---

## Step 11 — Coverage Estimation (CoverM)

```bash
# Installation
conda create -n coverm_env -c bioconda -c conda-forge coverm -y
conda activate coverm_env

# Convert .fa to .fna
cd HQ_MAGs
for f in *.fa; do mv "$f" "${f%.fa}.fna"; done
cd ..

# Run CoverM
coverm genome \
    --coupled ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz \
    --genome-fasta-directory HQ_MAGs \
    --output-file ${SAMPLE}_coverage.tsv \
    --threads 16
```

---

## Author

**Gayathry Dinesh**
M.Sc. Genomic Science, Central University of Kerala, India
CSIR-UGC NET Qualified for PhD Eligibility in Life Sciences (2025)
gayathrydinesh@gmail.com | [LinkedIn](https://linkedin.com/in/gayathrydinesh)

---

## Citation

If you use this pipeline, please cite:
> Gayathry Dinesh.
Exploring Antimicrobial Resistance Gene Reservoirs in the Gut Microbiome of Type 1 Diabetes Patients through Metagenomics.
M.Sc. Dissertation, Central University of Kerala, 2026.

