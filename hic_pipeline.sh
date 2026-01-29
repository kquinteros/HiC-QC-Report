#!/bin/bash

# Hi-C QC Pipeline Bash Script
# Usage: ./hic_pipeline.sh <sample_name> <reference_genome> <r1_fastq> <r2_fastq> <

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# Check arguments (now 4, 5, OR 6 are acceptable)
if [ "$#" -lt 4 ] || [ "$#" -gt 6 ]; then
    echo "Usage: $0 <sample_name> <reference_genome> <r1_fastq> <r2_fastq> [threads] [hic_qc_reads]"
    echo "Example: $0 sample1 genome.fasta sample1_R1.fq.gz sample1_R2.fq.gz"
    echo "Example: $0 sample1 genome.fasta sample1_R1.fq.gz sample1_R2.fq.gz 16"
    echo "Example: $0 sample1 genome.fasta sample1_R1.fq.gz sample1_R2.fq.gz 16 2000000"
    echo ""
    echo "Defaults:"
    echo "  threads: 8"
    echo "  hic_qc_reads: 1000000"
    exit 1
fi

# Arguments
SAMPLE=$1
REFERENCE=$2
R1=$3
R2=$4
THREADS=${5:-8}              # Default: 8 threads
HIC_QC_READS=${6:-1000000}   # Default: 1
JSON="hic_qc/collateral/thresholds.json"

# Derived variables
REF_NAME=$(basename "$REFERENCE" | sed 's/\.[^.]*$//')
OUTDIR="results/${REF_NAME}/${SAMPLE}"
TRIMMED_R1="trimmed/${SAMPLE}.r1.trimmed.fq.gz"
TRIMMED_R2="trimmed/${SAMPLE}.r2.trimmed.fq.gz"
BAM="${OUTDIR}/${SAMPLE}.bam"

# Create output directories
mkdir -p "$OUTDIR"
mkdir -p trimmed
mkdir -p qc/fastp
mkdir -p logs

echo "========================================"
echo "Hi-C QC Pipeline"
echo "========================================"
echo "Sample: $SAMPLE"
echo "Reference: $REFERENCE"
echo "Output directory: $OUTDIR"
echo "========================================"

# Step 1: Index reference genome if not already indexed
echo "[$(date)] Step 1: Indexing reference genome..."
if [ ! -f "${REFERENCE}.bwt" ]; then
    bwa index "$REFERENCE" 2>&1 | tee "logs/bwa_index_${REF_NAME}.log"
    echo "[$(date)] Indexing complete"
else
    echo "[$(date)] Index already exists, skipping..."
fi

# Step 2: Trim adapters with fastp
echo "[$(date)] Step 2: Trimming adapters with fastp..."
conda run -n hic_tools fastp \
    -i "$R1" \
    -I "$R2" \
    -o "$TRIMMED_R1" \
    -O "$TRIMMED_R2" \
    --thread "$THREADS" \
    --detect_adapter_for_pe \
    --html "qc/fastp/${SAMPLE}.html" \
    --json "qc/fastp/${SAMPLE}.json" \
    2>&1 | tee "logs/fastp_${SAMPLE}.log"
echo "[$(date)] Adapter trimming complete"

# Step 3: Align with BWA-MEM, mark duplicates with samblaster, convert to BAM
echo "[$(date)] Step 3: Aligning reads with BWA-MEM + samblaster..."
conda run -n hic_tools \
bwa mem -t "$THREADS" -M "$REFERENCE" "$TRIMMED_R1" "$TRIMMED_R2" | \
    samblaster 2>> "logs/samblaster_${SAMPLE}.log" | \
    samtools view -Sb - > "$BAM" 2>> "logs/bwa_mem_${SAMPLE}.log"
echo "[$(date)] Alignment complete"

# Step 4: Run Hi-C QC
echo "[$(date)] Step 4: Running Hi-C QC..."
conda run -n hic_qc python hic_qc/hic_qc.py \
    -b "$BAM" \
    -n "$HIC_QC_READS" \
    -o "${OUTDIR}/${SAMPLE}" \
    --thresholds "${JSON}" \
    --sample_type genome \
    -r \
    2>&1 | tee "logs/hic_qc_${SAMPLE}.log"
echo "[$(date)] Hi-C QC complete"

echo "========================================"
echo "[$(date)] Pipeline complete!"
echo "Results in: $OUTDIR"
echo "  - ${SAMPLE}.bam"
echo "  - ${SAMPLE}_hic_qc_report.pdf"
echo "  - ${SAMPLE}_Read_mate_dist.pdf"
echo "Fastp QC: qc/fastp/${SAMPLE}.html"
echo "========================================"
