#!/bin/bash
# Author: Kevin Quinteros
# Hi-C QC Pipeline Bash Script
# Usage: ./hic_pipeline.sh <sample_name> <reference_genome> <r1_fastq> <r2_fastq> 

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ----------------------------
# Resolve script directory safely
# ----------------------------
if [[ -n "${BASH_SOURCE[0]:-}" && -f "${BASH_SOURCE[0]}" ]]; then
    SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
else
    echo "ERROR: Cannot determine script directory. Run the script as a file, not via stdin." >&2
    echo "Use: ./hic_pipeline.sh or bash hic_pipeline.sh" >&2
    exit 1
fi

# ----------------------------
# Required files
# ----------------------------
HIC_QC_SCRIPT="$SCRIPT_DIR/hic_qc.hic_qc"
THRESHOLDS_JSON="$SCRIPT_DIR/hic_qc/collateral/thresholds.json"

# ----------------------------
# Check arguments
# ----------------------------
# Check arguments (now 4, 5, OR 6 are acceptable)
if [[ "$#" -lt 4 || "$#" -gt 6 ]]; then
    echo "Usage: $0 <sample> <reference.fa> <r1.fq.gz> <r2.fq.gz> [threads] [hic_qc_reads]" >&2
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
HIC_QC_READS=${6:-1000000}   # Default: 1000000

# Derived variables
REF_NAME=$(basename "$REFERENCE" | sed 's/\.[^.]*$//')
OUTDIR="results/${REF_NAME}/${SAMPLE}"
TRIMMED_R1="trimmed/${SAMPLE}.r1.trimmed.fq.gz"
TRIMMED_R2="trimmed/${SAMPLE}.r2.trimmed.fq.gz"
BAM="${OUTDIR}/${SAMPLE}.bam"


#check that required file exist
for f in "$REFERENCE" "$R1" "$R2" "$HIC_QC_SCRIPT" "$THRESHOLDS_JSON"; do
    [[ -f "$f" ]] || { echo "ERROR: File not found: $f" >&2; exit 1; }
done

# Create output directories
mkdir -p "$OUTDIR" trimmed qc/fastp logs

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
fastp \
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
bwa mem -t "$THREADS" -M "$REFERENCE" "$TRIMMED_R1" "$TRIMMED_R2" | \
    samblaster 2>> "logs/samblaster_${SAMPLE}.log" | \
    samtools view -Sb - > "$BAM" 2>> "logs/bwa_mem_${SAMPLE}.log"
echo "[$(date)] Alignment complete"

# Step 4: Run Hi-C QC
echo "[$(date)] Step 4: Running Hi-C QC..." 
python -m "$HIC_QC_SCRIPT" \
    -b "$BAM" \
    -n "$HIC_QC_READS" \
    -o "${OUTDIR}/${SAMPLE}" \
    --thresholds "$THRESHOLDS_JSON" \
    --sample_type genome \
    -r \
    2>&1 | tee "logs/hic_qc_${SAMPLE}.log"
echo "[$(date)] Hi-C QC complete"


# ----------------------------
# Done
# ----------------------------

echo "========================================"
echo "[$(date)] Pipeline complete!"
echo "Results in: $OUTDIR"
echo "  - ${SAMPLE}.bam"
echo "  - ${SAMPLE}_hic_qc_report.pdf"
echo "  - ${SAMPLE}_Read_mate_dist.pdf"
echo "Fastp QC: qc/fastp/${SAMPLE}.html"
echo "========================================"
