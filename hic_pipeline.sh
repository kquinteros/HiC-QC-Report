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
HIC_QC_SCRIPT="$SCRIPT_DIR/hic_qc/hic_qc.py"

#You can pass your own Threshold values by editinng this JSON file
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
FASTP_LOG="logs/fastp_${SAMPLE}.log"
FASTP_HTML="qc/fastp/${SAMPLE}.html"
FASTP_JSON="qc/fastp/${SAMPLE}.json"

echo "[$(date)] Step 2: Trimming adapters with fastp..."
if [[ -s "$TRIMMED_R1" && -s "$TRIMMED_R2" && -s "$FASTP_HTML" && -s "$FASTP_JSON" && -s "$FASTP_LOG" ]]; then
    echo "[$(date)] fastp outputs + log exist, skipping..."
else
    fastp \
        -i "$R1" \
        -I "$R2" \
        -o "$TRIMMED_R1" \
        -O "$TRIMMED_R2" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --html "$FASTP_HTML" \
        --json "$FASTP_JSON" \
        2>&1 | tee "$FASTP_LOG"
    echo "[$(date)] Adapter trimming complete"
fi

# Step 3: Align with BWA-MEM, mark duplicates with samblaster, convert to BAM
BWA_LOG="logs/bwa_mem_${SAMPLE}.log"
SBL_LOG="logs/samblaster_${SAMPLE}.log"

echo "[$(date)] Step 3: Aligning reads with BWA-MEM + samblaster..."
if [[ -s "$BAM" && -s "$BWA_LOG" && -s "$SBL_LOG" ]]; then
    echo "[$(date)] BAM + logs exist, skipping alignment..."
else
    bwa mem -t "$THREADS" -M "$REFERENCE" "$TRIMMED_R1" "$TRIMMED_R2" | \
        samblaster 2>> "$SBL_LOG" | \
        samtools view -Sb - > "$BAM" 2>> "$BWA_LOG"
    echo "[$(date)] Alignment complete"
fi

# Step 4: Run Hi-C QC
echo "[$(date)] Step 4: Running Hi-C QC..." 

ln -s "$(pwd)/hic_qc/collateral" "results/${REF_NAME}/${SAMPLE}/collateral"

python "$HIC_QC_SCRIPT" \
    -b "$BAM" \
    -n "$HIC_QC_READS" \
    -o "${OUTDIR}/${SAMPLE}" \
    --thresholds "$THRESHOLDS_JSON" \
    --sample_type genome \
    2>&1 | tee "logs/hic_qc_${SAMPLE}.log"
echo "[$(date)] Hi-C QC complete"


# Step 5: De-absolute report paths + convert to PDF with pandoc
echo "[$(date)] Step 5: Cleaning report paths and converting to PDF with pandoc..."

# Input report produced by HiC-QC (adjust if your filename differs)
REPORT_IN="${OUTDIR}/${SAMPLE}_qc_report.html"
REPORT_CLEAN_MD="${OUTDIR}/${SAMPLE}_qc_report.cleaned.md"
REPORT_PDF="${OUTDIR}/${SAMPLE}_qc_report.pandoc.pdf"

[[ -f "$REPORT_IN" ]] || {
  echo "ERROR: Report not found: $REPORT_IN" >&2
  echo "Adjust REPORT_IN in the script to match your HiC-QC output name." >&2
  exit 1
}

# Prefixes to remove
PREFIX1="${PWD}/hic_qc/"
PREFIX2="${PWD}/results/${REF_NAME}/${SAMPLE}/"

# Remove the prefixes everywhere in the report text
sed -e "s|${PREFIX1}||g" \
    -e "s|${PREFIX2}||g" \
    "$REPORT_IN" > "$REPORT_CLEAN_MD"

# Convert cleaned markdown-ish file to PDF
command -v pandoc >/dev/null 2>&1 || { echo "ERROR: pandoc not found in PATH" >&2; exit 1; }

pandoc "$REPORT_CLEAN_MD" \
  -o "$REPORT_PDF" \
  --resource-path=".:${OUTDIR}" \
  --css "${OUTDIR}/collateral/style.css" \
  --pdf-engine weasyprint

echo "[$(date)] Report cleaned: $REPORT_CLEAN_MD"
echo "[$(date)] Pandoc PDF written: $REPORT_PDF"

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
