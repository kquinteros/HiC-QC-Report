#!/bin/bash

# Batch Hi-C QC Pipeline
# Process multiple samples from a config file
# Usage: ./batch_hic_pipeline.sh samples.txt

set -e
set -u
set -o pipefail

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <samples_file>"
    echo ""
    echo "Sample file format (tab-separated):"
    echo "sample_name  reference_genome.fasta  R1.fastq.gz  R2.fastq.gz"
    echo ""
    echo "Example samples.txt:"
    echo "sample1  genome1.fasta  data/s1_R1.fq.gz  data/s1_R2.fq.gz"
    echo "sample2  genome1.fasta  data/s2_R1.fq.gz  data/s2_R2.fq.gz"
    exit 1
fi

SAMPLES_FILE=$1

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file not found: $SAMPLES_FILE"
    exit 1
fi

# Count total samples
TOTAL=$(grep -v '^#' "$SAMPLES_FILE" | grep -v '^$' | wc -l)
CURRENT=0

echo "========================================"
echo "Batch Hi-C QC Pipeline"
echo "Total samples: $TOTAL"
echo "========================================"

# Process each sample
while IFS=$'\t' read -r SAMPLE REFERENCE R1 R2 || [ -n "$SAMPLE" ]; do
    # Skip comments and empty lines
    [[ "$SAMPLE" =~ ^#.*$ ]] && continue
    [[ -z "$SAMPLE" ]] && continue
    
    CURRENT=$((CURRENT + 1))
    
    echo ""
    echo "========================================"
    echo "Processing sample $CURRENT of $TOTAL: $SAMPLE"
    echo "========================================"
    
    # Run pipeline for this sample
    ./hic_pipeline.sh "$SAMPLE" "$REFERENCE" "$R1" "$R2"
    
done < "$SAMPLES_FILE"

echo ""
echo "========================================"
echo "Batch processing complete!"
echo "Processed $TOTAL samples"
echo "========================================"
