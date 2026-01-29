#!/bin/bash

# Set up script for Hi-C QC pipeline
# This scripts downloads and installs necessary files and 
# assumes that user has git and codna installed. 

set -e
set -u
set -o pipefail

# Initialize conda for this script
eval "$(conda shell.bash hook)"

# git clone hic.py repo from phase genomics 
git clone https://github.com/phasegenomics/hic_qc.git
cd hic_qc/
conda env create -n hic_qc --file env.minimal.yml
conda activate hic_qc
conda install conda-forge::wkhtmltopdf
pip install --no-deps -e .
conda deactivate

# Install conda environment
conda env create --file env-hic-tools.yaml

# Create directories
mkdir -p data
mkdir -p references

# Make scripts executable
chmod +x batch_hic_pipeline.sh hic_pipeline.sh
