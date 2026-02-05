# Hi-C QC Pipeline - Bash Scripts

Simple bash scripts to run Hi-C QC pipeline.

## Pipeline Steps

1. **Trim adapters** with fastp
2. **Index reference** genome with BWA (if needed)
3. **Align reads** with BWA-MEM → samblaster → samtools
4. **Run Hi-C QC** analysis

## Dependencies 

This script is dependent on the following tools: bwa, samtools, samblaster, and fastp. 
The `hic_qc.py` script by Phase genomics has additional dependencies, so please check their [git repo page](https://github.com/phasegenomics/hic_qc/tree/master) for more information.

This pipeline requires the following external tools: bwa, samtools, samblaster, and fastp.
The hic_qc.py script from Phase Genomics has additional dependencies; please refer to the [Phase Genomics GitHub repository](https://github.com/phasegenomics/hic_qc/tree/master) for full details.

A Conda YAML file is provided to create an environment with the required packages.


```bash
git clone https://github.com/phasegenomics/hic_qc.git
conda env create -f env-hic-tools.yaml
conda activate hic_tools
conda install conda-forge::wkhtmltopdf
pip install --no-deps -e hic_qc
chmod +x batch_hic_pipeline.sh hic_pipeline.sh
```


## Quick Start

### Single Sample

```bash
bash hic_pipeline.sh sample1 genome.fasta sample1_R1.fq.gz sample1_R2.fq.gz
```

### Multiple Samples

1. Edit `samples.txt` with your samples (tab-separated):
```
sample1  genome.fasta  data/s1_R1.fq.gz  data/s1_R2.fq.gz
sample2  genome.fasta  data/s2_R1.fq.gz  data/s2_R2.fq.gz
```

2. Run batch script:
```bash
 bash batch_hic_pipeline.sh samples.txt
```

## Output Structure

```
results/
└── genome_name/
    └── sample_name/
        ├── sample_name.bam
        ├── sample_name_hic_qc_report.pdf
        └── sample_name_Read_mate_dist.pdf

qc/
└── fastp/
    ├── sample_name.html
    └── sample_name.json

logs/
├── bwa_index_*.log
├── fastp_*.log
├── bwa_mem_*.log
├── samblaster_*.log
└── hic_qc_*.log
```

## Configuration

You can specify the number of threads and reads to sample for QC by passing optional parameters to hic_pipeline.sh.

- `THREADS=8` - Number of CPU threads
- `HIC_QC_READS=1000000` - Number of reads to sample for QC

## Notes

- Reference genome is indexed once and reused
- All logs saved to `logs/` directory
- Fastp QC reports in `qc/fastp/`
- Results organized by reference genome name

## Parallel Processing

To run multiple samples in parallel:
```bash
# Process 4 samples at once
cat samples.txt | parallel -j 4 --colsep '\t' ./hic_pipeline.sh {1} {2} {3} {4} 
```

## Troubleshooting

**Command not found errors:**
- Install missing tools with conda

**Permission denied:**
- Run: `chmod +x *.sh`

**Out of memory:**
- Reduce `THREADS` in the script

**BWA index fails:**
- Ensure reference FASTA file exists and is readable
