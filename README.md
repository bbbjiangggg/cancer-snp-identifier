# Pan-Cancer SNP Discovery and Gene-Level Prioritization Pipeline

## Overview

This repository contains the computational pipeline for **"Pan-Cancer SNP Discovery and Gene-Level Prioritization Using Machine Learning and Functional Annotation"** by Bowen Jiang and Roberto Aguilar.

The pipeline identifies and functionally prioritizes shared SNPs across whole-genome sequencing data from seven cancer types: breast, lung, prostate, gastric, liver, ovarian, and pancreatic cancers.

## Repository Structure

```
cancer-snp-identifier/
├── snp_analysis/                    # Core SNP discovery pipeline
│   └── snp_analysis_pipeline_v3/    # Latest version (v3)
│       ├── main_pipeline.py         # Main pipeline orchestrator
│       ├── data_processing.py       # Data download and preprocessing
│       ├── file_handling.py         # File I/O operations
│       ├── command_execution.py     # Shell command execution
│       ├── logging_module.py        # Logging functionality
│       ├── path_management.py       # Directory management
│       └── utils.py                 # Utility functions
│
├── final_analysis_comb/             # VCF combination and processing
│   ├── final_combo.py               # Combine VCF files across samples
│   ├── vcf_2_comb_v2.3.1.py        # VCF processing utilities
│   ├── gene_id_v1.1.py              # Gene ID mapping
│   ├── snp_effect.py                # SNP functional effect annotation
│   ├── normal.py                    # Normal control data download
│   └── genomic/                     # Genomic-specific processing
│       ├── genomic_split_vcf_v1.0.py
│       ├── genomic_vcf_2_comb_v1.1.py
│       └── vcf_2_comb_v2.4.py
│
├── indexing/                        # Reference genome indexing
│   ├── all_chrom_bowtie_v1.0.py     # Bowtie2 index generation
│   ├── all_chrom_bwa_v1.0.py        # BWA index generation
│   ├── bowtie_ind_v1.1.py           # Individual chromosome indexing
│   └── bwa_ind_ref_chrom_v1.1.py    # BWA reference indexing
│
├── installation/                    # Software installation scripts
│   ├── install_linux_prog_v1.2.py   # Linux installation (latest)
│   ├── install_mac_prog_v1.0.py     # macOS installation
│   ├── install_mac_prog_v1.0_brew.py # macOS with Homebrew
│   └── install_mac_prog_v1.0_port.py # macOS with MacPorts
│
├── snpnexus/                        # SNP functional annotation
│   └── rawdata_2_snpnexus_v1.1.py   # SNPnexus integration
│
└── tools/                           # Utility scripts
    ├── effect_categories.py         # SNP effect categorization
    ├── snp_effect_count.py          # Effect counting utilities
    ├── comb_2_csv_5of6.py           # Data format conversion
    ├── comp_male_fem_snps.py        # Gender-specific SNP analysis
    ├── final_vcf_transfer_v1.0.py   # VCF file management
    ├── remove_vcf_zero.py           # Clean empty VCF files
    └── position_count.sh            # Position counting utilities
```

## Pipeline Workflow

### 1. Data Acquisition and Preprocessing
- Downloads whole-genome sequencing data from SRA
- Processes FASTQ files with quality control (fastp)
- Supports both single-end and paired-end reads

### 2. Reference Genome Alignment
- Aligns reads to GRCh38 reference genome
- Uses Bowtie2 for fast local alignment
- Processes chromosomes 1-22, X, and Y individually

### 3. Variant Calling and Filtering
- Generates BAM files using SAMtools
- Calls variants using BCFtools
- Applies quality filters (QUAL ≥ 30, DP ≥ 10)
- Normalizes multiallelic variants to biallelic format

### 4. VCF Harmonization
- Combines sample-level VCFs by chromosome and cancer type
- Compresses and indexes VCF files using bgzip and tabix
- Creates multisample VCF files for downstream analysis

### 5. Functional Annotation
- Annotates SNPs using Ensembl VEP API
- Classifies molecular consequences (missense, regulatory, etc.)
- Maps SNPs to gene identifiers using GTF reference files

## Requirements

### System Dependencies
- **Bowtie2** (≥2.5.1) - Fast read alignment
- **BWA** (≥0.7.17) - Alternative aligner
- **SAMtools** (≥1.15) - BAM file processing
- **BCFtools** (≥1.15) - Variant calling and processing
- **tabix** - VCF indexing
- **fastp** (≥0.23.2) - Read preprocessing

### Python Dependencies
- **Python** (≥3.10)
- **pandas** - Data manipulation
- **requests** - API interactions
- **termcolor** - Colored terminal output
- **tqdm** - Progress bars
- **pyfiglet** - ASCII art banners

## Usage

### 1. Environment Setup
```bash
# Install system dependencies (choose your platform)
python installation/install_linux_prog_v1.2.py    # Linux
python installation/install_mac_prog_v1.0.py      # macOS

# Create reference genome indexes
python indexing/all_chrom_bowtie_v1.0.py          # Bowtie2 indexes
python indexing/all_chrom_bwa_v1.0.py             # BWA indexes
```

### 2. Prepare Input Data
Create a text file with SRA accession numbers (one per line):
```
SRR12345678
SRR12345679
SRR12345680
```

### 3. Run SNP Discovery Pipeline
```bash
cd snp_analysis/snp_analysis_pipeline_v3/
python main_pipeline.py
```

The pipeline will prompt for:
- Accession list file selection
- Chromosomes to analyze (e.g., "1,2,3" or "all")
- Number of samples to process

### 4. Combine and Process VCFs
```bash
cd final_analysis_comb/
python final_combo.py                    # Combine VCFs across samples
python gene_id_v1.1.py                  # Map to gene identifiers
python snp_effect.py                    # Annotate functional effects
```

### 5. Functional Annotation
```bash
cd snpnexus/
python rawdata_2_snpnexus_v1.1.py       # Additional annotation via SNPnexus
```

## Input/Output Formats

### Input
- **SRA accession numbers** (text file)
- **Reference genome** (GRCh38 FASTA files)
- **Gene annotation** (GTF files from Ensembl)

### Output
- **VCF files** - Variant calls per sample/chromosome
- **Combined VCFs** - Merged variants across samples
- **Annotated CSV** - SNPs with functional effects and gene mappings
- **Log files** - Processing logs and QC metrics

## Data Sources

The pipeline uses data from the following BioProjects:

| Cancer Type | BioProject | Description |
|-------------|------------|-------------|
| Breast | PRJNA734808 | Breast cancer WGS |
| Gastric | PRJNA307236 | Gastric cancer WGS |
| Lung | PRJNA448888 | Lung cancer WGS |
| Liver | PRJNA504942 | Liver cancer WGS |
| Ovarian | PRJEB28664, PRJEB47696 | Ovarian cancer WGS |
| Pancreatic | PRJNA646156 | Pancreatic cancer WGS |
| Prostate | PRJNA513939, PRJNA554329 | Prostate cancer WGS (including cfDNA) |
| Normal | PRJEB13738 | Normal control samples |

## ⚠️ Missing Components

**IMPORTANT**: The following critical components from the manuscript are currently missing from this repository:

### Machine Learning Pipeline
- **XGBoost + SVM ensemble classifier** implementation
- **Binary encoding** of SNP matrices for ML
- **SMOTE** for class imbalance handling
- **Model training and validation** scripts
- **Performance evaluation** (accuracy, ROC AUC, F1 score)

### Statistical Analysis
- **Binomial testing** for SNP significance
- **Manhattan plot** generation
- **UpSet plot** visualization
- **Circos plot** generation for genome-wide visualization

### Gene Classification
- **COSMIC/OncoKB** gene annotation
- **ImmPort** immune gene classification
- **Novel candidate** identification

Without these components, the current pipeline can only process raw sequencing data to annotated SNPs but cannot reproduce the main findings of the manuscript.

## Publication

If you use this pipeline, please cite:

```
Jiang, B., & Aguilar, R. Pan-Cancer SNP Discovery and Gene-Level Prioritization Using Machine Learning and Functional Annotation. [Journal], [Year].
```

## Contact

- **Bowen Jiang** - b.jiang@duke.edu
- **Roberto Aguilar** - Western Reserve Academy

## License

This project is part of academic research. Please contact the authors for usage permissions.

---

**Note**: This repository contains the data processing pipeline. The machine learning analysis code will be added in a future update to fully reproduce the manuscript results.