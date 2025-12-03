# Metabarcoding Processing Pipelines

This repository contains scripts for quality checking, preprocessing, clustering, and taxonomic annotation of amplicon sequencing data. The pipelines include: `1.1-quality_check_fastp.sh`, `1.2-check_primers.py`, `1.3-primer_removal_cutadapt.sh`, `2-dada2_pipeline.R`, and `3-taxa_annot.R`, written in Bash, Python, and R.

## Repository structure

```
.
├── LICENSE                                # License file
├── README.md                              # This file
├── environment.yml                        # Conda/Mamba environment specification
└── modules/                               # Main analysis scripts
    ├── 1.1-quality_check_fastp.sh         # Quality check with fastp
    ├── 1.2-check_primers.py               # Primer checking utility
    ├── 1.3-primer_removal_cutadapt.sh     # Primer removal with cutadapt
    ├── 2-dada2_pipeline.R                 # DADA2 ASV generation pipeline
    ├── 3-taxa_annot.R                     # Taxonomic annotation pipeline
    ├── toolbox.R                          # Utility functions
    └── conf.sh                            # Configuration file
```

## Installation instructions

Clone the repository:

```bash
git clone https://github.com/ocean-biometrics/metabarcoding-processing-obm.git
cd metabarcoding-processing-obm
```

### Installation with mamba

All dependencies can be installed using mamba (or conda).

First, check if mamba is installed:

```bash
command -v mamba
```

If mamba is not installed, you can install it via:

- **Miniforge (recommended)**: https://github.com/conda-forge/miniforge
- **Mambaforge**: https://github.com/conda-forge/miniforge#mambaforge
- Or install mamba into an existing conda installation:
  ```bash
  conda install -n base -c conda-forge mamba
  ```

Once mamba is installed, create a new environment using the provided `environment.yml` file:

```bash
mamba env create -f environment.yml
```

Then activate the environment:

```bash
conda activate metabarcoding-processing-obm
```

## How to use

### 1.1-quality_check_fastp.sh

**1.1-quality_check_fastp.sh**: This Bash script generates quality control reports for paired-end FASTQ files using fastp. The main tasks include:

- Analyzing read quality metrics for R1 and R2 reads
- Generating comprehensive HTML and JSON reports
- Computing statistics on base quality, read counts, and filtering metrics
- Running in report-only mode (no reads are modified)

The output consists of:

- `<sample_name>_fastp.html`: Interactive HTML report with quality plots and statistics
- `<sample_name>_fastp.json`: JSON format report with detailed metrics
- `<sample_name>_fastp.log`: Processing log file
- `summary.tsv`: Summary statistics table for all samples
- `summary_report.txt`: Human-readable summary report

To see the help run `bash modules/1.1-quality_check_fastp.sh --help`

```
Usage: 1.1-quality_check_fastp.sh <options>

Options:
    --help
        Print this help message and exit

    --input_dir=CHAR
        Directory containing input FASTQ files (required)

    --output_dir=CHAR
        Directory to output generated data (QC reports and plots) (required)

    --pattern_r1=CHAR
        Pattern for R1 FASTQ files [default=_R1_001.fastq.gz]

    --pattern_r2=CHAR
        Pattern for R2 FASTQ files [default=_R2_001.fastq.gz]

    --nslots=NUM
        Number of threads to use [default=12]

    --html_report=t|f
        Generate HTML report [default=t]

    --json_report=t|f
        Generate JSON report [default=t]

    --overwrite=t|f
        Overwrite previous directory [default=f]
```

### 1.2-check_primers.py

**1.2-check_primers.py**: This Python script checks for primer presence in FASTQ files with support for IUPAC ambiguity codes. The main tasks include:

- Subsampling reads from paired-end FASTQ files
- Detecting primers in all orientations (forward, complement, reverse, reverse-complement)
- Counting primer hits with IUPAC ambiguity code support
- Generating statistics on primer detection rates

The output consists of:

- `<sample_name>_primer_check.csv`: CSV file with primer detection percentages (or counts) for each orientation per sample. Columns include primer orientation combinations (e.g., FwdReads.FwdPrimer.Forward, RevReads.RevPrimer.ReverseComplement)
- `primer_sequences.csv`: CSV file listing the primer sequences in all four orientations (Forward, Complement, Reverse, ReverseComplement) for both forward and reverse primers
- Helps determine if primer removal is needed before downstream processing

To see the help run `python modules/1.2-check_primers.py --help`

```
usage: check_primers.py [-h] -i INPUT_DIR -o OUTPUT_DIR --suffix_r1 SUFFIX_R1 
                        --suffix_r2 SUFFIX_R2 --primer_fwd PRIMER_FWD 
                        --primer_rev PRIMER_REV [-s SUBSAMPLE_SIZE] [--counts]

Subsample reads and count IUPAC‐aware primer hits

options:
  -h, --help            show this help message and exit
  -i INPUT_DIR, --input_dir INPUT_DIR
                        Path to input directory
  -o OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Path to output directory
  --suffix_r1 SUFFIX_R1
                        Glob pattern for forward reads (R1)
  --suffix_r2 SUFFIX_R2
                        Glob pattern for reverse reads (R2)
  --primer_fwd PRIMER_FWD
                        Forward primer sequence
  --primer_rev PRIMER_REV
                        Reverse primer sequence
  -s SUBSAMPLE_SIZE, --subsample_size SUBSAMPLE_SIZE
                        Reads to subsample per file (default: 1000)
  --counts              Write raw counts instead of percentages (default: percentages)
```

### 1.3-primer_removal_cutadapt.sh

**1.3-primer_removal_cutadapt.sh**: This Bash script removes primer sequences from paired-end FASTQ files using cutadapt. The main tasks include:

- Trimming forward primers from R1 reads
- Trimming reverse primers (reverse complement) from R2 reads
- Filtering reads by minimum length after trimming
- Support for IUPAC ambiguity codes in primers
- Optional discarding of reads without detected primers
- Computing trimming statistics

The output consists of:

- `<sample_name>_R1_trimmed.fastq.gz`: Trimmed R1 reads (optionally compressed)
- `<sample_name>_R2_trimmed.fastq.gz`: Trimmed R2 reads (optionally compressed)
- `<sample_name>_cutadapt.log`: Detailed cutadapt log per sample
- `summary.tsv`: Summary statistics table for all samples
- `summary_report.txt`: Human-readable summary report with primer sequences and parameters

To see the help run `bash modules/1.3-primer_removal_cutadapt.sh --help`

```
Usage: 1.3-primer_removal_cutadapt.sh <options>

Options:
    --help
        Print this help message and exit

    --input_dir=CHAR
        Directory containing input FASTQ files (required)

    --output_dir=CHAR
        Directory to output trimmed FASTQ files (required)

    --pattern_r1=CHAR
        Pattern for R1 FASTQ files [default=_R1_001.fastq.gz]

    --pattern_r2=CHAR
        Pattern for R2 FASTQ files [default=_R2_001.fastq.gz]

    --primer_fwd=CHAR
        Forward primer sequence (5' to 3') (required)

    --primer_rev=CHAR
        Reverse primer sequence (5' to 3') (required)

    --nslots=NUM
        Number of threads to use [default=12]

    --error_rate=NUM
        Maximum allowed error rate (mismatches/length) [default=0.1]

    --min_overlap=NUM
        Minimum overlap between primer and read [default=3]

    --min_length=NUM
        Discard reads shorter than this after trimming [default=50]

    --discard_untrimmed=t|f
        Discard reads where no primer was found [default=f]

    --compress=t|f
        Compress output files with gzip [default=t]

    --overwrite=t|f
        Overwrite previous directory [default=f]
```

### 2-dada2_pipeline.R

**2-dada2_pipeline.R**: This R pipeline processes raw Illumina paired-end reads from amplicon sequencing samples and generates Amplicon Sequence Variants (ASVs) using the DADA2 algorithm. The main tasks include:

- Quality filtering and trimming of paired-end reads
- Learning error rates from the data
- Denoising and sample inference to identify ASVs
- Merging paired-end reads
- Chimera detection and removal
- Generating quality and error plots
- Computing sequence statistics (counts and lengths)

The output consists of:

- `asv_table.csv`: CSV file with ASV abundance table (sequences as rows, samples as columns)
- `plots/quality_plot_R1.pdf` and `plots/quality_plot_R2.pdf`: Quality profiles for forward and reverse reads
- `plots/error_plot_R1.pdf` and `plots/error_plot_R2.pdf`: Error rate learning plots
- `plots/nseq_barplot.pdf`: Barplot showing sequence counts through each processing step
- `plots/seq_length_hist.pdf`: Boxplots of merged read lengths per sample
- `tables/nseq_counts.csv`: Table with read counts at each processing step
- `tables/seq_length_stats.csv`: Table with length statistics (mean, sd, min, max) per sample
- `tables/session_info.txt`: Complete session information including input files, parameters used, and package versions

To see the help run `Rscript modules/2-dada2_pipeline.R --help`

```
Usage: ./dada2_pipeline.R <options>
--help                          print this help
--input_dir CHAR                directory with the input raw fastq files (required)
--output_dir CHAR               directory to output generated data (required)
--nslots NUM                    number of threads used (default: 12)
--trunc_r1 NUM                  number of nuc to remove in R1 from the 3' end (default: 250)
--trunc_r2 NUM                  number of nuc to remove in R2 from the 3' end (default: 200)
--pattern_r1 CHAR               pattern of R1 reads to load fastq files (default: _L001_R1_001.fastq.gz)
--pattern_r2 CHAR               pattern of R2 reads to load fastq files (default: _L001_R2_001.fastq.gz)
--bimeras_method CHAR           method to check bimeras: pooled, consensus, per-sample (default: consensus)
--min_overlap NUM               minimum number of nucleotides to overlap in merging (default: 12)
--pooled                        use pooled option when running dada2 (default: TRUE)
--no_pooled                     disable pooled option
--qual_plot                     create quality plots (default: TRUE)
--no_qual_plot                  disable quality plots
--err_plot                      create error plots (default: TRUE)
--no_err_plot                   disable error plots
--save_workspace                save R workspace image (default: TRUE)
--no_save_workspace             disable saving workspace
--overwrite                     overwrite previous directory (default: FALSE)
```

### 3-taxa_annot.R

**3-taxa_annot.R**: This R script assigns taxonomy to ASVs generated by the DADA2 pipeline. It supports three annotation methods:

- **NBC (Naive Bayes Classifier)**: Uses a trained classifier on a reference database
- **NBCandEM (NBC + Exact Matching)**: Combines NBC with exact matching for species-level assignment
- **BLAST**: Uses BLAST search against a nucleotide database with taxonomy mapping

The main tasks include:

- Reading ASV table and sequences
- Running taxonomic classification using the selected method
- Mapping taxonomic assignments to ASVs
- Generating annotated ASV table with taxonomy information

The output consists of:

- Annotated ASV table (CSV format) with taxonomic assignments
  - For NBC/NBCandEM: includes Kingdom, Phylum, Class, Order, Family, Genus, (Species) with bootstrap values
  - For BLAST: includes best hit information (accession, organism name, taxonomy path, percent identity)
- `session_info_taxa_annot.txt`: Session information including input files, parameters used, and package versions

To see the help run `Rscript modules/3-taxa_annot.R --help`

```
Usage: ./taxa_annot.R <options>
--help                          print this help
--input_asv_table CHAR          asv table generated with DADA2 (required)
--input_fasta CHAR              fasta file with sequences to be annotated (optional, will be generated from ASV table if not provided)
--output_asv_table CHAR         output asv table with taxonomic annotation (required)
--method CHAR                   method used to annotate: NBC, NBCandEM, BLAST (default: NBC)
                                NBC: Naive Bayes Classifier; EM: Exact Matching
--evalue NUM                    evalue used in BLAST search (default: 1e-10)
--min_identity NUM              minimum identity used in BLAST search (default: 97)
--train_db CHAR                 training database to run NBC (default: silva_nr_v138_train_seq.fa.gz)
--ref_db CHAR                   reference database to run EM (default: silva_species_assignment_v138.fa.gz)
--blast_db CHAR                 blast formatted database to run BLAST (default: SILVA_138_SSURef_NR99_tax_silva.fasta)
--taxa_map CHAR                 tsv file mapping silva acc with taxonomy (used when running BLAST)
--nslots NUM                    number of threads used (default: 12)
--save_workspace                save R workspace image (default: TRUE)
--no_save_workspace             disable saving workspace
--overwrite                     overwrite previous output (default: FALSE)
```

## Complete workflow example

```bash
# Step 1: Quality check (optional but recommended)
bash modules/1.1-quality_check_fastp.sh \
  --input_dir=raw_data/ \
  --output_dir=results/qc_reports \
  --nslots=16 \
  --overwrite=t

# Step 2: Check primers (optional but recommended)
python modules/1.2-check_primers.py \
  --input_dir raw_data/ \
  --output_dir results/primer_check \
  --suffix_r1 _R1_001.fastq.gz \
  --suffix_r2 _R2_001.fastq.gz \
  --primer_fwd GTGYCAGCMGCCGCGGTAA \
  --primer_rev CCGYCAATTYMTTTRAGTTT

# Step 3: Remove primers
bash modules/1.3-primer_removal_cutadapt.sh \
  --input_dir=raw_data/ \
  --output_dir=results/trimmed \
  --primer_fwd=GTGYCAGCMGCCGCGGTAA \
  --primer_rev=CCGYCAATTYMTTTRAGTTT \
  --nslots=16 \
  --overwrite=t

# Step 4: Generate ASVs
Rscript modules/2-dada2_pipeline.R \
  --input_dir results/trimmed/trimmed/ \
  --output_dir results/asvs \
  --nslots 16 \
  --overwrite

# Step 5: Annotate taxonomy
Rscript modules/3-taxa_annot.R \
  --input_asv_table results/asvs/asv_table.csv \
  --output_asv_table results/asv_table_annotated.csv \
  --method NBC \
  --train_db databases/silva_nr_v138_train_seq.fa.gz \
  --nslots 16
```

## Dependencies

- [Python 3](https://www.python.org/)
- [BioPython](https://biopython.org/)
- [R](https://www.r-project.org/)
- [fastp](https://github.com/OpenGene/fastp) (for quality control)
- [cutadapt](https://cutadapt.readthedocs.io/) (for primer removal)
- [DADA2](https://benjjneb.github.io/dada2/) R/Bioconductor package
- [tidyverse](https://www.tidyverse.org/) R package
- [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) R/Bioconductor package
- [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) R/Bioconductor package
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download) (for BLAST-based annotation)

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.





