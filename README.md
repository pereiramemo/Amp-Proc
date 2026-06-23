# Amplicon Processing Pipelines

This repository contains scripts for quality checking, preprocessing, clustering, and taxonomic annotation of amplicon sequencing data. The pipeline is written in Bash, Python, and R and includes: `1.1-quality_check_fastp.sh`, `1.2-check_primers.py`, `1.3-primer_removal_cutadapt.sh`, `2.1-dada2_pipeline.R`, `2.2-vsearch_pipeline.sh`, and `3-taxa_annot.R`.

## Repository structure

```
.
├── LICENSE                                # License file
├── README.md                              # This file
├── environment.yml                        # Conda/Mamba environment specification
├── modules/                               # Main analysis scripts
│   ├── 1.1-quality_check_fastp.sh         # Quality check with fastp
│   ├── 1.2-check_primers.py               # Primer checking utility
│   ├── 1.3-primer_removal_cutadapt.sh     # Primer removal with cutadapt
│   ├── 2.1-dada2_pipeline.R               # DADA2 ASV generation pipeline
│   ├── 2.2-vsearch_pipeline.sh            # VSEARCH OTU generation pipeline
│   ├── 3-taxa_annot.R                     # Taxonomic annotation pipeline
│   ├── toolbox.R                          # Utility functions
│   └── conf.sh                            # Configuration file
└── tests/
    ├── data/                              # Test FASTQ files (3 samples)
    └── test_commands.sh                   # End-to-end test runner
```

## Installation

Clone the repository:

```bash
git clone https://github.com/pereiramemo/Amp-Proc.git
cd Amp-Proc
```

### Environment setup with mamba

All dependencies are installed via a single conda/mamba environment. Install mamba if needed:

- **Miniforge (recommended)**: https://github.com/conda-forge/miniforge
- Or add mamba to an existing conda base:
  ```bash
  conda install -n base -c conda-forge mamba
  ```

Create and activate the environment:

```bash
mamba env create -f environment.yml
conda activate amp-proc-env
```

## Scripts

### 1.1-quality_check_fastp.sh

Generates quality control reports for paired-end FASTQ files using fastp. Runs in report-only mode — reads are not modified.

Output:
- `reports/<sample>_fastp.html` — interactive HTML report
- `reports/<sample>_fastp.json` — JSON metrics
- `reports/<sample>_fastp.log` — processing log
- `stats/summary.tsv` — per-sample statistics table
- `summary_report.txt` — human-readable summary

```
Usage: 1.1-quality_check_fastp.sh <options>

Options:
    --help
    --input_dir=CHAR                 Directory containing input FASTQ files (required)
    --output_dir=CHAR                Output directory for reports (required)
    --pattern_r1=CHAR                R1 filename pattern [default=_R1_001.fastq.gz]
    --pattern_r2=CHAR                R2 filename pattern [default=_R2_001.fastq.gz]
    --nslots=NUM                     Threads [default=12]
    --min_length=NUM                 Minimum read length for reporting [default=50]
    --qualified_quality_phred=NUM    Minimum base quality (Phred) [default=20]
    --unqualified_percent_limit=NUM  Max % unqualified bases [default=40]
    --disable_adapter_trimming=t|f   Disable adapter trimming in report [default=t]
    --html_report=t|f                Generate HTML report [default=t]
    --json_report=t|f                Generate JSON report [default=t]
    --overwrite=t|f                  Overwrite previous output [default=f]
```

### 1.2-check_primers.py

Checks primer presence in FASTQ files with IUPAC ambiguity code support. Subsamples reads and reports detection rates across all primer orientations.

Output:
- `<sample>_primer_check.csv` — primer hit percentages (or counts) per orientation
- `primer_sequences.csv` — primer sequences in all four orientations

```
usage: check_primers.py [-h] -i INPUT_DIR -o OUTPUT_DIR --suffix_r1 SUFFIX_R1
                        --suffix_r2 SUFFIX_R2 --primer_fwd PRIMER_FWD
                        --primer_rev PRIMER_REV [-s SUBSAMPLE_SIZE] [--counts]

options:
  -i, --input_dir       Path to input directory (required)
  -o, --output_dir      Path to output directory (required)
  --suffix_r1           Filename suffix for R1 reads (required)
  --suffix_r2           Filename suffix for R2 reads (required)
  --primer_fwd          Forward primer sequence (required)
  --primer_rev          Reverse primer sequence (required)
  -s, --subsample_size  Reads to subsample per file [default=1000]
  --counts              Write raw counts instead of percentages
```

### 1.3-primer_removal_cutadapt.sh

Removes primer sequences from paired-end FASTQ files using cutadapt. Handles IUPAC ambiguity codes and optionally discards reads without detected primers.

Output:
- `trimmed/<sample>_R1_trimmed.fastq.gz` — trimmed R1 reads
- `trimmed/<sample>_R2_trimmed.fastq.gz` — trimmed R2 reads
- `logs/<sample>_cutadapt.log` — per-sample cutadapt log
- `stats/summary.tsv` — trimming statistics table
- `summary_report.txt` — human-readable summary

```
Usage: 1.3-primer_removal_cutadapt.sh <options>

Options:
    --help
    --input_dir=CHAR          Directory containing input FASTQ files (required)
    --output_dir=CHAR         Output directory for trimmed reads (required)
    --pattern_r1=CHAR         R1 filename pattern [default=_R1_001.fastq.gz]
    --pattern_r2=CHAR         R2 filename pattern [default=_R2_001.fastq.gz]
    --primer_fwd=CHAR         Forward primer sequence 5'→3' (required)
    --primer_rev=CHAR         Reverse primer sequence 5'→3' (required)
    --nslots=NUM              Threads [default=12]
    --error_rate=NUM          Maximum error rate for primer matching [default=0.1]
    --min_overlap=NUM         Minimum primer–read overlap [default=3]
    --min_length=NUM          Discard reads shorter than this after trimming [default=50]
    --discard_untrimmed=t|f   Discard reads where no primer was found [default=f]
    --compress=t|f            Compress output with gzip [default=t]
    --overwrite=t|f           Overwrite previous output [default=f]
```

### 2.1-dada2_pipeline.R

Processes primer-trimmed paired-end reads and generates Amplicon Sequence Variants (ASVs) using DADA2. Steps: quality filtering → error learning → denoising → PE merging → chimera removal.

Output:
- `asv_table.csv` — ASV abundance table (sequences × samples)
- `track_n_seqs.csv` — read counts at each processing step
- `plots/quality_plot_R1.pdf`, `plots/quality_plot_R2.pdf` — read quality profiles
- `plots/error_plot_R1.pdf`, `plots/error_plot_R2.pdf` — error model plots
- `plots/nseq_barplot.pdf` — read retention barplot across steps
- `plots/seq_length_hist.pdf` — merged read length distribution
- `tables/seq_length_stats.csv` — per-sample length statistics
- `session_info_dada2.txt` — parameters and package versions

```
Usage: ./2.1-dada2_pipeline.R <options>

Options:
    --help
    --input_dir CHAR        Directory with primer-trimmed FASTQ files (required)
    --output_dir CHAR       Output directory (required)
    --nslots NUM            Threads [default=12]
    --trunc_r1 NUM          Truncate R1 reads to this length [default=250]
    --trunc_r2 NUM          Truncate R2 reads to this length [default=200]
    --pattern_r1 CHAR       R1 filename pattern [default=_L001_R1_001.fastq.gz]
    --pattern_r2 CHAR       R2 filename pattern [default=_L001_R2_001.fastq.gz]
    --min_overlap NUM       Minimum overlap for PE merging [default=12]
    --bimeras_method CHAR   Chimera method: pooled, consensus, per-sample [default=consensus]
    --pooled / --no_pooled  Enable/disable pooled denoising [default=pooled]
    --qual_plot / --no_qual_plot   Generate quality plots [default=on]
    --err_plot / --no_err_plot     Generate error plots [default=on]
    --save_workspace / --no_save_workspace   Save .RData [default=on]
    --overwrite             Overwrite previous output [default=FALSE]
```

### 2.2-vsearch_pipeline.sh

Processes primer-trimmed paired-end reads and generates OTUs using VSEARCH. Steps: fastp quality filtering → vsearch PE merging → expected-error filtering → global dereplication → chimera removal → OTU clustering → read mapping → OTU table.

Output:
- `01_fastp/` — quality-filtered reads per sample
- `02_merged/` — merged FASTQ per sample
- `03_filtered/` — per-sample FASTA (quality filtered, sample-labeled) and concatenated `all_samples.fasta`
- `04_derep/derep.fasta` — globally dereplicated sequences
- `05_chimera/nochimeras.fasta` — chimera-filtered sequences
- `06_otus/otus.fasta` — OTU representative sequences
- `06_otus/otu_table.tsv` — OTU abundance table (OTUs × samples)
- `stats/fastp_summary.tsv`, `stats/merge_summary.tsv` — per-step statistics
- `summary_report.txt` — human-readable summary

```
Usage: 2.2-vsearch_pipeline.sh <options>

Options:
    --help
    --input_dir=CHAR          Directory with primer-trimmed FASTQ files (required)
    --output_dir=CHAR         Output directory (required)
    --pattern_r1=CHAR         R1 filename pattern [default=_R1_001.fastq.gz]
    --pattern_r2=CHAR         R2 filename pattern [default=_R2_001.fastq.gz]
    --nslots=NUM              Threads [default=12]
    --min_quality=NUM         Minimum base Phred quality (fastp) [default=20]
    --min_length=NUM          Minimum read length after trimming and after merging [default=50]
    --fastq_minovlen=NUM      Minimum overlap for PE merging [default=20]
    --fastq_maxdiffs=NUM      Maximum mismatches in overlap region [default=5]
    --fastq_maxee=NUM         Maximum expected errors per merged read [default=1.0]
    --min_size=NUM            Minimum abundance after dereplication [default=2]
    --chimera_method=CHAR     Chimera detection: denovo or ref [default=denovo]
    --ref_db=CHAR             Reference FASTA for chimera checking (required if chimera_method=ref)
    --identity=NUM            OTU clustering identity threshold [default=0.97]
    --overwrite=t|f           Overwrite previous output [default=f]
```

### 3-taxa_annot.R

Assigns taxonomy to ASVs or OTUs. Supports three annotation methods: Naive Bayes Classifier (NBC), NBC combined with exact species matching (NBCandEM), or BLAST.

Output:
- Annotated ASV/OTU table (CSV) with taxonomy columns
  - NBC/NBCandEM: Kingdom → Genus (± Species) with bootstrap confidence values
  - BLAST: best hit accession, organism name, taxonomy path, percent identity
- `session_info_taxa_annot.txt` — parameters and package versions

```
Usage: ./3-taxa_annot.R <options>

Options:
    --help
    --input_asv_table CHAR    ASV/OTU table from DADA2 or VSEARCH (required)
    --input_fasta CHAR        FASTA of sequences to annotate (generated from table if omitted)
    --output_asv_table CHAR   Output annotated table path (required)
    --method CHAR             NBC, NBCandEM, or BLAST [default=NBC]
    --evalue NUM              E-value threshold for BLAST [default=1e-10]
    --min_identity NUM        Minimum identity for BLAST [default=97]
    --train_db CHAR           Training database for NBC [default=silva_nr_v138_train_seq.fa.gz]
    --ref_db CHAR             Reference database for EM [default=silva_species_assignment_v138.fa.gz]
    --blast_db CHAR           BLAST database for BLAST method [default=SILVA_138_SSURef_NR99_tax_silva.fasta]
    --taxa_map CHAR           TSV mapping SILVA accessions to taxonomy (BLAST only)
    --nslots NUM              Threads [default=12]
    --save_workspace / --no_save_workspace   Save .RData [default=on]
    --overwrite               Overwrite previous output [default=FALSE]
```

## Complete workflow examples

### Workflow A: DADA2 (ASVs)

```bash
PRIMER_FWD="GTGYCAGCMGCCGCGGTAA"
PRIMER_REV="CCGTCAATTCMTTTRAGTTT"
INPUT_DIR="data/raw/"

# Step 1: Quality check
bash modules/1.1-quality_check_fastp.sh \
  --input_dir="${INPUT_DIR}" \
  --output_dir=results/01_qc \
  --nslots=16 \
  --overwrite=t

# Step 2: Check primers
python modules/1.2-check_primers.py \
  --input_dir "${INPUT_DIR}" \
  --output_dir results/02_primer_check \
  --suffix_r1 _R1_001.fastq.gz \
  --suffix_r2 _R2_001.fastq.gz \
  --primer_fwd "${PRIMER_FWD}" \
  --primer_rev "${PRIMER_REV}"

# Step 3: Remove primers
bash modules/1.3-primer_removal_cutadapt.sh \
  --input_dir="${INPUT_DIR}" \
  --output_dir=results/03_cutadapt \
  --primer_fwd="${PRIMER_FWD}" \
  --primer_rev="${PRIMER_REV}" \
  --nslots=16 \
  --discard_untrimmed=t \
  --overwrite=t

# Step 4a: Generate ASVs with DADA2
Rscript modules/2.1-dada2_pipeline.R \
  --input_dir results/03_cutadapt/trimmed/ \
  --output_dir results/04_dada2 \
  --pattern_r1 _R1_trimmed.fastq.gz \
  --pattern_r2 _R2_trimmed.fastq.gz \
  --trunc_r1 250 \
  --trunc_r2 200 \
  --nslots 16 \
  --overwrite

# Step 5: Annotate taxonomy
Rscript modules/3-taxa_annot.R \
  --input_asv_table results/04_dada2/asv_table.csv \
  --output_asv_table results/asv_table_annotated.csv \
  --method NBC \
  --train_db databases/silva_nr_v138_train_seq.fa.gz \
  --nslots 16
```

### Workflow B: VSEARCH (OTUs)

Steps 1–3 are identical to Workflow A. Replace step 4 with:

```bash
# Step 4b: Generate OTUs with VSEARCH
bash modules/2.2-vsearch_pipeline.sh \
  --input_dir=results/03_cutadapt/trimmed/ \
  --output_dir=results/04_vsearch \
  --pattern_r1=_R1_trimmed.fastq.gz \
  --pattern_r2=_R2_trimmed.fastq.gz \
  --chimera_method=denovo \
  --identity=0.97 \
  --nslots=16 \
  --overwrite=t

# Step 5: Annotate taxonomy
Rscript modules/3-taxa_annot.R \
  --input_asv_table results/04_vsearch/06_otus/otu_table.tsv \
  --output_asv_table results/otu_table_annotated.csv \
  --method NBC \
  --train_db databases/silva_nr_v138_train_seq.fa.gz \
  --nslots 16
```

## Testing

An end-to-end test script is provided in `tests/`. It runs all scripts in pipeline order using the three-sample dataset in `tests/data/` and reports pass/fail for each step.

```bash
./tests/test_commands.sh --clean
```

The `--clean` flag removes any previous test output before running. Taxa annotation (step 6) is skipped automatically as it requires external reference databases.

## Dependencies

| Tool | Purpose |
|------|---------|
| [fastp](https://github.com/OpenGene/fastp) | Read quality control and filtering |
| [cutadapt](https://cutadapt.readthedocs.io/) | Primer removal |
| [vsearch](https://github.com/torognes/vsearch) | PE merging, dereplication, chimera checking, OTU clustering |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | BLAST-based taxonomic annotation |
| [Python 3](https://www.python.org/) + [BioPython](https://biopython.org/) | Primer checking |
| [R](https://www.r-project.org/) ≥ 4.0 | DADA2 pipeline and taxonomic annotation |
| [DADA2](https://benjjneb.github.io/dada2/) | ASV inference |
| [tidyverse](https://www.tidyverse.org/) | Data wrangling in R |
| [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) | FASTQ handling in R |
| [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) | Sequence handling in R |

All dependencies are available through the provided `environment.yml`.

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
