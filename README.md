# Amplicon Processing Pipelines

This repository contains scripts for quality checking, preprocessing, denoising/clustering,
and taxonomic annotation of amplicon sequencing data. The analysis scripts live in `bin/`
and are written in Python and R: `1.1-quality-check.py`, `1.2-check-primers.py`,
`1.3-remove-primers.py`, `2.1-dada2-pipeline.R`, `2.2.1-vsearch-pipeline.py`,
`2.2.2-vsearch-pipeline.py`, `2.2.3-otu-to-seqtable.py`, and `3-taxa_annot.R`.

The pipeline can be run in two ways:

1. **Standalone scripts** — run each step directly from `bin/`, using a single conda/mamba
   environment (`environment.yml`). See [Installation](#installation) and [Scripts](#scripts).
2. **Nextflow workflow** — a containerized, end-to-end workflow (`main.nf`) where each step
   runs in its own Docker image. See [Nextflow pipeline](#nextflow-pipeline).

Two denoising/clustering strategies are available and can be used interchangeably:
**DADA2** (Amplicon Sequence Variants, ASVs) and **VSEARCH** (Operational Taxonomic
Units, OTUs).

## Repository structure

```
.
├── LICENSE                                # License file
├── README.md                              # This file
├── environment.yml                        # Conda/Mamba environment (standalone route)
├── main.nf                                # Nextflow workflow entry point
├── nextflow.config                        # Nextflow parameters and Docker settings
├── bin/                                   # Analysis scripts (Python + R)
│   ├── 1.1-quality-check.py               # Quality check with fastp
│   ├── 1.2-check-primers.py               # IUPAC-aware primer check
│   ├── 1.3-remove-primers.py             # Primer removal with cutadapt
│   ├── 2.1-dada2-pipeline.R             # DADA2 ASV pipeline
│   ├── 2.2.1-vsearch-pipeline.py       # VSEARCH per-sample processing
│   ├── 2.2.2-vsearch-pipeline.py       # VSEARCH OTU clustering
│   ├── 2.2.3-otu-to-seqtable.py        # OTU -> sequence-keyed table (for taxonomy)
│   ├── 3-taxa_annot.R                    # Taxonomic annotation
│   ├── toolbox.R                         # Shared R utility functions
│   └── old/                             # Legacy bash/R scripts (deprecated)
├── modules/                             # Nextflow process definitions (*.nf)
├── docker/                              # Per-module Dockerfiles + build script
│   ├── Dockerfile.module-*
│   ├── dockerbuild_commands.sh
│   └── resources/requirements-module-*.yml
└── tests/
    ├── data/                           # Test FASTQ files (3 samples)
    └── test_commands.sh                # End-to-end test runner
```

## Installation

Clone the repository:

```bash
git clone https://github.com/pereiramemo/Amp-Proc.git
cd Amp-Proc
```

### Environment setup with mamba (standalone route)

All standalone-script dependencies are installed via a single conda/mamba environment.
Install mamba if needed:

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

> For the Nextflow route you do **not** need this environment — see
> [Nextflow pipeline → Installation](#installation-1) for the Nextflow + Docker setup.

## Scripts

The preprocessing scripts (`1.1`, `1.2`, `1.3`, `2.2.1`) operate on **one paired-end
sample at a time** (`--reads1` / `--reads2`); run them in a loop over your samples (see
[Workflow examples](#complete-workflow-examples)). DADA2 (`2.1`) and OTU clustering
(`2.2.2`) operate on the whole dataset at once.

### 1.1-quality-check.py

Generates quality-control reports for a paired-end sample using fastp. Runs in
report-only mode — reads are not modified.

Output (under `--output_dir`):
- `reports/<sample>_fastp.html` — interactive HTML report
- `reports/<sample>_fastp.json` — JSON metrics
- `reports/<sample>_fastp.log` — processing log
- `summary_report.txt` — human-readable summary

```
Usage: 1.1-quality-check.py [options]

  --reads1 FILE                    R1 FASTQ file (required)
  --reads2 FILE                    R2 FASTQ file (required)
  -o, --output_dir DIR             Output directory (required)
  --nslots INT                     Threads [default=12]
  --min_length INT                 Min read length, reporting only [default=50]
  --qualified_quality_phred INT    Phred for a qualified base [default=20]
  --unqualified_percent_limit INT  Max % unqualified bases [default=40]
  --disable_adapter_trimming t|f   Disable adapter trimming [default=t]
  --html_report t|f                Generate HTML report [default=t]
  --json_report t|f                Keep JSON report [default=t]
  --overwrite t|f                  Overwrite existing output [default=f]
```

### 1.2-check-primers.py

Checks primer presence in a paired-end sample with IUPAC ambiguity-code support.
Subsamples reads and reports detection rates across all four primer orientations. Run it
both before and after primer removal to confirm trimming worked.

Output (under `--output_dir`):
- `<sample>_primer_check.csv` — primer hit percentages (or counts) per orientation
- `primer_sequences.csv` — primer sequences in all four orientations
- `summary_report.txt` — key primer-hit statistics

```
Usage: 1.2-check-primers.py [options]

  --reads1 FILE              R1 FASTQ file (required)
  --reads2 FILE              R2 FASTQ file (required)
  -o, --output_dir DIR       Output directory (required)
  --primer_fwd CHAR          Forward primer sequence (required)
  --primer_rev CHAR          Reverse primer sequence (required)
  -s, --subsample_size INT   Reads to subsample per file [default=1000]
  --raw_counts               Write raw counts instead of percentages
```

### 1.3-remove-primers.py

Removes primers from a paired-end sample using cutadapt. Handles IUPAC ambiguity codes
and optionally discards reads without detected primers.

Output:
- `<trimmed_dir>/<sample>_R1_trimmed.fastq.gz` — trimmed R1 reads
- `<trimmed_dir>/<sample>_R2_trimmed.fastq.gz` — trimmed R2 reads
- `<output_dir>/logs/<sample>_cutadapt.log` — per-sample cutadapt log
- `<output_dir>/stats/summary.tsv` — trimming statistics
- `<output_dir>/summary_report.txt` — human-readable summary

```
Usage: 1.3-remove-primers.py [options]

  --reads1 FILE              R1 FASTQ file (required)
  --reads2 FILE              R2 FASTQ file (required)
  -o, --output_dir DIR       Per-sample directory for logs/stats (required)
  --trimmed_dir DIR          Directory for trimmed FASTQ [default: {output_dir}/trimmed]
  --primer_fwd CHAR          Forward primer 5'->3' (required)
  --primer_rev CHAR          Reverse primer 5'->3' (required)
  --nslots INT               Threads [default=12]
  --error_rate NUM           Maximum error rate for primer matching [default=0.1]
  --min_overlap INT          Minimum primer-read overlap [default=3]
  --min_length INT           Discard reads shorter than this after trimming [default=50]
  --discard_untrimmed t|f    Discard reads where no primer was found [default=f]
  --compress t|f             Compress output with gzip [default=t]
  --overwrite t|f            Overwrite previous output [default=f]
```

### 2.1-dada2-pipeline.R

Processes primer-trimmed paired-end reads and generates ASVs with DADA2. Steps: quality
filtering → error learning → denoising → PE merging → chimera removal. Operates on a
directory of trimmed reads.

Output (under `--output_dir`):
- `asv_table.csv` — ASV abundance table (sequences × samples)
- `filtered/` — quality-filtered reads
- `plots/` — quality, error, read-retention, and length plots
- `tables/` — read counts and per-sample length statistics
- `session_info.txt` — parameters and package versions

```
Usage: 2.1-dada2-pipeline.R [options]

  --input_dir CHAR        Directory with primer-trimmed FASTQ files (required)
  --output_dir CHAR       Output directory (required)
  --nslots NUM            Threads [default=12]
  --trunc_r1 NUM          Truncate R1 reads to this length [default=250]
  --trunc_r2 NUM          Truncate R2 reads to this length [default=200]
  --pattern_r1 CHAR       R1 filename pattern [default=_L001_R1_001.fastq.gz]
  --pattern_r2 CHAR       R2 filename pattern [default=_L001_R2_001.fastq.gz]
  --min_overlap NUM       Minimum overlap for PE merging [default=12]
  --bimeras_method CHAR   Chimera method: pooled, consensus, per-sample [default=consensus]
  --pooled / --no_pooled                   Enable/disable pooled denoising [default=pooled]
  --qual_plot / --no_qual_plot             Generate quality plots [default=on]
  --err_plot / --no_err_plot               Generate error plots [default=on]
  --save_workspace / --no_save_workspace   Save .RData [default=on]
  --overwrite             Overwrite previous output [default=FALSE]
```

### 2.2.1-vsearch-pipeline.py

Per-sample VSEARCH preprocessing for the OTU branch. Steps: PE merging → expected-error
filtering → dereplication → de novo chimera detection. Each sample is processed
independently; the outputs are then pooled by `2.2.2`.

Output (under `--output_dir`, one directory per sample):
- `01-merged/<sample>-01-merged.fastq.gz`
- `02-filtered/<sample>-02-filtered.fasta.gz`
- `03-derep/<sample>-03-derep.fasta.gz`
- `04-chimera-checked/<sample>-04-chimera-checked.fasta.gz`
- `logs/`, `stats/01..04-*-summary.tsv`, `summary_report.txt`

```
Usage: 2.2.1-vsearch-pipeline.py [options]

  --reads1 FILE          R1 FASTQ file (required)
  --reads2 FILE          R2 FASTQ file (required)
  -o, --output_dir DIR   Per-sample output directory (required)
  --nslots INT           Threads [default=12]
  --min_length INT       Minimum merged-read length [default=50]
  --fastq_minovlen INT   Minimum overlap for PE merging [default=5]
  --fastq_maxdiffs INT   Maximum mismatches in overlap region [default=2]
  --fastq_maxee NUM      Maximum expected errors per merged read [default=1.0]
  --min_size INT         Minimum abundance to keep after dereplication [default=1]
  --abskew NUM           Minimum parent/child abundance ratio for chimeras [default=2.0]
  --overwrite t|f        Overwrite previous output [default=f]
```

### 2.2.2-vsearch-pipeline.py

Pools the per-sample chimera-checked sequences produced by `2.2.1`, clusters them into
OTUs with `vsearch --cluster_size`, and writes the OTU table. Because the sequences are
sample-labeled and already dereplicated per sample, no global dereplication or extra
read-mapping step is needed — `--otutabout` builds the per-sample OTU table directly.

Output (under `--output_dir`):
- `all_samples.fasta.gz` — pooled, sample-labeled sequences
- `otus/otus.fasta.gz` — OTU representative (centroid) sequences
- `otus/otu_table.tsv` — OTU abundance table (OTUs × samples)
- `stats/otu_summary.tsv`, `logs/cluster.log`, `summary_report.txt`

```
Usage: 2.2.2-vsearch-pipeline.py [options]

  --samples_dir DIR   Parent directory containing per-sample 2.2.1 outputs (required)
  --output_dir DIR    Output directory (required)
  --nslots INT        Threads [default=12]
  --identity NUM      OTU clustering identity threshold [default=0.97]
  --overwrite t|f     Overwrite previous output [default=f]
```

### 2.2.3-otu-to-seqtable.py

Bridges the VSEARCH OTU outputs into a DADA2-style, **sequence-keyed** count table so the
OTU centroids can be annotated with `3-taxa_annot.R` (which expects DNA sequences as row
identifiers). Maps each `OTU_N` to its centroid sequence and writes a table whose first
column is the sequence and whose remaining columns are the per-sample counts.

```
Usage: 2.2.3-otu-to-seqtable.py [options]

  --otus_fasta FILE   OTU centroid FASTA (gzip or plain) (required)
  --otu_table FILE    OTU count table from 2.2.2 (required)
  -o, --output FILE   Output sequence-keyed CSV (required)
```

### 3-taxa_annot.R

Assigns taxonomy to ASVs or OTUs. Supports three methods: Naive Bayes Classifier (NBC),
NBC combined with exact species matching (NBCandEM), or BLAST. Reference databases are not
included and must be supplied (e.g. SILVA).

Output:
- Annotated table (CSV) with taxonomy columns
  - NBC/NBCandEM: Kingdom → Genus (± Species) with bootstrap confidence values
  - BLAST: best-hit accession, organism, taxonomy path, percent identity

```
Usage: 3-taxa_annot.R [options]

  --input_asv_table CHAR    Sequence-keyed ASV/OTU table (required)
  --input_fasta CHAR        FASTA of sequences to annotate (generated from table if omitted)
  --output_asv_table CHAR   Output annotated table path (required)
  --method CHAR             NBC, NBCandEM, or BLAST [default=NBC]
  --evalue NUM              E-value threshold for BLAST [default=1e-10]
  --min_identity NUM        Minimum identity for BLAST [default=97]
  --train_db CHAR           Training database for NBC [default=silva_nr_v138_train_seq.fa.gz]
  --ref_db CHAR             Reference database for EM [default=silva_species_assignment_v138.fa.gz]
  --blast_db CHAR           BLAST database for the BLAST method
  --taxa_map CHAR           TSV mapping SILVA accessions to taxonomy (BLAST only)
  --nslots NUM              Threads [default=12]
  --save_workspace / --no_save_workspace   Save .RData [default=on]
  --overwrite               Overwrite previous output [default=FALSE]
```

> For the OTU branch, run `2.2.3-otu-to-seqtable.py` first and pass its output to
> `--input_asv_table`. The DADA2 `asv_table.csv` is already sequence-keyed and can be used
> directly.

## Complete workflow examples

The preprocessing steps run per sample. The examples below assume raw reads named
`<sample>_R1_001.fastq.gz` / `<sample>_R2_001.fastq.gz` in `data/raw/`.

### Workflow A: DADA2 (ASVs)

```bash
PRIMER_FWD="GTGYCAGCMGCCGCGGTAA"
PRIMER_REV="CCGYCAATTYMTTTRAGTTT"
INPUT_DIR="data/raw"
TRIMMED_DIR="results/03_primer_removal/trimmed"
NSLOTS=16

for R1 in "${INPUT_DIR}"/*_R1_001.fastq.gz; do
  SAMPLE=$(basename "${R1}" _R1_001.fastq.gz)
  R2="${INPUT_DIR}/${SAMPLE}_R2_001.fastq.gz"

  # Step 1: Quality check
  bin/1.1-quality-check.py \
    --reads1 "${R1}" --reads2 "${R2}" \
    --output_dir "results/01_qc/${SAMPLE}" \
    --nslots "${NSLOTS}" --overwrite t

  # Step 2: Check primers (before trimming)
  bin/1.2-check-primers.py \
    --reads1 "${R1}" --reads2 "${R2}" \
    --output_dir "results/02_primer_check/${SAMPLE}" \
    --primer_fwd "${PRIMER_FWD}" --primer_rev "${PRIMER_REV}"

  # Step 3: Remove primers
  bin/1.3-remove-primers.py \
    --reads1 "${R1}" --reads2 "${R2}" \
    --output_dir "results/03_primer_removal/${SAMPLE}" \
    --trimmed_dir "${TRIMMED_DIR}" \
    --primer_fwd "${PRIMER_FWD}" --primer_rev "${PRIMER_REV}" \
    --nslots "${NSLOTS}" --discard_untrimmed t --overwrite t
done

# Step 4: Generate ASVs with DADA2 (whole dataset)
bin/2.1-dada2-pipeline.R \
  --input_dir "${TRIMMED_DIR}" \
  --output_dir results/04_dada2 \
  --pattern_r1 _R1_trimmed.fastq.gz \
  --pattern_r2 _R2_trimmed.fastq.gz \
  --trunc_r1 250 --trunc_r2 200 \
  --nslots "${NSLOTS}" --overwrite

# Step 5: Annotate taxonomy
bin/3-taxa_annot.R \
  --input_asv_table results/04_dada2/asv_table.csv \
  --output_asv_table results/asv_table_annotated.csv \
  --method NBC \
  --train_db databases/silva_nr_v138_train_seq.fa.gz \
  --nslots "${NSLOTS}"
```

### Workflow B: VSEARCH (OTUs)

Steps 1–3 are identical to Workflow A. Replace steps 4–5 with:

```bash
# Step 4a: Per-sample VSEARCH processing
for R1 in "${TRIMMED_DIR}"/*_R1_trimmed.fastq.gz; do
  SAMPLE=$(basename "${R1}" _R1_trimmed.fastq.gz)
  R2="${TRIMMED_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"

  bin/2.2.1-vsearch-pipeline.py \
    --reads1 "${R1}" --reads2 "${R2}" \
    --output_dir "results/05_vsearch/${SAMPLE}" \
    --nslots "${NSLOTS}" --overwrite t
done

# Step 4b: Pool samples and cluster into OTUs
bin/2.2.2-vsearch-pipeline.py \
  --samples_dir results/05_vsearch \
  --output_dir results/05_vsearch/otu \
  --identity 0.97 --nslots "${NSLOTS}" --overwrite t

# Step 5a: Build a sequence-keyed table from the OTU centroids
bin/2.2.3-otu-to-seqtable.py \
  --otus_fasta results/05_vsearch/otu/otus/otus.fasta.gz \
  --otu_table  results/05_vsearch/otu/otus/otu_table.tsv \
  --output     results/05_vsearch/otu/otus/otu_seqtable.csv

# Step 5b: Annotate taxonomy
bin/3-taxa_annot.R \
  --input_asv_table results/05_vsearch/otu/otus/otu_seqtable.csv \
  --output_asv_table results/otu_table_annotated.csv \
  --method NBC \
  --train_db databases/silva_nr_v138_train_seq.fa.gz \
  --nslots "${NSLOTS}"
```

## Testing

An end-to-end test script is provided in `tests/`. It runs all steps in pipeline order
using the three-sample dataset in `tests/data/`.

```bash
bash tests/test_commands.sh
```

Taxonomic annotation is skipped automatically, as it requires external reference
databases.

## Dependencies

| Tool | Purpose |
|------|---------|
| [fastp](https://github.com/OpenGene/fastp) | Read quality control |
| [cutadapt](https://cutadapt.readthedocs.io/) | Primer removal |
| [vsearch](https://github.com/torognes/vsearch) | PE merging, dereplication, chimera checking, OTU clustering |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | BLAST-based taxonomic annotation |
| [Python 3](https://www.python.org/) + [BioPython](https://biopython.org/) | Primer checking, OTU table bridging |
| [R](https://www.r-project.org/) ≥ 4.0 | DADA2 pipeline and taxonomic annotation |
| [DADA2](https://benjjneb.github.io/dada2/) | ASV inference and NBC taxonomy |
| [tidyverse](https://www.tidyverse.org/) | Data wrangling in R |
| [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) | FASTQ handling in R |
| [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) | Sequence handling in R |

For the standalone route, all dependencies are available through the provided
`environment.yml`. For the Nextflow route they are pinned per module in
`docker/resources/requirements-module-*.yml`.

## Nextflow pipeline

In addition to the standalone scripts, the pipeline is available as a containerized
[Nextflow](https://www.nextflow.io/) workflow that mirrors the structure and
conventions of [Mg-Clust](https://github.com/pereiramemo/Mg-Clust). Each step is a
process defined under `modules/`, wrapping the corresponding script in `bin/`, and runs
in its own Docker image built from `docker/Dockerfile.module-*`.

### Installation

The Nextflow route does **not** require the conda/mamba environment described above —
every step runs inside its own Docker image, so the only prerequisites are:

- **Java** 11 or later (required by Nextflow)
- **[Nextflow](https://www.nextflow.io/)** 23.04 or later:
  ```bash
  curl -s https://get.nextflow.io | bash
  sudo mv nextflow /usr/local/bin/    # or any directory on your PATH
  ```
- **[Docker](https://docs.docker.com/get-docker/)** (the daemon must be running; the
  invoking user must be able to run `docker`)

Then clone the repository and build the per-module images:

```bash
git clone https://github.com/pereiramemo/Amp-Proc.git
cd Amp-Proc
bash docker/dockerbuild_commands.sh
```

This builds and tags one image per module (`ghcr.io/epereira/amp-proc/module-*`).
For taxonomic annotation, place the SILVA reference databases under `~/.amp-proc/db/`
(the directory is mounted into the `MODULE_3` container — see `nextflow.config`).

### Layout

```
.
├── main.nf                                 # Workflow entry point
├── nextflow.config                         # Parameters and Docker settings
├── bin/                                    # Step scripts (auto-staged onto PATH)
├── modules/
│   ├── 1.1-quality-check.nf                # MODULE_1_1   — fastp QC
│   ├── 1.2-check-primers.nf                # MODULE_1_2   — primer check (before/after)
│   ├── 1.3-remove-primers.nf              # MODULE_1_3   — cutadapt primer removal
│   ├── 2.1-dada2-pipeline.nf             # MODULE_2_1   — DADA2 ASV inference
│   ├── 2.2.1-vsearch-pipeline.nf        # MODULE_2_2_1 — VSEARCH per-sample
│   ├── 2.2.2-vsearch-pipeline.nf        # MODULE_2_2_2 — VSEARCH OTU clustering
│   ├── 2.2.3-otu-to-seqtable.nf         # MODULE_2_2_3 — OTU -> sequence-keyed table
│   └── 3-taxa-annot.nf                    # MODULE_3     — taxonomic annotation
└── docker/
    ├── Dockerfile.module-*                 # One image per module (micromamba)
    ├── dockerbuild_commands.sh             # Build + tag all images
    └── resources/requirements-module-*.yml # Conda environment per module
```

### Workflow

After primer removal (`MODULE_1_3`), the `--method` parameter selects the denoising
branch:

- `dada2`   → `MODULE_2_1` (ASV table)
- `vsearch` → `MODULE_2_2_1` + `MODULE_2_2_2` (OTU table)
- `both`    → both branches in parallel (default)

`MODULE_1_1` (fastp QC) and `MODULE_1_2` (primer check, before and after) are
diagnostic and always run. Taxonomic annotation (`MODULE_3`) runs on the ASV table
and/or the OTU centroids when `--skip_tax_annot false` is set; for the OTU branch,
`MODULE_2_2_3` first rebuilds a sequence-keyed count table from the OTU centroids so
the same `3-taxa_annot.R` script can be reused unchanged.

### Run

```bash
# Quick test with the bundled data (denoising only; taxonomy needs reference DBs)
nextflow run main.nf

# Choose a single branch
nextflow run main.nf --method vsearch

# Enable taxonomic annotation (place SILVA databases under ~/.amp-proc/db/)
nextflow run main.nf --skip_tax_annot false \
    --train_db ~/.amp-proc/db/silva_nr_v138_train_seq.fa.gz \
    --ref_db   ~/.amp-proc/db/silva_species_assignment_v138.fa.gz

# On your own data
nextflow run main.nf \
    --input_dir     /path/to/fastq \
    --reads_pattern '*_R{1,2}_001.fastq.gz' \
    --output_dir    /path/to/results \
    --primer_fwd    GTGYCAGCMGCCGCGGTAA \
    --primer_rev    CCGYCAATTYMTTTRAGTTT \
    --nslots        16

# Full parameter listing
nextflow run main.nf --help
```

Reference databases for `MODULE_3` are mounted into the container from `~/.amp-proc`
(configured via `containerOptions` in `nextflow.config`).

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
