# Amplicon Processing Pipelines

This repository provides a containerized [Nextflow](https://www.nextflow.io/)
pipeline for quality checking, preprocessing, denoising/clustering, and taxonomic
annotation of amplicon sequencing data. Each step is a Nextflow process (under
`modules/`) that wraps a Python or R script in `bin/` and runs in its own Docker
image.

Two denoising/clustering strategies are available and can be used interchangeably
(or together), selected with `--method`: **DADA2** (Amplicon Sequence Variants,
ASVs) and **VSEARCH** (Operational Taxonomic Units, OTUs).

## Repository structure

```
.
├── LICENSE                                 # License file
├── README.md                               # This file
├── amp-proc.nf                                 # Nextflow workflow entry point
├── nextflow.config                         # Nextflow parameters and Docker settings
├── bin/                                    # Step scripts (auto-staged onto PATH)
│   ├── 1.1-quality-check.py                # Quality check with fastp
│   ├── 1.2-primers-check.py                # IUPAC-aware primer check
│   ├── 1.3-primers-removal.py              # Primer removal with cutadapt
│   ├── 2.1-dada2-pipeline.R                # DADA2 ASV pipeline
│   ├── 2.2.1-vsearch-pipeline.py           # VSEARCH per-sample processing
│   ├── 2.2.2-vsearch-pipeline.py           # VSEARCH OTU clustering
│   ├── 3-taxa-annot.R                      # Taxonomic annotation
│   ├── toolbox.py                          # Shared Python helpers
│   └── toolbox.R                           # Shared R helpers
├── modules/                                # Nextflow process definitions (*.nf)
├── docker/                                 # Per-module Dockerfiles + build script
    ├── *.Dockerfile
    ├── dockerbuild_commands.sh
    └── resources/*.requirements.yml

```

## Installation

The pipeline runs entirely in containers, so the only prerequisites are:

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

This builds and tags one image per module (`ghcr.io/epereira/amp-proc/*`).
For taxonomic annotation, the SILVA reference databases live under `~/.amp-proc/db/`
(the directory is mounted into the `MODULE_3_TAXA_ANNOT` container — see `nextflow.config`).
You do not have to download them manually: `3-taxa-annot.R` fetches the DADA2-formatted
SILVA v138.1 files into that directory automatically on first use if they are missing —
`silva_nr99_v138.1_train_set.fa.gz` (NBC) and `silva_species_assignment_v138.1.fa.gz`
(NBCandEM). Place them there beforehand to avoid a runtime download (or if the container
has no network access).

## Pipeline steps

| Module | Script | Purpose |
|--------|--------|---------|
| `MODULE_1_1_QUALITY_CHECK`     | `1.1-quality-check.py`     | fastp QC report (report-only; diagnostic, always runs) |
| `MODULE_1_2_PRIMERS_CHECK`     | `1.2-primers-check.py`     | IUPAC-aware primer detection (before & after trimming) |
| `MODULE_1_3_PRIMERS_REMOVAL`   | `1.3-primers-removal.py`   | cutadapt primer removal |
| `MODULE_2_1_DADA2_PIPELINE`    | `2.1-dada2-pipeline.R`     | DADA2 ASV inference (filter → denoise → merge → de-chimera) |
| `MODULE_2_2_1_VSEARCH_PIPELINE`| `2.2.1-vsearch-pipeline.py`| Per-sample merge → EE filter → derep → chimera check |
| `MODULE_2_2_2_VSEARCH_PIPELINE`| `2.2.2-vsearch-pipeline.py`| Pool samples → cluster OTUs → OTU table |
| `MODULE_3_TAXA_ANNOT`          | `3-taxa-annot.R`           | Taxonomy (NBC / NBCandEM) for ASVs and/or OTUs |

Each step writes a standardized layout under its publish directory: `output/`
(main results), `logs/` (a log file with a general-info header followed by any
third-party tool output), and `stats/` (TSV statistics). The output, logging, and
naming conventions are documented in `.claude/CLAUDE.md`.

## Workflow

After primer removal (`MODULE_1_3_PRIMERS_REMOVAL`), the `--method` parameter selects the
denoising branch:

- `dada2`   → `MODULE_2_1_DADA2_PIPELINE` (ASV table)
- `vsearch` → `MODULE_2_2_1_VSEARCH_PIPELINE` + `MODULE_2_2_2_VSEARCH_PIPELINE` (OTU table)
- `both`    → both branches in parallel (default)

`MODULE_1_1_QUALITY_CHECK` (fastp QC) and `MODULE_1_2_PRIMERS_CHECK` (primer check, before and
after) are diagnostic and always run. Taxonomic annotation (`MODULE_3_TAXA_ANNOT`) runs on
the ASV table and/or the OTU table when `--skip_tax_annot false` is set. Both are
sequence-keyed count tables (the VSEARCH OTU table is relabelled by sequence via
`--relabel_self`), so the same `3-taxa-annot.R` script handles either unchanged.

## Run

```bash
# Quick test with the bundled data (denoising only; taxonomy needs reference DBs)
nextflow run amp-proc.nf

# Choose a single branch
nextflow run amp-proc.nf --method vsearch

# Enable taxonomic annotation (SILVA v138.1 DBs auto-download to ~/.amp-proc/db/ if missing)
nextflow run amp-proc.nf --skip_tax_annot false \
    --train_db ~/.amp-proc/db/silva_nr99_v138.1_train_set.fa.gz \
    --ref_db   ~/.amp-proc/db/silva_species_assignment_v138.1.fa.gz

# On your own data
nextflow run amp-proc.nf \
    --input_dir     /path/to/fastq \
    --reads_pattern '*_R{1,2}_001.fastq.gz' \
    --output_dir    /path/to/results \
    --primer_fwd    GTGYCAGCMGCCGCGGTAA \
    --primer_rev    CCGYCAATTYMTTTRAGTTT \
    --nslots        16

# Full parameter listing
nextflow run amp-proc.nf --help
```

Reference databases for `MODULE_3_TAXA_ANNOT` are mounted into the container from `~/.amp-proc`
(configured via `containerOptions` in `nextflow.config`). If a requested database is absent,
`3-taxa-annot.R` downloads the recognized SILVA v138.1 file into `~/.amp-proc/db/` before
annotating (a missing but *unrecognized* filename is a fatal error). The download runs inside
the process container, so it needs network access at runtime — pre-populate `~/.amp-proc/db/`
to skip it.

## Dependencies

Dependencies are pinned per module in `docker/resources/*.requirements.yml` and built
into the per-module images by `docker/dockerbuild_commands.sh` — there is nothing to
install manually beyond Nextflow and Docker (see [Installation](#installation)).

| Tool | Purpose |
|------|---------|
| [fastp](https://github.com/OpenGene/fastp) | Read quality control |
| [cutadapt](https://cutadapt.readthedocs.io/) | Primer removal |
| [vsearch](https://github.com/torognes/vsearch) | PE merging, dereplication, chimera checking, OTU clustering |
| [Python 3](https://www.python.org/) + [BioPython](https://biopython.org/) | Primer checking, OTU table bridging |
| [R](https://www.r-project.org/) ≥ 4.0 | DADA2 pipeline and taxonomic annotation |
| [DADA2](https://benjjneb.github.io/dada2/) | ASV inference and NBC taxonomy |
| [tidyverse](https://www.tidyverse.org/) | Data wrangling in R |
| [ShortRead](https://bioconductor.org/packages/release/bioc/html/ShortRead.html) | FASTQ handling in R |
| [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) | Sequence handling in R |

## License

This project is licensed under the GNU General Public License v3.0 - see the LICENSE file for details.

Copyright (C) 2025 Emiliano Pereira

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
