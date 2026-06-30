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
├── amp-proc.nf                             # Nextflow workflow entry point
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

Then clone the repository:

```bash
git clone https://github.com/pereiramemo/Amp-Proc.git
cd Amp-Proc
```

That is all that is required to run the pipeline. The per-module Docker images are
published publicly at `ghcr.io/pereiramemo/amp-proc/*`, and Nextflow pulls them
automatically on the first `nextflow run` (`docker.enabled = true` in `nextflow.config`).
You only need to build images yourself if you change a Dockerfile or a pinned dependency —
see [Building & publishing the images](#building--publishing-the-images).

For taxonomic annotation, the SILVA reference databases live under `~/.amp-proc/db/`
(the directory is mounted into the `MODULE_3_TAXA_ANNOT` container — see `nextflow.config`).
You do not have to download them manually: `3-taxa-annot.R` fetches the DADA2-formatted
SILVA files into that directory automatically on first use if they are missing. The
defaults are SILVA v138.2 — `silva_nr99_v138.2_toGenus_trainset.fa.gz` (NBC) and
`silva_v138.2_assignSpecies.fa.gz` (NBCandEM); the v138.1 files
(`silva_nr99_v138.1_train_set.fa.gz`, `silva_species_assignment_v138.1.fa.gz`) are also
recognized and can be selected via `--train_db`/`--ref_db`. Place them there beforehand
to avoid a runtime download (or if the container has no network access).

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
./amp-proc.nf

# Choose a single branch
./amp-proc.nf --method vsearch

# Enable taxonomic annotation (SILVA v138.2 DBs auto-download to ~/.amp-proc/db/ if missing)
./amp-proc.nf --skip_tax_annot false \
    --train_db ~/.amp-proc/db/silva_nr99_v138.2_toGenus_trainset.fa.gz \
    --ref_db   ~/.amp-proc/db/silva_v138.2_assignSpecies.fa.gz

# On your own data
./amp-proc.nf \
    --input_dir     /path/to/fastq \
    --reads_pattern '*_R{1,2}_001.fastq.gz' \
    --output_dir    /path/to/results \
    --primer_fwd    GTGYCAGCMGCCGCGGTAA \
    --primer_rev    CCGYCAATTYMTTTRAGTTT \
    --nslots        16

# Full parameter listing
./amp-proc.nf --help
```

Reference databases for `MODULE_3_TAXA_ANNOT` are mounted into the container from `~/.amp-proc`
(configured via `containerOptions` in `nextflow.config`). If a requested database is absent,
`3-taxa-annot.R` downloads the recognized SILVA v138.1/v138.2 file into `~/.amp-proc/db/` before
annotating (a missing but *unrecognized* filename is a fatal error). The download runs inside
the process container, so it needs network access at runtime — pre-populate `~/.amp-proc/db/`
to skip it.

## Parameters

All parameters have defaults in `nextflow.config` and can be overridden on the command
line (e.g. `--nslots 16`). The full list (output of `nextflow run amp-proc.nf --help`):

```text
Amp-Proc: amplicon processing from paired-end reads to ASV/OTU tables

Usage: nextflow run amp-proc.nf [options]

General:
  --input_dir       DIR   Input directory with paired-end FASTQ files (default: ./tests/data)
  --reads_pattern   STR   Glob pattern for fromFilePairs (default: *_R{1,2}_001_redu.fastq.gz)
  --output_dir      DIR   Output directory (default: ./tests/output_nf)
  --nslots          INT   CPU threads per tool (default: 12)
  --method          STR   Denoising branch: dada2 | vsearch | both (default: both)
  --full_output     BOOL  Publish all intermediate outputs (default: true)
  --skip_tax_annot  BOOL  Skip MODULE_3_TAXA_ANNOT taxonomic annotation (default: false)
  --maxForks        INT   Max parallel process instances (default: 3)
  --container_tag   STR   Tag of the ghcr.io/pereiramemo/amp-proc/* images to pull (default: latest)

Primers (MODULE_1_2_PRIMERS_CHECK, MODULE_1_3_PRIMERS_REMOVAL):
  --primer_fwd      STR   Forward primer 5'->3' (default: GTGYCAGCMGCCGCGGTAA)
  --primer_rev      STR   Reverse primer 5'->3' (default: CCGYCAATTYMTTTRAGTTT)

MODULE_1_2_PRIMERS_CHECK — primer check:
  --subsample_size  INT   Reads to subsample per file (default: 1000)

MODULE_1_3_PRIMERS_REMOVAL — cutadapt primer removal:
  --error_rate        NUM  Max allowed error rate (default: 0.1)
  --min_overlap       INT  Min primer-read overlap (default: 5)
  --min_length        INT  Discard reads shorter than this (default: 50)
  --discard_untrimmed STR  Discard reads with no primer, t/f (default: t)

MODULE_2_1_DADA2_PIPELINE — DADA2 (ASV):
  --trunc_r1          INT  Truncate R1 from 3' end (default: 250)
  --trunc_r2          INT  Truncate R2 from 3' end (default: 200)
  --dada2_min_overlap INT  Min overlap when merging (default: 12)
  --bimeras_method    STR  pooled | consensus | per-sample (default: consensus)

MODULE_2_2_1_VSEARCH_PIPELINE — VSEARCH per-sample:
  --fastq_minovlen     INT  Min overlap for PE merging (default: 5)
  --fastq_maxdiffs     INT  Max mismatches in overlap (default: 2)
  --fastq_maxee        NUM  Max expected errors per read (default: 1.0)
  --min_size           INT  Min abundance after derep (default: 1)
  --abskew             NUM  Min parent/child ratio, chimeras (default: 2.0)
  --vsearch_min_length INT  Min merged-read length (default: 50)

MODULE_2_2_2_VSEARCH_PIPELINE — VSEARCH OTU construction:
  --identity        NUM   OTU clustering identity 0-1 (default: 0.97)

MODULE_3_TAXA_ANNOT — taxonomic annotation (SILVA):
  --taxa_method     STR   NBC | NBCandEM (default: NBC)
  --train_db        PATH  NBC training database (default: $HOME/.amp-proc/db/silva_nr99_v138.2_toGenus_trainset.fa.gz)
  --ref_db          PATH  EM reference database (default: $HOME/.amp-proc/db/silva_v138.2_assignSpecies.fa.gz)
```

## Building & publishing the images

End users do **not** need this section — the published images pull automatically. It is
only for rebuilding and republishing after changing a Dockerfile or a pinned dependency
in `docker/resources/*.requirements.yml`. The images are built from the per-module
Dockerfiles in `docker/` by `docker/dockerbuild_commands.sh` (run from the repository
root):

```bash
# Build + tag :latest locally
bash docker/dockerbuild_commands.sh

# Build, tag with a version, and push (:latest and :v1.0.0) to the registry
echo "$GHCR_PAT" | docker login ghcr.io -u pereiramemo --password-stdin   # PAT needs write:packages
PUSH=1 VERSION=v1.0.0 bash docker/dockerbuild_commands.sh
```

The script honours two environment variables: `VERSION` (adds an extra immutable tag
alongside `:latest`) and `PUSH=1` (pushes after building). For air-gapped hosts,
`docker save`/`docker load` the images instead of pulling.

Newly pushed packages are **private by default**; make each one public (GitHub → your
profile → **Packages** → select the package → **Package settings** → **Change visibility**
→ **Public**) so machines can pull them anonymously. Otherwise every host must run
`docker login ghcr.io` before its first `nextflow run`.

### Reproducible installs (image version pinning)

Every module pulls `ghcr.io/pereiramemo/amp-proc/<module>:${params.container_tag}`.
`container_tag` is an ordinary pipeline parameter: it defaults to `latest` (always tracks
the most recent build) in `nextflow.config`, and like any parameter it can be overridden on
the command line with `--container_tag`. For a reproducible install — where a given checkout
always resolves to the same immutable image set — pin a published version, either per run
(`./amp-proc.nf --container_tag v1.0.0`) or by changing the default in
`nextflow.config`.

To cut a versioned release:

1. Build and push the versioned tag (this also updates `:latest`):
   ```bash
   PUSH=1 VERSION=v1.0.0 bash docker/dockerbuild_commands.sh
   ```
2. Pin `container_tag` to that version (in `nextflow.config` or via `--container_tag`).

A pinned tag **must already be published**, or the pull fails.

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
