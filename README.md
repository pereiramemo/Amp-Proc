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
│   ├── 1.2-quality-check.R                 # Comparative QC across samples (plots, PhiX)
│   ├── 1.3-primers-check.py                # IUPAC-aware primer check
│   ├── 1.4-primers-removal.py              # Primer removal with cutadapt
│   ├── 2.1-dada2-pipeline.R                # DADA2 ASV pipeline
│   ├── 2.2.1-vsearch-pipeline.py           # VSEARCH per-sample processing
│   ├── 2.2.2-vsearch-pipeline.py           # VSEARCH OTU clustering
│   ├── 3-taxa-annot.R                      # Taxonomic annotation
│   ├── utils.py                            # Shared Python helpers
│   └── utils.R                             # Shared R helpers
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
| `MODULE_1_2_QUALITY_CHECK`     | `1.2-quality-check.R`      | Comparative cross-sample QC: quality-vs-count plots, count histograms, PhiX (diagnostic, always runs) |
| `MODULE_1_3_PRIMERS_CHECK`     | `1.3-primers-check.py`     | IUPAC-aware primer detection (before & after trimming) |
| `MODULE_1_4_PRIMERS_REMOVAL`   | `1.4-primers-removal.py`   | cutadapt primer removal |
| `MODULE_2_1_DADA2_PIPELINE`    | `2.1-dada2-pipeline.R`     | DADA2 ASV inference (filter → denoise → merge → de-chimera) |
| `MODULE_2_2_1_VSEARCH_PIPELINE`| `2.2.1-vsearch-pipeline.py`| Per-sample merge → EE filter → derep → chimera check |
| `MODULE_2_2_2_VSEARCH_PIPELINE`| `2.2.2-vsearch-pipeline.py`| Pool samples → cluster OTUs → OTU table |
| `MODULE_3_TAXA_ANNOT`          | `3-taxa-annot.R`           | Taxonomy (NBC / NBCandEM) for ASVs and/or OTUs |

Each step writes a standardized layout under its publish directory: `output/`
(main results), `logs/` (a log file with a general-info header followed by any
third-party tool output), and `stats/` (TSV statistics). The output, logging, and
naming conventions are documented in `.claude/CLAUDE.md`.

## Workflow

After primer removal (`MODULE_1_4_PRIMERS_REMOVAL`), the `--method` parameter selects the
denoising branch:

- `dada2`   → `MODULE_2_1_DADA2_PIPELINE` (ASV table)
- `vsearch` → `MODULE_2_2_1_VSEARCH_PIPELINE` + `MODULE_2_2_2_VSEARCH_PIPELINE` (OTU table)
- `both`    → both branches in parallel (default)

`MODULE_1_1_QUALITY_CHECK` (fastp QC), `MODULE_1_2_QUALITY_CHECK` (comparative cross-sample QC)
and `MODULE_1_3_PRIMERS_CHECK` (primer check, before and after) are diagnostic and always run. Taxonomic annotation (`MODULE_3_TAXA_ANNOT`) runs on
the ASV table and/or the OTU table when `--skip_tax_annot false` is set. Both are
sequence-keyed count tables (the VSEARCH OTU table is relabelled by sequence via
`--relabel_self`), so the same `3-taxa-annot.R` script handles either unchanged.

## Run

```bash
# Quick test with the bundled data (default: both denoising branches + taxonomy)
nextflow run amp-proc.nf

# Choose a single branch
nextflow run amp-proc.nf --method vsearch

# Enable taxonomic annotation (SILVA v138.2 DBs auto-download to ~/.amp-proc/db/ if missing)
nextflow run amp-proc.nf --skip_tax_annot false \
    --train_db ~/.amp-proc/db/silva_nr99_v138.2_toGenus_trainset.fa.gz \
    --ref_db   ~/.amp-proc/db/silva_v138.2_assignSpecies.fa.gz

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
`3-taxa-annot.R` downloads the recognized SILVA v138.1/v138.2 file into `~/.amp-proc/db/` before
annotating (a missing but *unrecognized* filename is a fatal error). The download runs inside
the process container, so it needs network access at runtime — pre-populate `~/.amp-proc/db/`
to skip it.

## Outputs

Every module publishes into `--output_dir` under a directory named after the module
(`<module>-out/`), and each of those directories follows the same three-way layout:

- `output/` — the module's actual results
- `logs/` — a log file with a general-info header (date, sample, inputs, parameters,
  outputs, command, exit status) followed by the raw output of any third-party tool the
  module ran
- `stats/` — a tab-delimited table, samples as rows and statistics as columns

Per-sample modules insert a `<sample_name>/` level between the module directory and that
layout; modules that process all samples at once do not.

```
<output_dir>/
├── 1.1-quality-check-out/<sample>/{output,logs,stats}/
├── 1.2-quality-check-out/{output,logs,stats}/
├── 1.3-primers-check-before-out/<sample>/{output,logs,stats}/
├── 1.3-primers-check-after-out/<sample>/{output,logs,stats}/
├── 1.4-primers-removal-out/<sample>/{output,logs,stats}/
├── 2.1-dada2-pipeline-out/{output,logs,stats}/                # --method dada2 | both
├── 2.2.1-vsearch-pipeline-out/<sample>/{output,logs,stats}/   # --method vsearch | both
├── 2.2.2-vsearch-pipeline-out/{output,logs,stats}/            # --method vsearch | both
├── 3-taxa-annot-asv-out/{output,logs,stats}/                  # DADA2 branch, unless --skip_tax_annot true
└── 3-taxa-annot-otu-out/{output,logs,stats}/                  # VSEARCH branch, unless --skip_tax_annot true
```

> **Note on `--full_output false`.** It suppresses publication of modules 1.1, 1.2, 1.3,
> 1.4, 2.2.1 **and 2.1** — so `2.1-dada2-pipeline-out/` and with it `asv_table.csv` are
> not written either. Modules 2.2.2 and 3 always publish. Keep the default (`true`) unless
> you only want the OTU table and the annotated tables.

### Main results per module

| Module directory | Key files under `output/` | Contents |
|---|---|---|
| `1.1-quality-check-out/<sample>/` | `<sample>_fastp.html`, `<sample>_fastp.json` | fastp per-sample QC report — per-cycle quality, base composition, adapter content, duplication — in browser-readable and machine-readable form |
| `1.2-quality-check-out/` | `r1_mean_q_vs_nseq.png`, `r2_mean_q_vs_nseq.png`, `samples_hist.png`, `samples_hist_log.png`, `samples_perc_phix_barplot.png` | Cross-sample comparison: mean quality vs read count per sample (R1 and R2), read-count histograms (linear and log), and % PhiX per sample (bars above 0.005 % highlighted) |
| `1.3-primers-check-before-out/<sample>/`<br>`1.3-primers-check-after-out/<sample>/` | `<sample>_primer-check.tsv` | 4 × 4 matrix: rows are read/primer combinations (`FwdReads.FwdPrimer`, `RevReads.FwdPrimer`, `FwdReads.RevPrimer`, `RevReads.RevPrimer`), columns are orientations (`Forward`, `Complement`, `Reverse`, `RevComp`), values are the % of subsampled reads carrying that primer in that orientation (IUPAC-aware). Run before and after trimming — the "after" table should be ~0 everywhere |
| `1.4-primers-removal-out/<sample>/` | `<sample>_R1_trimmed.fastq.gz`, `<sample>_R2_trimmed.fastq.gz` | Primer-trimmed reads. These are the input to both denoising branches |
| `2.1-dada2-pipeline-out/` | `tables/asv_table.csv`, `filtered/<sample>_R{1,2}_filt.fastq.gz` | **ASV table** — rows are ASV sequences, columns are samples, values are read counts (the first header field is empty, since sequences are row names) — plus the quality-filtered reads DADA2 worked from |
| `2.2.1-vsearch-pipeline-out/<sample>/` | `01-merged/`, `02-filtered/`, `03-derep/`, `04-chimera-checked/` | Per-sample VSEARCH intermediates in pipeline order: merged pairs (FASTQ) → expected-error-filtered FASTA → dereplicated FASTA plus its `.uc` cluster map → chimera-checked FASTA plus the sequences flagged as chimeric |
| `2.2.2-vsearch-pipeline-out/` | `otu_table.tsv`, `otus.fasta.gz`, `all_samples.fasta.gz` | **OTU table** — the first column `#OTU ID` holds the centroid sequence (`--relabel_self`), the remaining columns are samples — plus the OTU centroid sequences and the pooled per-sample FASTA that was clustered |
| `3-taxa-annot-asv-out/`<br>`3-taxa-annot-otu-out/` | `tables/asv_table_annot.csv`<br>`tables/otu_table_annot.csv` | **The final annotated count table** — taxonomy, bootstrap support and per-sample counts in one file. Detailed below |

The `stats/` table of each module carries these columns (sample names in the first column):

| Stats file | Columns |
|---|---|
| `1.1-quality-check-<sample>-stats.tsv` | `total_reads`, `total_bases`, `q20_bases`, `q20_rate`, `q30_bases`, `q30_rate`, `read1_mean_length`, `read2_mean_length`, `gc_content` |
| `1.2-quality-check-stats.tsv` | `nseq`, `mean_q_r1`, `mean_q_r2`, `phix_pct` |
| `1.3-primers-check-<sample>-stats.tsv` | the four expected orientations, in %: `FwdReads.FwdPrimer.Forward`, `FwdReads.RevPrimer.RevComp`, `RevReads.RevPrimer.Forward`, `RevReads.FwdPrimer.RevComp` |
| `1.4-primers-removal-<sample>-stats.tsv` | `total_pairs`, `trimmed_pairs`, `percent_trimmed` |
| `2.1-dada2-pipeline-stats.tsv` | DADA2 read tracking: `raw`, `filtered`, `denoisedR1`, `denoisedR2`, `merged`, `nobim`, plus merged-length summaries `mean_length`, `sd_length`, `max_length`, `min_length` |
| `2.2.1-vsearch-pipeline-<sample>-stats.tsv` | `pairs_in`, `pairs_merged`, `percent_merged`, `reads_passed`, `percent_passed`, `seqs_unique`, `percent_unique`, `chimeras`, `nonchimeras`, `pct_chimeric_seqs`, `abund_in`, `abund_nonchimeric`, `pct_abund_retained` |
| `2.2.2-vsearch-pipeline-stats.tsv` | one `all_samples` row: `pooled_seqs`, `otus` |
| `3-taxa-annot-stats.tsv` | `n_asvs` (ASV branch) or `n_otus` (OTU branch), then `mean_<rank>_boot` / `sd_<rank>_boot` for phylum, class, order, family and genus, then `perc_spec_annot`. One row per sample (sequences with count > 0 in it) plus a pooled `all_samples` row |

### The final table: `<unit>_table_annot.csv`

This is the end product of the pipeline — taxonomy, classifier confidence and abundances
in a single file:

- `3-taxa-annot-asv-out/output/tables/asv_table_annot.csv` — DADA2 ASVs
- `3-taxa-annot-otu-out/output/tables/otu_table_annot.csv` — VSEARCH OTUs

> Both branches run the same `3-taxa-annot.R`, which takes the unit name from
> `--taxa_unit` (the module passes the branch label, `asv` or `otu`). That name drives
> both the filename and the key column, so the two files are identical in structure and
> differ only in whether the first column is called `asv` or `otu`. Code that reads either
> one should key off the first column by position rather than by name.

**File format.** Comma-separated with a single header line, written by R's `write.csv()`:
character fields are quoted with `"`, numeric fields are unquoted, missing values are
written as bare `NA`, line endings are `\n`, and there is no row-name column. One row per
ASV/OTU.

**Columns, in order:**

| Position | Column | Type | Description |
|---|---|---|---|
| 1 | `asv` or `otu` | string (DNA) | The full sequence, named for `--taxa_unit`. This is the row key — not an accession or an `ASV_1`-style ID — and it is the same string that keys `asv_table.csv` / `otu_table.tsv`, which is what lets one script annotate either table |
| 2–7 | `tax.Kingdom`, `tax.Phylum`, `tax.Class`, `tax.Order`, `tax.Family`, `tax.Genus` | string or `NA` | SILVA taxonomy from DADA2's naive Bayes classifier (`assignTaxonomy` against `--train_db`, bootstrap threshold `minBoot = 50`). A rank is `NA` for either of two reasons: bootstrap support below 50 %, **or** the SILVA lineage simply not naming that rank. Once a rank is `NA`, every rank below it is `NA` too |
| 8 | `tax.Species` | string or `NA` | **Only with `--taxa_method NBCandEM`.** Species from exact matching (`addSpecies` against `--ref_db`); `NA` when the sequence has no exact match in the reference. The column is absent under the default `NBC` |
| next 6 | `boot.Kingdom`, `boot.Phylum`, `boot.Class`, `boot.Order`, `boot.Family`, `boot.Genus` | integer, 0–100 | Bootstrap support for that rank, as a percentage of 100 replicates. Always read it next to the matching `tax.` column: a high `boot.` on an `NA` rank means the classifier was confident but the reference has no name there — do not treat the bootstrap value alone as evidence of an assignment. There is deliberately no `boot.Species`: exact matching is not bootstrapped |
| last *N* | one column per sample | integer | Read count of that sequence in that sample |

So the table is `1 + 6 + 6 + N` = **13 + N** columns wide with the default `--taxa_method NBC`,
and **14 + N** with `NBCandEM`, where *N* is the number of samples.

Sample column names come straight from the count table being annotated, i.e. the
`fromFilePairs` prefix of the input FASTQ files. Note that the VSEARCH branch replaces
`-` with `_` (VSEARCH requires sample names without hyphens), so a sample that appears as
`1-samo1_S1_L001` in the ASV table appears as `1_samo1_S1_L001` in the OTU table.

**Example** — the ASV branch (default `NBC`, three samples; sequences truncated for
display). The OTU branch is identical except the first column is headed `"otu"`:

```csv
"asv","tax.Kingdom","tax.Phylum","tax.Class","tax.Order","tax.Family","tax.Genus","boot.Kingdom","boot.Phylum","boot.Class","boot.Order","boot.Family","boot.Genus","1-samo1_S1_L001","2-samo2_S2_L001","3-samo3_S3_L001"
"TACGAAGGGACCTAGCGTAGTTCGG…CGCAAGGTTA","Bacteria","Pseudomonadota","Alphaproteobacteria","Pelagibacterales","Clade I","Clade Ia",100,100,100,100,100,100,918,59,1239
"TACGGAGGGGGTTAGCGTTGTTCGG…CGCAAGATTA","Bacteria","Pseudomonadota","Alphaproteobacteria","Rhodobacterales","Paracoccaceae","Amylibacter",100,100,100,100,100,100,5,1452,357
"TACCGGCAGCTCAAGTGGTCGTCGC…CGCAAGGCTG","Archaea","Thermoplasmatota","Thermoplasmata","Marine Group II",NA,NA,100,100,100,100,100,100,11,936,5
```

The third row shows the common `NA` case: SILVA has no family or genus name under
*Marine Group II*, so `tax.Family` and `tax.Genus` are `NA` even though the bootstrap
columns read 100. In the bundled three-sample test run, 325 of the 357 `NA` ranks are of
this kind and only 32 come from support below the 50 % threshold.

**Loading it.** Split the three column blocks by prefix, so the code does not depend on the
number of ranks or samples:

```r
annot  <- read.csv("3-taxa-annot-otu-out/output/tables/otu_table_annot.csv",
                   check.names = FALSE)          # keep sample names verbatim
key    <- names(annot)[1]                        # "asv" or "otu"
tax    <- annot[, grep("^tax\\.",  names(annot)), drop = FALSE]
boot   <- annot[, grep("^boot\\.", names(annot)), drop = FALSE]
counts <- annot[, setdiff(names(annot), c(key, names(tax), names(boot)))]
rownames(tax) <- rownames(counts) <- annot[[key]]  # ready for phyloseq
```

```python
import pandas as pd

annot  = pd.read_csv("3-taxa-annot-asv-out/output/tables/asv_table_annot.csv",
                     index_col=0)   # first column: "asv" or "otu"
tax    = annot.filter(regex=r"^tax\.")
boot   = annot.filter(regex=r"^boot\.")
counts = annot.drop(columns=tax.columns.union(boot.columns))
```

`check.names = FALSE` matters in R: sample names often start with a digit, and the default
would silently rename `1-samo1_S1_L001` to `X1.samo1_S1_L001`.

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

Primers (MODULE_1_3_PRIMERS_CHECK, MODULE_1_4_PRIMERS_REMOVAL):
  --primer_fwd      STR   Forward primer 5'->3' (default: GTGYCAGCMGCCGCGGTAA)
  --primer_rev      STR   Reverse primer 5'->3' (default: CCGYCAATTYMTTTRAGTTT)

MODULE_1_3_PRIMERS_CHECK — primer check:
  --subsample_size  INT   Reads to subsample per file (default: 1000)

MODULE_1_4_PRIMERS_REMOVAL — cutadapt primer removal:
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

## Threads and parallelism

Two parameters control CPU usage:

- `--nslots` — threads given to one tool invocation (a single task).
- `--maxForks` — how many per-sample tasks run at the same time.

Modules fall into two groups, handled differently:

| Module group | Modules | Threads used | Concurrency |
|---|---|---|---|
| Per-sample | quality check, primer check/removal, per-sample VSEARCH | `nslots` each | up to `maxForks` at once |
| Aggregate (all samples together) | `MODULE_2_1_DADA2_PIPELINE`, `MODULE_2_2_2_VSEARCH_PIPELINE`, `MODULE_3_TAXA_ANNOT` | `nslots * maxForks` | one at a time, alone |

This is enforced by CPU reservations in `nextflow.config`: the local-executor budget is
`executor.cpus = nslots * maxForks`; per-sample processes reserve `cpus = nslots` (so
`maxForks` of them fit within the budget), while the aggregate processes reserve the entire
budget. Because an aggregate task needs every CPU, Nextflow will not start it until the
other tasks finish and will not run anything beside it — so it uses all threads without
oversubscription. Those modules pass `${task.cpus}` as their `--nslots`, so the thread
count always matches the reservation.

> **Keep `nslots * maxForks` at or below the machine's physical core count.** The executor
> budget is pinned to that product, so a larger value oversubscribes the CPUs (slower,
> higher risk of running out of memory) — Nextflow will not clamp it for you. Example on a
> 48-core host: `nextflow run amp-proc.nf --nslots 16 --maxForks 3` (16 × 3 = 48).

With `--method both`, the DADA2 and VSEARCH-OTU stages therefore run sequentially, as do
the ASV and OTU taxonomy steps. This is deliberate: those steps are thread- and
memory-heavy, so running them one at a time with all threads is faster and avoids doubling
peak memory.

## Building & publishing the images

End users do **not** need this section — the published images pull automatically. It is
only for rebuilding and republishing after changing a Dockerfile or a pinned dependency
in `docker/resources/*.requirements.yml`. The images are built from the per-module
Dockerfiles in `docker/` by `docker/dockerbuild_commands.sh` (run from the repository
root):

```bash
# Build/tag/push the current module images
echo "$GHCR_PAT" | docker login ghcr.io -u pereiramemo --password-stdin   # PAT needs write:packages
bash docker/dockerbuild_commands.sh

# Build/tag/push plus an additional immutable version tag (:v1.0.0)
VERSION=v1.0.0 bash docker/dockerbuild_commands.sh
```

The script currently honours `VERSION` (adds an extra immutable tag alongside
`:latest`). At the moment it is configured to always push after building
(`PUSH=1` inside `docker/dockerbuild_commands.sh`). If you need build-only behavior,
set `PUSH` back to environment-driven logic in that script first. For air-gapped
hosts, `docker save`/`docker load` the images instead of pulling.

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
(`nextflow run amp-proc.nf --container_tag v1.0.0`) or by changing the default in
`nextflow.config`.

To cut a versioned release:

1. Build and push the versioned tag (this also updates `:latest`):
   ```bash
  VERSION=v1.0.0 bash docker/dockerbuild_commands.sh
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
