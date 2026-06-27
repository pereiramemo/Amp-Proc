#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# DEV ONLY — comment out before production use
""" 
reads1       = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/03_primer_removal/trimmed/1-samo1_S1_L001_R1_trimmed.fastq.gz"
reads2       = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/03_primer_removal/trimmed/1-samo1_S1_L001_R2_trimmed.fastq.gz"
output_dir   = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch/1-samo1_S1_L001"
filtered_dir = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch/filtered"
nslots       = 4
overwrite    = True
min_length   = 50
minovlen     = 5
maxdiffs     = 2
maxee        = 1.0
min_size = 1
 """
################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Per-sample vsearch pre-processing: PE merging → EE filtering → dereplication"
    )
    p.add_argument("--reads1",          required=True,              help="R1 FASTQ file (required)")
    p.add_argument("--reads2",          required=True,              help="R2 FASTQ file (required)")
    p.add_argument("-o", "--output_dir", required=True,             help="Per-sample output directory for logs and stats (required)")
    p.add_argument("--sample_name",     default=None,               help="Sample name [default: derived from R1 filename]")
    p.add_argument("--nslots",          type=int,   default=12,     help="Number of threads [default=12]")
    p.add_argument("--min_length",      type=int,   default=50,     help="Minimum merged-read length [default=50]")
    p.add_argument("--fastq_minovlen",  type=int,   default=5,      help="Minimum overlap for PE merging [default=5]")
    p.add_argument("--fastq_maxdiffs",  type=int,   default=2,      help="Maximum mismatches in overlap region [default=2]")
    p.add_argument("--fastq_maxee",     type=float, default=1.0,    help="Maximum expected errors per merged read [default=1.0]")
    p.add_argument("--min_size",        type=int,   default=1,      help="Minimum abundance to keep after dereplication [default=1]")
    p.add_argument("--abskew",          type=float, default=2.0,    help="Minimum abundance ratio parent/child for chimera detection [default=2.0]")
    p.add_argument("--overwrite",       choices=["t", "f"], default="f", help="Overwrite existing per-sample output directory [default=f]")
    return p.parse_args()


def log(msg):
    print(f"[INFO] {msg}")

def log_warn(msg):
    print(f"\033[1;33m[WARN]\033[0m {msg}", file=sys.stderr)

def log_error(msg):
    print(f"\033[0;31m[ERROR]\033[0m {msg}", file=sys.stderr)

def run(cmd, log_path=None, check=True):
    result = subprocess.run(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    if log_path:
        Path(log_path).write_text(result.stderr)
    if check and result.returncode != 0:
        log_error(f"Command failed: {' '.join(str(c) for c in cmd)}")
        if result.stderr:
            print(result.stderr, file=sys.stderr)
        sys.exit(1)
    return result

def compress(src: Path, dst: Path):
    with open(src, "rb") as f_in, gzip.open(dst, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    src.unlink()

def sum_fasta_sizes(fasta_gz: Path) -> int:
    """Sum size= abundance annotations across all headers in a gzip FASTA."""
    total = 0
    with gzip.open(fasta_gz, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                m = re.search(r'size=(\d+)', line)
                if m:
                    total += int(m.group(1))
    return total

def derive_sample_name(reads1: str) -> str:
    name = Path(reads1).name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    stripped = re.sub(r'[_\-][Rr]1(?:[_\-].*)?$', '', name)
    safe = (stripped if stripped else name).replace("-", "_")
    return safe

def fmt_tsv(rows):
    col_w = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    return "\n".join(
        "  ".join(r[i].ljust(col_w[i]) for i in range(len(r))) for r in rows
    )

################################################################################
# 3. Define the main function
################################################################################

def main():
    opts = parse_args()
    reads1       = opts.reads1
    reads2       = opts.reads2
    output_dir   = Path(opts.output_dir)
    sample_name  = opts.sample_name
    nslots       = opts.nslots
    min_length   = opts.min_length
    minovlen     = opts.fastq_minovlen
    maxdiffs     = opts.fastq_maxdiffs
    maxee        = opts.fastq_maxee
    min_size     = opts.min_size
    abskew       = opts.abskew
    overwrite    = opts.overwrite == "t"

    ###########################################################################
    # Step 0: Validate inputs, create output directories, and check dependencies
    ###########################################################################

    # Validate input files
    if not os.path.isfile(reads1):
        log_error(f"R1 file does not exist: {reads1}")
        sys.exit(1)
    if not os.path.isfile(reads2):
        log_error(f"R2 file does not exist: {reads2}")
        sys.exit(1)

    # Check dependencies
    log("Checking dependencies...")
    for tool in ("vsearch",):
        if not shutil.which(tool):
            log_error(f"{tool} not found. Please install it or add it to PATH.")
            sys.exit(1)

    # Prepare per-sample output directory
    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    for sub in ("01-merged", "02-filtered", "03-derep", "04-chimera-checked", "logs", "stats"):
        (output_dir / sub).mkdir(parents=True, exist_ok=True)

    # Derive sample name
    if sample_name is None:
        sample_name = derive_sample_name(reads1)
    log(f"Processing sample: {sample_name}")

    ###########################################################################
    # Step 1: Merge paired-end reads
    ###########################################################################

    log("Step 1/4: Merging with vsearch...")

    merged_fastq = output_dir / "01-merged" / f"{sample_name}-01-merged.fastq"
    merged_fastq_gz = output_dir / "01-merged" / f"{sample_name}-01-merged.fastq.gz"

    merge_cmd = [
        "vsearch",
        "--fastq_mergepairs", reads1,
        "--reverse",          reads2,
        "--fastqout",         merged_fastq,
        "--fastq_minovlen",   minovlen,
        "--fastq_maxdiffs",   maxdiffs,
        "--threads",          nslots,
    ]
    merge_result = run(merge_cmd,
                       log_path=output_dir / "logs" / f"{sample_name}_merge.log",
                       check=False)

    if merge_result.returncode != 0:
        log_error(f"vsearch merge failed for {sample_name}.")
        sys.exit(1)

    compress(merged_fastq, merged_fastq_gz)

    m = re.search(r'([\d,]+)\s+Pairs', merge_result.stderr)
    pairs_merged_in = int(m.group(1).replace(',', '')) if m else 0
    m = re.search(r'([\d,]+)\s+Merged', merge_result.stderr)
    pairs_merged_out = int(m.group(1).replace(',', '')) if m else 0
    pct_merge_out = f"{pairs_merged_out / pairs_merged_in * 100:.2f}" if pairs_merged_in > 0 else "0.00"

    merge_rows = [
        ["sample", "pairs_in", "pairs_merged", "percent_merged"],
        [sample_name, str(pairs_merged_in), str(pairs_merged_out), pct_merge_out],
    ]

    ###########################################################################
    # Step 3: Filter by expected errors → FASTA, label reads with sample name
    ###########################################################################

    log("Step 2/4: Filtering by expected errors with vsearch...")
    filtered_fasta = output_dir / "02-filtered" / f"{sample_name}-02-filtered.fasta"
    filtered_fasta_gz = output_dir / "02-filtered" / f"{sample_name}-02-filtered.fasta.gz"

    filter_cmd = [
        "vsearch",
        "--fastq_filter", merged_fastq_gz,
        "--fastq_maxee",  maxee,
        "--fastq_minlen", min_length,
        "--fastaout",     filtered_fasta,
        "--relabel",      f"{sample_name}.",
        "--fasta_width",  0,
    ]
    filter_result = run(filter_cmd,
                        log_path=output_dir / "logs" / f"{sample_name}_filter.log",
                        check=False)

    if filter_result.returncode != 0:
        log_error(f"vsearch filter failed for {sample_name}.")
        sys.exit(1)

    reads_filter_in = pairs_merged_out
    # Count the number of reads that passed the filter
    reads_filter_out = 0
    with open(filtered_fasta, "r") as f:
        for line in f:
            if line.startswith(">"):
                reads_filter_out += 1
    pct_filter = f"{reads_filter_out / reads_filter_in * 100:.2f}" if reads_filter_in > 0 else "0.00"
    

    compress(filtered_fasta, filtered_fasta_gz)

    filter_rows = [
        ["sample", "reads_in", "reads_passed", "percent_passed"],
        [sample_name, str(reads_filter_in), str(reads_filter_out), pct_filter],
    ]

    log(f"  Filtered FASTA written to: {filtered_fasta_gz}")

    ###########################################################################
    # Step 4: Dereplication
    ###########################################################################

    log("Step 3/4: Dereplicating with vsearch...")

    derep_fasta = output_dir / "03-derep" / f"{sample_name}-03-derep.fasta"
    derep_fasta_gz = output_dir / "03-derep" / f"{sample_name}-03-derep.fasta.gz"
    derep_uc = output_dir / "03-derep" / f"{sample_name}-03-derep.uc"

    derep_cmd = [
        "vsearch",
        "--derep_fulllength", filtered_fasta_gz,
        "--output",           derep_fasta,
        "--uc",               derep_uc,
        "--minuniquesize",    min_size,
        "--sizeout",
        "--fasta_width",      0,
    ]
    derep_result = run(derep_cmd,
                       log_path=output_dir / "logs" / f"{sample_name}_derep.log",
                       check=False)

    if derep_result.returncode != 0:
        log_error(f"vsearch dereplication failed for {sample_name}.")
        sys.exit(1)

    compress(derep_fasta, derep_fasta_gz)

    m = re.search(r'(\d+) unique sequences', derep_result.stderr)
    n_unique = int(m.group(1)) if m else 0
    pct_unique = f"{n_unique / reads_filter_out * 100:.2f}" if reads_filter_out > 0 else "0.00"

    log(f"  Unique sequences after dereplication: {n_unique}")

    derep_rows = [
        ["sample", "seqs_filtered", "seqs_unique", "percent_unique"],
        [sample_name, str(reads_filter_out), str(n_unique), pct_unique],
    ]

    ###########################################################################
    # Step 4: Chimera detection (de novo)
    ###########################################################################

    log("Step 4/4: Chimera detection with vsearch...")

    chimera_checked_fasta = output_dir / "04-chimera-checked" / f"{sample_name}-04-chimera-checked.fasta"
    chimera_checked_fasta_gz = output_dir / "04-chimera-checked" / f"{sample_name}-04-chimera-checked.fasta.gz"
    chimera_fasta = output_dir / "04-chimera-checked" / f"{sample_name}-04-chimeras.fasta"
    chimera_fasta_gz = output_dir / "04-chimera-checked" / f"{sample_name}-04-chimeras.fasta.gz"

    chimera_cmd = [
        "vsearch",
        "--uchime_denovo",  derep_fasta_gz,
        "--sizein",
        "--sizeout",
        "--abskew",         abskew,
        "--nonchimeras",    chimera_checked_fasta,
        "--chimeras",       chimera_fasta,
        "--fasta_width",    0,
    ]
    chimera_result = run(chimera_cmd,
                         log_path=output_dir / "logs" / f"{sample_name}_chimera.log",
                         check=False)

    if chimera_result.returncode != 0:
        log_error(f"vsearch chimera detection failed for {sample_name}.")
        sys.exit(1)

    compress(chimera_checked_fasta, chimera_checked_fasta_gz)
    compress(chimera_fasta, chimera_fasta_gz)

    m = re.search(r'Found ([\d,]+)', chimera_result.stderr)
    n_chimeras = int(m.group(1).replace(',', '')) if m else 0
    m = re.search(r'([\d,]+).*?non-chimeras', chimera_result.stderr)
    n_nonchimeras = int(m.group(1).replace(',', '')) if m else 0
    pct_chimera = f"{n_chimeras / n_unique * 100:.2f}" if n_unique > 0 else "0.00"

    abund_derep = sum_fasta_sizes(derep_fasta_gz)
    abund_chimera_checked = sum_fasta_sizes(chimera_checked_fasta_gz)
    pct_abund_retained = f"{abund_chimera_checked / abund_derep * 100:.2f}" if abund_derep > 0 else "0.00"

    log(f"  Chimeras: {n_chimeras}  Non-chimeras: {n_nonchimeras}")
    log(f"  Abundance retained after chimera removal: {abund_chimera_checked}/{abund_derep} ({pct_abund_retained}%)")

    chimera_rows = [
        ["sample", "seqs_in", "chimeras", "nonchimeras", "pct_chimeric_seqs",
         "abund_in", "abund_nonchimeric", "pct_abund_retained"],
        [sample_name, str(n_unique), str(n_chimeras), str(n_nonchimeras), pct_chimera,
         str(abund_derep), str(abund_chimera_checked), pct_abund_retained],
    ]

    if abund_derep != reads_filter_out:
        log_warn(f"Abundance mismatch: {abund_derep} != {reads_filter_out}")

    # Write a single consolidated stats table for the sample
    stats_out = output_dir / "stats" / f"{sample_name}_stats.tsv"
    stats_header = [
        "sample", "pairs_in", "pairs_merged", "percent_merged",
        "reads_passed", "percent_passed", "seqs_unique", "percent_unique",
        "chimeras", "nonchimeras", "pct_chimeric_seqs",
        "abund_in", "abund_nonchimeric", "pct_abund_retained",
    ]
    stats_row = [
        sample_name, str(pairs_merged_in), str(pairs_merged_out), pct_merge_out,
        str(reads_filter_out), pct_filter, str(n_unique), pct_unique,
        str(n_chimeras), str(n_nonchimeras), pct_chimera,
        str(abund_derep), str(abund_chimera_checked), pct_abund_retained,
    ]
    stats_out.write_text("\t".join(stats_header) + "\n" + "\t".join(stats_row) + "\n")

    ###########################################################################
    # Step 5: Write summary report
    ###########################################################################

    log("Generating summary report...")
    report_out = output_dir / f"{sample_name}_summary_report.txt"

    report = (
        f"{'=' * 80}\n"
        f"VSEARCH Per-Sample Processing Report\n"
        f"{'=' * 80}\n"
        f"Date: {datetime.now()}\n"
        f"Sample: {sample_name}\n"
        f"R1: {reads1}\n"
        f"R2: {reads2}\n"
        f"Output directory: {output_dir}\n"
        f"\n"
        f"Parameters:\n"
        f"-----------\n"
        f"  Threads:             {nslots}\n"
        f"  Min read length:     {min_length} bp\n"
        f"  Min merge overlap:   {minovlen} bp\n"
        f"  Max merge diffs:     {maxdiffs}\n"
        f"  Max expected errors: {maxee}\n"
        f"  Min unique size:     {min_size}\n"
        f"  Chimera abskew:      {abskew}\n"
        f"\n"
        f"Output files:\n"
        f"-------------\n"
        f"  Merged reads:     {output_dir}/01-merged/\n"
        f"  Filtered FASTA:   {output_dir}/02-filtered/\n"
        f"  Dereplicated:     {output_dir}/03-derep/\n"
        f"  Chimera filter:   {output_dir}/04-chimera/\n"
        f"  Logs:             {output_dir}/logs/\n"
        f"  Stats:            {output_dir}/stats/\n"
        f"\n"
        f"{'=' * 80}\n"
        f"\n"
        f"Merge statistics:\n"
        f"\n"
        f"{fmt_tsv(merge_rows)}\n"
        f"\n"
        f"Filter statistics:\n"
        f"\n"
        f"{fmt_tsv(filter_rows)}\n"
        f"\n"
        f"Dereplication statistics:\n"
        f"\n"
        f"{fmt_tsv(derep_rows)}\n"
        f"\n"
        f"Chimera statistics:\n"
        f"\n"
        f"{fmt_tsv(chimera_rows)}\n"
    )

    report_out.write_text(report)
    # print(report)
    log("\033[0;32m2.2.1-vsearch-pipeline.py completed successfully\033[0m")

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
