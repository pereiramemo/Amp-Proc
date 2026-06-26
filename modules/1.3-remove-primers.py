#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# DEV ONLY — comment out before production use
""" reads1      = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
reads2      = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
output_dir  = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/03_primer_removal/sample1"
trimmed_dir = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/03_primer_removal/trimmed"
primer_fwd  = "GTGYCAGCMGCCGCGGTAA"
primer_rev  = "CCGYCAATTYMTTTRAGTTT"
nslots      = 4
overwrite   = True
compress    = True
error_rate = 0.1
min_overlap = 5
min_length = 50
 """

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Remove primers from a single paired-end sample using cutadapt"
    )
    p.add_argument("--reads1",            required=True,              help="R1 FASTQ file (required)")
    p.add_argument("--reads2",            required=True,              help="R2 FASTQ file (required)")
    p.add_argument("-o", "--output_dir",  required=True,              help="Per-sample output directory for logs and stats (required)")
    p.add_argument("--trimmed_dir",       default=None,               help="Directory for trimmed FASTQ output [default: {output_dir}/trimmed]")
    p.add_argument("--primer_fwd",        required=True,              help="Forward primer sequence (5' to 3') (required)")
    p.add_argument("--primer_rev",        required=True,              help="Reverse primer sequence (5' to 3') (required)")
    p.add_argument("--nslots",            type=int,   default=12,     help="Number of threads [default=12]")
    p.add_argument("--error_rate",        type=float, default=0.1,    help="Maximum allowed error rate [default=0.1]")
    p.add_argument("--min_overlap",       type=int,   default=3,      help="Minimum primer-read overlap [default=3]")
    p.add_argument("--min_length",        type=int,   default=50,     help="Discard reads shorter than this after trimming [default=50]")
    p.add_argument("--discard_untrimmed", choices=["t", "f"], default="f", help="Discard reads where no primer was found [default=f]")
    p.add_argument("--compress",          choices=["t", "f"], default="t", help="Compress output files with gzip [default=t]")
    p.add_argument("--overwrite",         choices=["t", "f"], default="f", help="Overwrite existing per-sample output directory [default=f]")
    return p.parse_args()


def log(msg):
    print(f"[INFO] {msg}")

def log_warn(msg):
    print(f"\033[1;33m[WARN]\033[0m {msg}", file=sys.stderr)

def log_error(msg):
    print(f"\033[0;31m[ERROR]\033[0m {msg}", file=sys.stderr)


def derive_sample_name(reads1: str) -> str:
    name = Path(reads1).name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    # strip trailing _R1 / -R1 and everything after (e.g. _R1_001_redu → stripped)
    stripped = re.sub(r'[_\-][Rr]1(?:[_\-].*)?$', '', name)
    return stripped if stripped else name

################################################################################
# 3. Define the main function
################################################################################

def main():
    opts = parse_args()
    reads1            = opts.reads1
    reads2            = opts.reads2
    output_dir        = Path(opts.output_dir)
    primer_fwd        = opts.primer_fwd
    primer_rev        = opts.primer_rev
    nslots            = opts.nslots
    error_rate        = opts.error_rate
    min_overlap       = opts.min_overlap
    min_length        = opts.min_length
    discard_untrimmed = opts.discard_untrimmed == "t"
    compress          = opts.compress == "t"
    overwrite         = opts.overwrite == "t"

    # Validate input files
    if not os.path.isfile(reads1):
        log_error(f"R1 file does not exist: {reads1}")
        sys.exit(1)
    if not os.path.isfile(reads2):
        log_error(f"R2 file does not exist: {reads2}")
        sys.exit(1)

    # Check cutadapt is available
    log("Checking dependencies...")
    if not shutil.which("cutadapt"):
        log_error("cutadapt not found. Please install it or add it to PATH.")
        sys.exit(1)

    # Prepare per-sample output directory (logs, stats)
    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    logs_dir  = output_dir / "logs"
    stats_dir = output_dir / "stats"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # Shared trimmed output directory (not deleted on overwrite)
    trimmed_dir = Path(opts.trimmed_dir) if opts.trimmed_dir else output_dir / "trimmed"
    trimmed_dir.mkdir(parents=True, exist_ok=True)

    # Derive sample name
    sample_name = derive_sample_name(reads1)
    log(f"Processing sample: {sample_name}")
    log(f"Forward primer: {primer_fwd}")
    log(f"Reverse primer: {primer_rev}")

    # Define output file paths
    ext    = ".fastq.gz" if compress else ".fastq"
    r1_out = trimmed_dir / f"{sample_name}_R1_trimmed{ext}"
    r2_out = trimmed_dir / f"{sample_name}_R2_trimmed{ext}"
    log_out = logs_dir   / f"{sample_name}_cutadapt.log"

    # Build cutadapt command
    cmd = [
        "cutadapt",
        "-g", primer_fwd,
        "-G", primer_rev,
        "-o", str(r1_out),
        "-p", str(r2_out),
        "-j", str(nslots),
        "-e", str(error_rate),
        "-O", str(min_overlap),
        "-m", str(min_length),
    ]
    if discard_untrimmed:
        cmd.append("--discard-untrimmed")
    cmd += [reads1, reads2]

    # Run cutadapt
    log(f"Running cutadapt...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, 
                            stderr=subprocess.STDOUT, text=True)
    log_out.write_text(result.stdout)

    if result.returncode != 0:
        log_error(f"cutadapt failed for sample {sample_name}")
        log_error(f"Check log file: {log_out}")
        sys.exit(1)

    # Extract summary statistics from log
    log_text = result.stdout
    m = re.search(r'Total read pairs processed:\s+([\d,]+)', log_text)
    total_pairs = int(m.group(1).replace(',', '')) if m else 0
    m = re.search(r'Pairs written \(passing filters\):\s+([\d,]+)', log_text)
    pairs_written = int(m.group(1).replace(',', '')) if m else 0

    percent_trimmed = f"{pairs_written / total_pairs * 100:.2f}" if total_pairs > 0 else "0.00"

    summary_file = stats_dir / "summary.tsv"
    summary_file.write_text(
        "sample\ttotal_pairs\ttrimmed_pairs\tpercent_trimmed\n"
        f"{sample_name}\t{total_pairs}\t{pairs_written}\t{percent_trimmed}\n"
    )

    # Generate human-readable summary report
    log("Generating summary report...")
    summary_txt  = output_dir / "summary_report.txt"
    discard_label = "yes" if discard_untrimmed else "no"
    compress_label = "yes" if compress else "no"

    tsv_lines = summary_file.read_text().splitlines()
    rows  = [line.split("\t") for line in tsv_lines]
    col_w = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    tsv_formatted = "\n".join(
        "  ".join(r[i].ljust(col_w[i]) for i in range(len(r)))
        for r in rows
    )

    report = (
        f"{'=' * 80}\n"
        f"Primer Removal Report - cutadapt\n"
        f"{'=' * 80}\n"
        f"Date: {datetime.now()}\n"
        f"Sample: {sample_name}\n"
        f"R1: {reads1}\n"
        f"R2: {reads2}\n"
        f"Output directory: {output_dir}\n"
        f"\n"
        f"Primers:\n"
        f"--------\n"
        f"Forward primer (5'-3'): {primer_fwd}\n"
        f"Reverse primer (5'-3'): {primer_rev}\n"
        f"\n"
        f"Parameters:\n"
        f"-----------\n"
        f"Threads: {nslots}\n"
        f"Error rate: {error_rate}\n"
        f"Minimum overlap: {min_overlap}\n"
        f"Minimum length: {min_length}\n"
        f"Discard untrimmed: {discard_label}\n"
        f"Compress output: {compress_label}\n"
        f"\n"
        f"Output files:\n"
        f"-------------\n"
        f"- Trimmed R1: {r1_out}\n"
        f"- Trimmed R2: {r2_out}\n"
        f"- Processing log: {log_out}\n"
        f"- Summary statistics: {summary_file}\n"
        f"- This report: {summary_txt}\n"
        f"\n"
        f"{'=' * 80}\n"
        f"\n"
        f"Summary Statistics:\n"
        f"\n"
        f"{tsv_formatted}\n"
    )

    summary_txt.write_text(report)
    print(report)

    log("\033[0;32m1.3-primer_removal_cutadapt.py completed successfully\033[0m")
    log(f"Trimmed reads available in: {trimmed_dir}")
    log(f"Summary statistics: {summary_file}")

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
