#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# DEV ONLY — comment out before production use
""" reads1 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
reads2 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
output_dir = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/01_fastp/sample1"
nslots = 4
disable_adapter_trimming = True
overwrite = True
 """
################################################################################
# 2. Define Functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Run fastp QC report on a single paired-end sample (no filtering applied)"
    )
    p.add_argument("--reads1",                    required=True,              help="R1 FASTQ file (required)")
    p.add_argument("--reads2",                    required=True,              help="R2 FASTQ file (required)")
    p.add_argument("-o", "--output_dir",          required=True,              help="Output directory (required)")
    p.add_argument("--nslots",                    type=int,   default=12,     help="Number of threads [default=12]")
    p.add_argument("--min_length",                type=int,   default=50,     help="Minimum read length (reporting only) [default=50]")
    p.add_argument("--qualified_quality_phred",   type=int,   default=20,     help="Phred score for qualified base (reporting only) [default=20]")
    p.add_argument("--unqualified_percent_limit", type=int,   default=40,     help="Max percent of unqualified bases (reporting only) [default=40]")
    p.add_argument("--disable_adapter_trimming",  choices=["t", "f"], default="t", help="Disable adapter trimming [default=t]")
    p.add_argument("--html_report",               choices=["t", "f"], default="t", help="Generate HTML report [default=t]")
    p.add_argument("--json_report",               choices=["t", "f"], default="t", help="Keep JSON report after stats extraction [default=t]")
    p.add_argument("--overwrite",                 choices=["t", "f"], default="f", help="Overwrite existing output directory [default=f]")
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
            return name[: -len(ext)]
    return name

def main():
    opts = parse_args()
    reads1 = opts.reads1
    reads2 = opts.reads2
    output_dir = Path(opts.output_dir)
    nslots = opts.nslots
    disable_adapter_trimming = opts.disable_adapter_trimming == "t"
    html_report = opts.html_report == "t"
    json_report = opts.json_report == "t"
    overwrite = opts.overwrite == "t"

    # Validate input files
    if not os.path.isfile(reads1):
        log_error(f"R1 file does not exist: {reads1}")
        sys.exit(1)
    if not os.path.isfile(reads2):
        log_error(f"R2 file does not exist: {reads2}")
        sys.exit(1)

    # Check fastp is available
    log("Checking dependencies...")
    if not shutil.which("fastp"):
        log_error("fastp not found. Please install it or add it to PATH.")
        sys.exit(1)

    # Prepare output directory
    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    reports_dir = output_dir / "reports"
    stats_dir   = output_dir / "stats"
    reports_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # Derive sample name
    sample_name = derive_sample_name(reads1)
    log(f"Processing sample: {sample_name}")

    # Paths for fastp outputs
    html_out = reports_dir / f"{sample_name}_fastp.html"
    json_out = reports_dir / f"{sample_name}_fastp.json"
    log_out  = reports_dir / f"{sample_name}_fastp.log"

    # Build fastp command
    cmd = [
        "fastp",
        "-i", reads1,
        "-I", reads2,
        "-w", str(nslots),
        "--disable_quality_filtering",
        "--disable_length_filtering",
        "--disable_trim_poly_g",
        "--json", str(json_out),
    ]

    if html_report:
        cmd += ["--html", str(html_out)]
    else:
        cmd += ["--html", "/dev/null"]

    if disable_adapter_trimming:
        cmd.append("--disable_adapter_trimming")

    # Run fastp
    log("Running fastp...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, 
                            stderr=subprocess.STDOUT, text=True)
    print(result.stdout, end="")
    with open(log_out, "w") as fh:
        fh.write(result.stdout)

    if result.returncode != 0:
        log_error(f"fastp failed for sample {sample_name}")
        sys.exit(1)

    # Extract summary statistics from JSON
    with open(json_out) as f:
        bf = json.load(f)["summary"]["before_filtering"]

    total_reads    = bf["total_reads"]
    total_bases    = bf["total_bases"]
    q20_bases      = bf["q20_bases"]
    q20_rate       = bf["q20_rate"]
    q30_bases      = bf["q30_bases"]
    q30_rate       = bf["q30_rate"]
    r1_mean_length = bf["read1_mean_length"]
    r2_mean_length = bf["read2_mean_length"]
    gc_content     = bf["gc_content"]

    # Write summary TSV
    summary_file = stats_dir / "summary.tsv"
    header = "\t".join([
        "sample", "total_reads", "total_bases",
        "q20_bases", "q20_rate", "q30_bases", "q30_rate",
        "read1_mean_length", "read2_mean_length", "gc_content",
    ])
    row = "\t".join(str(v) for v in [
        sample_name, total_reads, total_bases,
        q20_bases, q20_rate, q30_bases, q30_rate,
        r1_mean_length, r2_mean_length, gc_content,
    ])
    

    if not json_report:
        json_out.unlink(missing_ok=True)

    # Generate human-readable summary report
    log("Generating summary report...")
    summary_txt = output_dir / "summary_report.txt"
    adapter_label = "disabled" if disable_adapter_trimming else "enabled"

    tsv_rows = [header.split("\t"), row.split("\t")]
    col_w = [max(len(r[i]) for r in tsv_rows) for i in range(len(tsv_rows[0]))]
    tsv_formatted = "\n".join(
        "  ".join(r[i].ljust(col_w[i]) for i in range(len(r)))
        for r in tsv_rows
    )

    report = (
        f"{'=' * 80}\n"
        f"Quality Check Report - fastp\n"
        f"{'=' * 80}\n"
        f"Date: {datetime.now()}\n"
        f"Sample: {sample_name}\n"
        f"R1: {reads1}\n"
        f"R2: {reads2}\n"
        f"Output directory: {output_dir}\n"
        f"\n"
        f"Parameters:\n"
        f"-----------\n"
        f"Threads: {nslots}\n"
        f"Mode: Report only (no filtering applied)\n"
        f"Adapter trimming: {adapter_label}\n"
        f"\n"
        f"Output files:\n"
        f"-------------\n"
        f"- HTML report: {html_out}\n"
        f"- JSON report: {json_out}\n"
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
    # print(report) 
    log("\033[0;32m1.1-quality-check.py completed successfully\033[0m")
    
################################################################################
# 3. Execute main function
################################################################################

if __name__ == "__main__":
    main()
