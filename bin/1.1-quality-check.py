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
from pathlib import Path

# Import shared helpers from bin/toolbox.py (sibling module).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolbox import log, log_warn, log_error, derive_sample_name, build_log

# DEV ONLY — comment out before production use
"""
reads1 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
reads2 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
output_dir = Path("/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/01_fastp/sample1")
nslots = 4
overwrite = True
 """

SCRIPT_NAME = "1.1-quality-check.py"
SCRIPT_DESC = "Run a fastp QC report on a single paired-end sample (no filtering applied)."

################################################################################
# 2. Define Functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Run fastp QC report on a single paired-end sample (no filtering applied)"
    )
    p.add_argument("--reads1",                    required=True,              help="R1 FASTQ file (required)")
    p.add_argument("--reads2",                    required=True,              help="R2 FASTQ file (required)")
    p.add_argument("--output_dir",                required=True,              help="Output directory (required)")
    p.add_argument("--sample_name",               default=None,   help="Sample name [default: derived from R1 filename]")
    p.add_argument("--nslots",                    type=int,   default=12,     help="Number of threads [default=12]")
    p.add_argument("--html_report",               choices=["t", "f"], default="t", help="Generate HTML report [default=t]")
    p.add_argument("--json_report",               choices=["t", "f"], default="t", help="Keep JSON report after stats extraction [default=t]")
    p.add_argument("--overwrite",                 choices=["t", "f"], default="f", help="Overwrite existing output directory [default=f]")
    return p.parse_args()

################################################################################
# 3. Define main function
################################################################################

def main():

    ###########################################################################
    # Step 1: Define input variables
    ###########################################################################

    opts = parse_args()
    reads1 = opts.reads1
    reads2 = opts.reads2
    output_dir = Path(opts.output_dir)
    sample_name = opts.sample_name
    nslots = opts.nslots
    html_report = opts.html_report == "t"
    json_report = opts.json_report == "t"
    overwrite = opts.overwrite == "t"

    ###########################################################################
    # Step 2: Validate input files and tools, and define and create output directories
    ###########################################################################

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

    results_dir = output_dir / "output"
    stats_dir   = output_dir / "stats"
    logs_dir    = output_dir / "logs"
    results_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Derive sample name
    if sample_name is None:
        sample_name = derive_sample_name(reads1)
    log(f"Processing sample: {sample_name}")

    # Paths for fastp outputs
    html_out  = results_dir / f"{sample_name}_fastp.html"
    json_out  = results_dir / f"{sample_name}_fastp.json"
    log_out   = logs_dir    / f"1.1-quality-check-{sample_name}.log"
    stats_out = stats_dir   / f"1.1-quality-check-{sample_name}-stats.tsv"

    ###########################################################################
    # Step 3: Build and execute fastp command
    ###########################################################################

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

    params = [f"Threads: {nslots}", "Mode: report only (no filtering applied)"]

    # Run fastp
    log("Running fastp...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)
    print(result.stdout, end="")

    ###########################################################################
    # Step 4: Check fastp result
    ###########################################################################

    if result.returncode != 0:
        log_error(f"fastp failed for sample {sample_name}")
        log_out.write_text(build_log(
            SCRIPT_NAME, SCRIPT_DESC, sample_name,
            inputs=[f"R1: {reads1}", f"R2: {reads2}"],
            params=params,
            outputs=[f"HTML report: {html_out}", f"JSON report: {json_out}"],
            command=" ".join(str(c) for c in cmd),
            exit_status=result.returncode,
            tool_log=result.stdout,
        ))
        sys.exit(1)

    ###########################################################################
    # Step 5: Create log and stats files
    ###########################################################################

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

    # Write summary TSV (samples as rows, statistics as columns)
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
    stats_out.write_text(header + "\n" + row + "\n")

    if not json_report:
        json_out.unlink(missing_ok=True)

    log("\033[0;32m1.1-quality-check.py completed successfully\033[0m")

    # Write the standardized log file (general info + fastp log)
    log_out.write_text(build_log(
        SCRIPT_NAME, SCRIPT_DESC, sample_name,
        inputs=[f"R1: {reads1}", f"R2: {reads2}"],
        params=params,
        outputs=[
            f"HTML report: {html_out}",
            f"JSON report: {json_out if json_report else 'not kept'}",
            f"Statistics: {stats_out}",
        ],
        command=" ".join(str(c) for c in cmd),
        exit_status=0,
        tool_log=result.stdout,
    ))

################################################################################
# 4. Execute main function
################################################################################

if __name__ == "__main__":
    main()
