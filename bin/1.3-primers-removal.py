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
from pathlib import Path

# Import shared helpers from bin/toolbox.py (sibling module).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolbox import log, log_warn, log_error, derive_sample_name, build_log

# DEV ONLY — comment out before production use
""" 
reads1      = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
reads2      = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
output_dir  = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/03_primer_removal/sample1"
primer_fwd  = "GTGYCAGCMGCCGCGGTAA"
primer_rev  = "CCGYCAATTYMTTTRAGTTT"
nslots      = 4
overwrite   = True
compress    = True
error_rate = 0.1
min_overlap = 5
min_length = 50
 """

SCRIPT_NAME = "1.3-primers-removal.py"
SCRIPT_DESC = "Remove primers from a single paired-end sample using cutadapt."

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Remove primers from a single paired-end sample using cutadapt"
    )
    p.add_argument("--reads1",            required=True,              help="R1 FASTQ file (required)")
    p.add_argument("--reads2",            required=True,              help="R2 FASTQ file (required)")
    p.add_argument("--output_dir",  required=True,              help="Per-sample output directory for logs and stats (required)")
    p.add_argument("--sample_name",       default=None,               help="Sample name [default: derived from R1 filename]")
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

################################################################################
# 3. Define the main function
################################################################################

def main():

    ###########################################################################
    # Step 1: Define input variables
    ###########################################################################

    opts = parse_args()
    reads1            = opts.reads1
    reads2            = opts.reads2
    output_dir        = Path(opts.output_dir)
    sample_name       = opts.sample_name
    primer_fwd        = opts.primer_fwd
    primer_rev        = opts.primer_rev
    nslots            = opts.nslots
    error_rate        = opts.error_rate
    min_overlap       = opts.min_overlap
    min_length        = opts.min_length
    discard_untrimmed = opts.discard_untrimmed == "t"
    compress          = opts.compress == "t"
    overwrite         = opts.overwrite == "t"

    ###########################################################################
    # Step 2: Validate inputs and create output directories
    ###########################################################################

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

    # Prepare per-sample output directory (output, logs, stats)
    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    results_dir = output_dir / "output"
    logs_dir    = output_dir / "logs"
    stats_dir   = output_dir / "stats"
    results_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    # Derive sample name from R1 filename if not provided
    if sample_name is None:
        sample_name = derive_sample_name(reads1, strip_read_suffix=True)
    log(f"Processing sample: {sample_name}")
    log(f"Forward primer: {primer_fwd}")
    log(f"Reverse primer: {primer_rev}")

    # Define output file paths
    ext       = ".fastq.gz" if compress else ".fastq"
    r1_out    = results_dir / f"{sample_name}_R1_trimmed{ext}"
    r2_out    = results_dir / f"{sample_name}_R2_trimmed{ext}"
    log_out   = logs_dir    / f"1.3-primers-removal-{sample_name}.log"
    stats_out = stats_dir   / f"1.3-primers-removal-{sample_name}-stats.tsv"

    ###########################################################################
    # Step 3: Run cutadapt to remove primers
    ###########################################################################

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

    discard_label  = "yes" if discard_untrimmed else "no"
    compress_label = "yes" if compress else "no"
    params = [
        f"Threads: {nslots}",
        f"Error rate: {error_rate}",
        f"Minimum overlap: {min_overlap}",
        f"Minimum length: {min_length}",
        f"Discard untrimmed: {discard_label}",
        f"Compress output: {compress_label}",
    ]
    inputs = [f"R1: {reads1}", f"R2: {reads2}",
              f"Forward primer: {primer_fwd}", f"Reverse primer: {primer_rev}"]

    # Run cutadapt
    log("Running cutadapt...")
    result = subprocess.run(cmd, stdout=subprocess.PIPE,
                            stderr=subprocess.STDOUT, text=True)

    if result.returncode != 0:
        log_error(f"cutadapt failed for sample {sample_name}")
        log_error(f"Check log file: {log_out}")
        log_out.write_text(build_log(
            SCRIPT_NAME, SCRIPT_DESC, sample_name,
            inputs=inputs, params=params,
            outputs=[f"Trimmed R1: {r1_out}", f"Trimmed R2: {r2_out}"],
            command=" ".join([SCRIPT_NAME] + sys.argv[1:]),
            exit_status=result.returncode,
            tool_log=result.stdout,
        ))
        sys.exit(1)

    ###########################################################################
    # Step 5: Create log and stats files
    ###########################################################################

    log_text = result.stdout
    m = re.search(r'Total read pairs processed:\s+([\d,]+)', log_text)
    total_pairs = int(m.group(1).replace(',', '')) if m else 0
    m = re.search(r'Pairs written \(passing filters\):\s+([\d,]+)', log_text)
    pairs_written = int(m.group(1).replace(',', '')) if m else 0

    percent_trimmed = f"{pairs_written / total_pairs * 100:.2f}" if total_pairs > 0 else "0.00"

    # Stats table: samples as rows, statistics as columns
    stats_out.write_text(
        "sample\ttotal_pairs\ttrimmed_pairs\tpercent_trimmed\n"
        f"{sample_name}\t{total_pairs}\t{pairs_written}\t{percent_trimmed}\n"
    )

    log("\033[0;32m1.3-primers-removal.py completed successfully\033[0m")
    
    log_out.write_text(build_log(
        SCRIPT_NAME, SCRIPT_DESC, sample_name,
        inputs=inputs, params=params,
        outputs=[
            f"Trimmed R1: {r1_out}",
            f"Trimmed R2: {r2_out}",
            f"Statistics: {stats_out}",
        ],
        command=" ".join([SCRIPT_NAME] + sys.argv[1:]),
        exit_status=0,
        tool_log=result.stdout,
    ))

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
