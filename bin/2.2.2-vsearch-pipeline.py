#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import gzip
import os
import shutil
import subprocess
import sys
from pathlib import Path

# Import shared helpers from bin/toolbox.py (sibling module).
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from toolbox import log, log_warn, log_error, build_log

# DEV ONLY — comment out before production use
"""
samples_dir = Path("/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch")
output_dir  = Path("/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/06_otu")
nslots      = 4
identity    = 0.97
overwrite   = True
"""

SCRIPT_NAME = "2.2.2-vsearch-pipeline.py"
SCRIPT_DESC = ("OTU construction: pool per-sample chimera-checked FASTAs -> "
               "cluster -> OTU table.")

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="OTU construction: pool per-sample chimera-checked FASTAs → cluster → OTU table"
    )
    p.add_argument("--samples_dir",      required=True,              help="Parent directory containing per-sample 2.2.1 output subdirectories (required)")
    p.add_argument("--output_dir",       required=True,              help="Directory to write all output (required)")
    p.add_argument("--nslots",           type=int,   default=12,     help="Number of threads [default=12]")
    p.add_argument("--identity",         type=float, default=0.97,   help="OTU clustering identity threshold [default=0.97]")
    p.add_argument("--overwrite",        choices=["t", "f"], default="f", help="Overwrite existing output directory [default=f]")
    return p.parse_args()

def run(cmd, check=True):
    result = subprocess.run(
        [str(c) for c in cmd],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True)
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

def count_fasta_seqs(fasta_path: Path) -> int:
    return sum(1 for line in open(fasta_path) if line.startswith(">"))

################################################################################
# 3. Define the main function
################################################################################

def main():

    ###########################################################################
    # Step 1: Define input variables
    ###########################################################################

    opts = parse_args()
    samples_dir = Path(opts.samples_dir)
    output_dir  = Path(opts.output_dir)
    nslots      = opts.nslots
    identity    = opts.identity
    overwrite   = opts.overwrite == "t"

    ###########################################################################
    # Step 2: Validate inputs, create output directories, and check dependencies
    ###########################################################################

    if not samples_dir.is_dir():
        log_error(f"Samples directory does not exist: {samples_dir}")
        sys.exit(1)

    log("Checking dependencies...")
    if not shutil.which("vsearch"):
        log_error("vsearch not found. Please install it or add it to PATH.")
        sys.exit(1)

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
    (results_dir / "otus").mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    stats_dir.mkdir(parents=True, exist_ok=True)

    ###########################################################################
    # Step 3: Concatenate all chimera-checked FASTA files
    ###########################################################################

    log("Collecting chimera-checked FASTA files...")
    fasta_files = sorted(samples_dir.glob("*/output/04-chimera-checked/*-04-chimera-checked.fasta.gz"))
    all_fasta = results_dir / "all_samples.fasta"
    all_fasta_gz = results_dir / "all_samples.fasta.gz"

    if not fasta_files:
        log_error(
            f"No *-04-chimera-checked.fasta.gz files found under: "
            f"{samples_dir}/*/output/04-chimera-checked/"
        )
        sys.exit(1)

    log(f"  Found {len(fasta_files)} sample file(s):")
    for f in fasta_files:
        log(f"    {f}")

    # Concatenate all samples into a single FASTA (decompress on the fly)
    with open(all_fasta, "w") as out_fh:
        for fasta_gz in fasta_files:
            with gzip.open(fasta_gz, "rt") as f_in:
                shutil.copyfileobj(f_in, out_fh)

    n_total = count_fasta_seqs(all_fasta)

    compress(all_fasta, all_fasta_gz)

    log(f"  Total sequences pooled: {n_total}")

    ###########################################################################
    # Step 4: OTU clustering
    ###########################################################################

    log(f"Step 1/1: Clustering OTUs at {identity} identity...")

    otus_fasta    = results_dir / "otus" / "otus.fasta"
    otus_fasta_gz = results_dir / "otus" / "otus.fasta.gz"
    otu_table     = results_dir / "otus" / "otu_table.tsv"

    cluster_cmd = [
        "vsearch",
        "--cluster_size", all_fasta_gz,
        "--id",           identity,
        "--centroids",    otus_fasta,
        "--otutabout",    otu_table,
        "--sizein",
        "--sizeout",
        "--relabel",      "OTU_",
        "--fasta_width",  0,
        "--threads",      nslots,
    ]
    cluster_result = run(cluster_cmd)

    n_otus = count_fasta_seqs(otus_fasta)

    compress(otus_fasta, otus_fasta_gz)

    log(f"  OTUs generated: {n_otus}")
    log(f"  OTU table written: {otu_table}")

    ###########################################################################
    # Step 5: Create log and stats files
    ###########################################################################

    stats_out = stats_dir / "2.2.2-vsearch-pipeline-stats.tsv"
    stats_out.write_text(
        "sample\tpooled_seqs\totus\n"
        f"all_samples\t{n_total}\t{n_otus}\n"
    )

    log("\033[0;32m2.2.2-vsearch-pipeline.py completed successfully\033[0m")
    tool_log = (
        f"### OTU clustering\n# {' '.join(str(c) for c in cluster_cmd)}\n"
        f"{cluster_result.stderr.rstrip()}\n"
    )
    log_out = logs_dir / "2.2.2-vsearch-pipeline.log"
    log_out.write_text(build_log(
        SCRIPT_NAME, SCRIPT_DESC, "all_samples",
        inputs=[f"Samples directory: {samples_dir}",
                f"Chimera-checked sample files: {len(fasta_files)}"],
        params=[f"Threads: {nslots}", f"OTU identity: {identity}"],
        outputs=[
            f"Concatenated FASTA: {all_fasta_gz}",
            f"OTU representatives: {otus_fasta_gz}",
            f"OTU table: {otu_table}",
            f"Statistics: {stats_out}",
        ],
        command=" ".join([SCRIPT_NAME] + sys.argv[1:]),
        exit_status=0,
        tool_log=tool_log,
    ))

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
