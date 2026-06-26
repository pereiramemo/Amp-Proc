#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import gzip
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# DEV ONLY — comment out before production use
"""
samples_dir = Path("/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch")
output_dir  = Path("/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/06_otu")
nslots      = 4
identity    = 0.97
overwrite   = True
"""

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
        text=True)
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

def fmt_tsv(rows):
    col_w = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    return "\n".join(
        "  ".join(r[i].ljust(col_w[i]) for i in range(len(r))) for r in rows
    )

def count_fasta_seqs(fasta_path: Path) -> int:
    return sum(1 for line in open(fasta_path) if line.startswith(">"))

################################################################################
# 3. Define the main function
################################################################################

def main():
    opts = parse_args()
    samples_dir = Path(opts.samples_dir)
    output_dir  = Path(opts.output_dir)
    nslots      = opts.nslots
    identity    = opts.identity
    overwrite   = opts.overwrite == "t"

    ###########################################################################
    # Step 0: Validate inputs, create output directories, and check dependencies
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

    for sub in ("otus", "logs", "stats"):
        (output_dir / sub).mkdir(parents=True, exist_ok=True)

    ###########################################################################
    # Step 1: Concatenate all chimera-checked FASTA files
    ###########################################################################

    log("Collecting chimera-checked FASTA files...")
    fasta_files = sorted(samples_dir.glob("*/04-chimera-checked/*-04-chimera-checked.fasta.gz"))
    all_fasta = output_dir / "all_samples.fasta"
    all_fasta_gz = output_dir / "all_samples.fasta.gz"
  
    if not fasta_files:
        log_error(
            f"No *-04-chimera-checked.fasta.gz files found under: "
            f"{samples_dir}/*/04-chimera-checked/"
        )
        sys.exit(1)

    log(f"  Found {len(fasta_files)} sample file(s):")
    for f in fasta_files:
        log(f"    {f}")

    # Concatenate all samples into a single FASTA (decompress on the fly)
    all_fasta = output_dir / "all_samples.fasta"
    with open(all_fasta, "w") as out_fh:
        for fasta_gz in fasta_files:
            with gzip.open(fasta_gz, "rt") as f_in:
                shutil.copyfileobj(f_in, out_fh)

    n_total = count_fasta_seqs(all_fasta)

    compress(all_fasta, all_fasta_gz)

    log(f"  Total sequences pooled: {n_total}")

    ###########################################################################
    # Step 1: OTU clustering
    ###########################################################################

    log(f"Step 1/1: Clustering OTUs at {identity} identity...")

    otus_fasta    = output_dir / "otus" / "otus.fasta"
    otus_fasta_gz = output_dir / "otus" / "otus.fasta.gz"
    otu_table     = output_dir / "otus" / "otu_table.tsv"

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
    run(cluster_cmd, log_path=output_dir / "logs" / "cluster.log")

    n_otus = count_fasta_seqs(otus_fasta)

    compress(otus_fasta, otus_fasta_gz)

    log(f"  OTUs generated: {n_otus}")

    log(f"  OTU table written: {otu_table}")

    otu_rows = [
        ["step",          "count"],
        ["pooled_seqs",   str(n_total)],
        ["otus",          str(n_otus)],
    ]
    (output_dir / "stats" / "otu_summary.tsv").write_text(
        "\n".join("\t".join(r) for r in otu_rows) + "\n"
    )

    ###########################################################################
    # Summary report
    ###########################################################################

    log("Generating summary report...")
    summary_txt = output_dir / "summary_report.txt"

    report = (
        f"{'=' * 80}\n"
        f"VSEARCH OTU Construction Report\n"
        f"{'=' * 80}\n"
        f"Date: {datetime.now()}\n"
        f"Samples directory: {samples_dir}\n"
        f"Output directory:  {output_dir}\n"
        f"Samples processed: {len(fasta_files)}\n"
        f"\n"
        f"Parameters:\n"
        f"-----------\n"
        f"  Threads:      {nslots}\n"
        f"  OTU identity: {identity}\n"
        f"\n"
        f"Output files:\n"
        f"-------------\n"
        f"  Concatenated FASTA:  {all_fasta_gz}\n"
        f"  OTU representatives: {otus_fasta_gz}\n"
        f"  OTU table:           {otu_table}\n"
        f"  Processing logs:     {output_dir}/logs/\n"
        f"\n"
        f"{'=' * 80}\n"
        f"\n"
        f"OTU statistics:\n"
        f"\n"
        f"{fmt_tsv(otu_rows)}\n"
    )

    summary_txt.write_text(report)
    # print(report)
    log("\033[0;32m2.2.2-vsearch-pipeline.py completed successfully\033[0m")

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
