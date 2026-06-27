#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import csv
import gzip
import os
import re
import random
import sys
import shutil
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData

# DEV ONLY — comment out before production use
""" 
reads1 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
reads2 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
output_dir = Path("/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/02_check_primers_before/sample1")
sample_name = "1-samo1_S1"
primer_fwd = "GTGYCAGCMGCCGCGGTAA"
primer_rev = "CCGYCAATTYMTTTRAGTTT"
subsample_size = 100
 """
################################################################################
# 2. Define functions
################################################################################

_log_buffer = []
_ANSI_RE = re.compile(r'\033\[[0-9;]*m')

def log(msg):
    line = f"[INFO] {msg}"
    print(line)
    _log_buffer.append(_ANSI_RE.sub('', line))

def log_warn(msg):
    line = f'[WARN] {msg}'
    print(f"\033[1;33m{line}\033[0m", file=sys.stderr)
    _log_buffer.append(_ANSI_RE.sub('', line))

def log_error(msg):
    line = f"[ERROR] {msg}"
    print(f"\033[0;31m{line}\033[0m", file=sys.stderr)
    _log_buffer.append(_ANSI_RE.sub('', line))


def fmt_tsv(rows):
    col_w = [max(len(r[i]) for r in rows) for i in range(len(rows[0]))]
    return "\n".join(
        "  ".join(r[i].ljust(col_w[i]) for i in range(len(r))) for r in rows
    )


def parse_args():
    p = argparse.ArgumentParser(
        description="Subsample reads and count IUPAC-aware primer hits"
    )
    p.add_argument("--reads1",          required=True,  help="R1 FASTQ file (required)")
    p.add_argument("--reads2",          required=True,  help="R2 FASTQ file (required)")
    p.add_argument("--output_dir",      required=True, help="Path to output directory")
    p.add_argument("--sample_name",     default=None,   help="Sample name [default: derived from R1 filename]")
    p.add_argument("--primer_fwd",      required=True,  help="Forward primer sequence")
    p.add_argument("--primer_rev",      required=True,  help="Reverse primer sequence")
    p.add_argument("--subsample_size",  type=int,  default=1000, help="Reads to subsample per file")
    p.add_argument("--raw_counts",      action="store_true", help="Write raw counts instead of percentages (default: percentages)")
    p.add_argument("--overwrite",       choices=["t", "f"], default="f", help="Overwrite existing output directory [default=f]")
    return p.parse_args()

def all_orients(primer: str):
    """Return dict of all IUPAC orientations of a DNA primer."""
    seq = Seq(primer)
    return {
        "Forward":    str(seq),
        "Complement": str(seq.complement()),
        "Reverse":    str(seq[::-1]),
        "RevComp":    str(seq.reverse_complement()),
    }

iupac_regex = {
    **{k: k for k in "ACGT"},
    **{amb: "[" + "".join(v) + "]"
       for amb, v in IUPACData.ambiguous_dna_values.items()}
}

def primer_to_regex(primer: str):
    """Convert IUPAC primer sequence to a regex pattern (case-insensitive)."""
    pattern = "".join(iupac_regex.get(nt, nt) for nt in primer.upper())
    return re.compile(pattern, re.IGNORECASE)

def count_hits(seqs, regex):
    """Count how many sequences have >=1 match to regex."""
    return sum(1 for s in seqs if regex.search(str(s.seq)))

def open_maybe_gzip(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")

def write_table(data_dict, out_file):
    rows = {}
    for key, val in data_dict.items():
        row, col = key.rsplit(".", 1)
        rows.setdefault(row, {})[col] = val
    with open(out_file, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        header = [""] + list(next(iter(rows.values())).keys())
        writer.writerow(header)
        for row_name, coldict in rows.items():
            writer.writerow([row_name] + [coldict[c] for c in header[1:]])

################################################################################
# 3. Define the main function
################################################################################

def main():

    opts = parse_args()
    reads1 = opts.reads1
    reads2 = opts.reads2
    output_dir = Path(opts.output_dir)
    sample_name = opts.sample_name
    primer_fwd = opts.primer_fwd
    primer_rev = opts.primer_rev
    subsample_size = opts.subsample_size
    raw_counts = opts.raw_counts
    overwrite = opts.overwrite == "t"

    ###########################################################################
    # Step 0: Validate inputs and create output directories
    ###########################################################################

    # Validate input files
    if not os.path.isfile(reads1):
        log_error(f"R1 file does not exist: {reads1}")
        sys.exit(1)
    if not os.path.isfile(reads2):
        log_error(f"R2 file does not exist: {reads2}")
        sys.exit(1)

    # Derive sample name from R1 filename if not provided
    if sample_name is None:
        sample_name = Path(reads1).name
        for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
            if sample_name.endswith(ext):
                sample_name = sample_name[: -len(ext)]
                break

    # Check if output directory exists and handle overwrite option
    if output_dir.exists():
        if overwrite:
            log_warn(f"Overwriting existing directory: {output_dir}")
            shutil.rmtree(output_dir)
        else:
            log_error(f"Output directory exists: {output_dir}. Use --overwrite t to overwrite.")
            sys.exit(1)

    # Create output directories
    primer_check_dir = output_dir / "output"
    logs_dir  = output_dir / "logs"
    output_dir.mkdir(parents=True, exist_ok=True)
    primer_check_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    primer_check_out = primer_check_dir / f"{sample_name}_primer-check.tsv"
    log_out  = logs_dir  / f"{sample_name}_primers-check.log"

    ###########################################################################
    # Step 1: Obtain all orientations of the primers and compile regex patterns
    ###########################################################################

    log(f"Processing sample: {sample_name}")
    log(f"Forward primer: {primer_fwd}")
    log(f"Reverse primer: {primer_rev}")

    # Prepare primer orientations and regexes
    fwd_orients = all_orients(primer_fwd)
    rev_orients = all_orients(primer_rev)
    fwd_regex = {name: primer_to_regex(seq) for name, seq in fwd_orients.items()}
    rev_regex = {name: primer_to_regex(seq) for name, seq in rev_orients.items()}

    # Log primer sequences in all orientations (captured into the log file)
    log("Primer sequences (all orientations):")
    for orient_name, seq in fwd_orients.items():
        log(f"  FwdPrimer.{orient_name}: {seq}")
    for orient_name, seq in rev_orients.items():
        log(f"  RevPrimer.{orient_name}: {seq}")

    ###########################################################################
    # Step 2: Read and subsample
    ###########################################################################

    with open_maybe_gzip(reads1) as handle:
        seqs1 = list(SeqIO.parse(handle, "fastq"))
    with open_maybe_gzip(reads2) as handle:
        seqs2 = list(SeqIO.parse(handle, "fastq"))

    if len(seqs1) > subsample_size:
        random.seed(123)
        seqs1 = random.sample(seqs1, subsample_size)
    if len(seqs2) > subsample_size:
        random.seed(123)
        seqs2 = random.sample(seqs2, subsample_size)

    log(f"Subsampled reads: R1={len(seqs1)} R2={len(seqs2)} (cap {subsample_size})")

    ###########################################################################
    # Step 3: Count primer hits
    ###########################################################################

    hit_counts = {}
    hit_counts.update({
        f"FwdReads.FwdPrimer.{orient}": count_hits(seqs1, regex)
        for orient, regex in fwd_regex.items()
    })
    hit_counts.update({
        f"RevReads.FwdPrimer.{orient}": count_hits(seqs2, regex)
        for orient, regex in fwd_regex.items()
    })
    hit_counts.update({
        f"FwdReads.RevPrimer.{orient}": count_hits(seqs1, regex)
        for orient, regex in rev_regex.items()
    })
    hit_counts.update({
        f"RevReads.RevPrimer.{orient}": count_hits(seqs2, regex)
        for orient, regex in rev_regex.items()
    })

    ###########################################################################
    # Step 4: Calculate percentages
    ###########################################################################

    if raw_counts:
        write_table(hit_counts, primer_check_out)
        display_vals = hit_counts
    else:
        perc = {k: round(v / subsample_size * 100, 2) for k, v in hit_counts.items()}
        write_table(perc, primer_check_out)
        display_vals = perc

    ###########################################################################
    # Step 5: Write summary report
    ###########################################################################

    log("Generating summary report...")

    # Key orientations: expected signal in a well-prepared library
    key_keys = [
        "FwdReads.FwdPrimer.Forward",
        "FwdReads.RevPrimer.RevComp",
        "RevReads.RevPrimer.Forward",
        "RevReads.FwdPrimer.RevComp",
    ]
    unit = "counts" if raw_counts else "%"
    stat_rows = [["orientation", unit]] + [
        [k, str(display_vals[k])] for k in key_keys if k in display_vals
    ]

    report_out = output_dir / f"{sample_name}_summary_report.txt"
    report = (
        f"{'=' * 80}\n"
        f"Primer Check Report\n"
        f"{'=' * 80}\n"
        f"Date:             {datetime.now()}\n"
        f"Sample:           {sample_name}\n"
        f"R1:               {reads1}\n"
        f"R2:               {reads2}\n"
        f"Output directory: {output_dir}\n"
        f"\n"
        f"Parameters:\n"
        f"-----------\n"
        f"  Forward primer:  {primer_fwd}\n"
        f"  Reverse primer:  {primer_rev}\n"
        f"  Subsample size:  {subsample_size}\n"
        f"  Output format:   {'raw counts' if raw_counts else 'percentages'}\n"
        f"\n"
        f"Output files:\n"
        f"-------------\n"
        f"  Statistics:           {primer_check_out}\n"
        f"  Log (incl. primers):  {log_out}\n"
        f"\n"
        f"{'=' * 80}\n"
        f"\n"
        f"Key primer hit statistics:\n"
        f"\n"
        f"{fmt_tsv(stat_rows)}\n"
    )

    report_out.write_text(report)
    # print(report)
    log("\033[0;32m1.2-primers-check.py completed successfully\033[0m")
    log_out.write_text("\n".join(_log_buffer) + "\n")

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
