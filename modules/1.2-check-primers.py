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
from datetime import datetime
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData

# DEV ONLY — comment out before production use
""" reads1 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R1_001_redu.fastq.gz"
reads2 = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/data/1-samo1_S1_L001_R2_001_redu.fastq.gz"
output_dir = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/02_check_primers_before/sample1"
primer_fwd = "GTGYCAGCMGCCGCGGTAA"
primer_rev = "CCGYCAATTYMTTTRAGTTT"
subsample_size = 100
 """
################################################################################
# 2. Define functions
################################################################################

def log(msg):
    print(f"[INFO] {msg}")

def log_error(msg):
    print(f"\033[0;31m[ERROR]\033[0m {msg}", file=sys.stderr)


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
    p.add_argument("-o", "--output_dir", required=True, help="Path to output directory")
    p.add_argument("--primer_fwd",      required=True,  help="Forward primer sequence")
    p.add_argument("--primer_rev",      required=True,  help="Reverse primer sequence")
    p.add_argument("-s", "--subsample_size", type=int,  default=1000, help="Reads to subsample per file")
    p.add_argument("--raw_counts", action="store_true", help="Write raw counts instead of percentages (default: percentages)")
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
        writer = csv.writer(fh)
        header = [""] + list(next(iter(rows.values())).keys())
        writer.writerow(header)
        for row_name, coldict in rows.items():
            writer.writerow([row_name] + [coldict[c] for c in header[1:]])

def write_primer_sequences(fwd_orients, rev_orients, out_file):
    """Write primer sequences in all orientations to a CSV file."""
    with open(out_file, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["Primer", "Orientation", "Sequence"])
        for orient_name, seq in fwd_orients.items():
            writer.writerow(["Forward", orient_name, seq])
        for orient_name, seq in rev_orients.items():
            writer.writerow(["Reverse", orient_name, seq])

################################################################################
# 3. Define the main function
################################################################################

def main():

    opts = parse_args()
    reads1 = opts.reads1
    reads2 = opts.reads2
    output_dir = opts.output_dir
    primer_fwd = opts.primer_fwd
    primer_rev = opts.primer_rev
    subsample_size = opts.subsample_size
    raw_counts = opts.raw_counts

    # Validate input files
    if not os.path.isfile(reads1):
        log_error(f"R1 file does not exist: {reads1}")
        sys.exit(1)
    if not os.path.isfile(reads2):
        log_error(f"R2 file does not exist: {reads2}")
        sys.exit(1)

    # Derive sample name from R1 filename
    sample_name = Path(reads1).name
    for ext in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if sample_name.endswith(ext):
            sample_name = sample_name[: -len(ext)]
            break

    # Prepare primer orientations and regexes
    fwd_orients = all_orients(primer_fwd)
    rev_orients = all_orients(primer_rev)
    fwd_regex = {name: primer_to_regex(seq) for name, seq in fwd_orients.items()}
    rev_regex = {name: primer_to_regex(seq) for name, seq in rev_orients.items()}

    # Create output directory
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    # Write primer sequences (once per run)
    primer_file = output_dir_path / "primer_sequences.csv"
    write_primer_sequences(fwd_orients, rev_orients, primer_file)

    # Read and subsample
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

    # Count primer hits
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

    output_file = output_dir_path / f"{sample_name}_primer_check.csv"

    if raw_counts:
        write_table(hit_counts, output_file)
        display_vals = hit_counts
    else:
        perc = {k: round(v / subsample_size * 100, 2) for k, v in hit_counts.items()}
        write_table(perc, output_file)
        display_vals = perc

    ###########################################################################
    # Summary report
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

    summary_txt = output_dir_path / "summary_report.txt"
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
        f"  Primer hits table:    {output_file}\n"
        f"  Primer sequences:     {primer_file}\n"
        f"\n"
        f"{'=' * 80}\n"
        f"\n"
        f"Key primer hit statistics:\n"
        f"\n"
        f"{fmt_tsv(stat_rows)}\n"
    )

    summary_txt.write_text(report)
    # print(report)
    log("\033[0;32m1.2-check-primers.py completed successfully\033[0m")

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
