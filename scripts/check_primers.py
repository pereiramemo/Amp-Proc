#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import os
import re
import random
import csv
import gzip
import io
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import IUPACData
import numpy as np

# Dev only
# input_dir = "/home/epereira/workspace/nal-case-studies/Sildever_2023_NW_Pacific_Metabarcoding/data/DRA016526/"
# output_dir = "/home/epereira/workspace/nal-case-studies/Sildever_2023_NW_Pacific_Metabarcoding/results/DRA016526/"
# suffix_r1 = "_1.fastq.gz"
# suffix_r2 = "_2.fastq.gz"
# primer_fwd = "CAAGTACCATGAGGGAAAG"
# primer_rev = "GACTCCTTGGTCCGTGTTTC"
# subsample_size = 100

################################################################################
# 2. Define functions
################################################################################

def parse_args():
    p = argparse.ArgumentParser(
        description="Subsample reads and count IUPAC‐aware primer hits"
    )
    p.add_argument("-i", "--input_dir",    required=True,  help="Path to input directory")
    p.add_argument("-o", "--output_dir",   required=True,  help="Path to output directory")
    p.add_argument("--suffix_r1",          required=True,  help="Glob pattern for forward reads (R1)")
    p.add_argument("--suffix_r2",          required=True,  help="Glob pattern for reverse reads (R2)")
    p.add_argument("--primer_fwd",         required=True,  help="Forward primer sequence")
    p.add_argument("--primer_rev",         required=True,  help="Reverse primer sequence")
    p.add_argument("-s", "--subsample_size",type=int,     default=1000, help="Reads to subsample per file")
    p.add_argument("--counts", action="store_true", help="Write raw counts instead of percentages (default: percentages)")
    return p.parse_args()

def all_orients(primer: str):
    """Return dict of all IUPAC orientations of a DNA primer."""
    seq = Seq(primer)
    orients = {
        "Forward":    str(seq),
        "Complement": str(seq.complement()),
        "Reverse":    str(seq[::-1]),
        "RevComp":    str(seq.reverse_complement()),
    }
    return orients

# Build regex from IUPAC to match ambiguity codes
iupac_regex = {
    **{k: k for k in "ACGT"},  # exact
    **{amb: "[" + "".join(v) + "]" 
       for amb, v in IUPACData.ambiguous_dna_values.items()}
}

def primer_to_regex(primer: str):
    """Convert IUPAC primer sequence to a regex pattern (case-insensitive)."""
    pattern = "".join(iupac_regex.get(nt, nt) for nt in primer.upper())
    return re.compile(pattern, re.IGNORECASE)

def count_hits(seqs, regex):
    """Count how many sequences have ≥1 match to regex."""
    return sum(1 for s in seqs if regex.search(str(s.seq)))
  
def open_maybe_gzip(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    else:
        return open(path, "r")  

def write_table(data_dict, out_file):
    # split keys into row, col
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
  
    # Parse command line arguments
    opts = parse_args()
    input_dir = opts.input_dir
    output_dir = opts.output_dir
    primer_fwd = opts.primer_fwd
    primer_rev = opts.primer_rev
    suffix_r1 = opts.suffix_r1
    suffix_r2 = opts.suffix_r2
    subsample_size = opts.subsample_size
    counts = opts.counts
    
    # Validate input directory
    if not os.path.isdir(input_dir):
        raise ValueError(f"Input directory '{input_dir}' does not exist or is not a directory.")

    # list files
    r1_files = sorted(Path(input_dir).glob("*" + suffix_r1))
    r2_files = sorted(Path(input_dir).glob("*" + suffix_r2))
    sample_names = [f.name.replace(suffix_r1, "") for f in r1_files]

    # check that at least one file is found
    if len(r1_files) == 0 or len(r2_files) == 0:
        raise ValueError("No input files found with the specified suffixes.")
    if len(r1_files) != len(r2_files):
        raise ValueError("The number of R1 and R2 files do not match.")
    
    # prepare primer orientations and regexes
    fwd_orients = all_orients(primer_fwd)
    rev_orients = all_orients(primer_rev)
    fwd_regex = {name: primer_to_regex(seq) for name, seq in fwd_orients.items()}
    rev_regex = {name: primer_to_regex(seq) for name, seq in rev_orients.items()}

    # Create output dir
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)
    
    # Write primer sequences file (once for all samples)
    primer_file = output_dir_path / "primer_sequences.csv"
    write_primer_sequences(fwd_orients, rev_orients, primer_file)

    # loop samples
    for r1, r2, sample in zip(r1_files, r2_files, sample_names):
      
        # read all records
        with open_maybe_gzip(r1) as handle:
            seqs1 = list(SeqIO.parse(handle, "fastq"))
        with open_maybe_gzip(r2) as handle:
            seqs2 = list(SeqIO.parse(handle, "fastq"))

        # subsample
        if len(seqs1) > subsample_size:
            seqs1 = random.sample(seqs1, subsample_size)
        if len(seqs2) > subsample_size:
            seqs2 = random.sample(seqs2, subsample_size)

        # count occurrences
        counts = {
            f"FwdReads.FwdPrimer.{orient}": count_hits(seqs1, regex)
            for orient, regex in fwd_regex.items()
        }
        counts.update({
            f"RevReads.FwdPrimer.{orient}": count_hits(seqs2, regex)
            for orient, regex in fwd_regex.items()
        })
        counts.update({
            f"FwdReads.RevPrimer.{orient}": count_hits(seqs1, regex)
            for orient, regex in rev_regex.items()
        })
        counts.update({
            f"RevReads.RevPrimer.{orient}": count_hits(seqs2, regex)
            for orient, regex in rev_regex.items()
        })
        
        # Define output file name
        output_file = output_dir_path / f"{sample}_primer_check.csv"

        if counts is not True:
            # compute and write percentages
            perc = {k: round(v / subsample_size * 100, 2) for k, v in counts.items()}
            write_table(perc, output_file)
        else:
            # write raw counts
            write_table(counts, output_file)

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
