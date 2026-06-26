#!/usr/bin/env python3

################################################################################
# 1. Set env
################################################################################

import argparse
import csv
import gzip
import re
import sys
from pathlib import Path

# DEV ONLY — comment out before production use
"""
otus_fasta = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch/otu/otus/otus.fasta.gz"
otu_table  = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch/otu/otus/otu_table.tsv"
output     = "/home/epereira/workspace/repos/tools/Amp-Proc/tests/output/05_vsearch/otu/otus/seqtable.csv"
"""

################################################################################
# 2. Define functions
################################################################################

def log(msg):
    print(f"[INFO] {msg}")

def log_error(msg):
    print(f"\033[0;31m[ERROR]\033[0m {msg}", file=sys.stderr)

def parse_args():
    p = argparse.ArgumentParser(
        description="Bridge VSEARCH OTU outputs into a DADA2-style sequence-keyed "
                    "count table so it can be fed to 3-taxa_annot.R"
    )
    p.add_argument("--otus_fasta", required=True, help="OTU centroid FASTA (gzip or plain) (required)")
    p.add_argument("--otu_table",  required=True, help="OTU count table from vsearch (--otutabout) (required)")
    p.add_argument("-o", "--output", required=True, help="Output sequence-keyed CSV (required)")
    return p.parse_args()

def open_maybe_gzip(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "r")

def read_centroids(fasta_path):
    """Return {otu_id: sequence}. Header size= annotations are stripped from the id."""
    seqs = {}
    current = None
    chunks = []
    with open_maybe_gzip(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current is not None:
                    seqs[current] = "".join(chunks)
                current = re.sub(r";size=\d+;?$", "", line[1:].split()[0])
                chunks = []
            else:
                chunks.append(line)
    if current is not None:
        seqs[current] = "".join(chunks)
    return seqs

################################################################################
# 3. Define the main function
################################################################################

def main():
    opts = parse_args()

    if not Path(opts.otus_fasta).is_file():
        log_error(f"OTU FASTA does not exist: {opts.otus_fasta}")
        sys.exit(1)
    if not Path(opts.otu_table).is_file():
        log_error(f"OTU table does not exist: {opts.otu_table}")
        sys.exit(1)

    centroids = read_centroids(opts.otus_fasta)
    log(f"Loaded {len(centroids)} OTU centroid sequences")

    with open(opts.otu_table, "r", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]

        # DADA2-style: empty first header, then sample columns; rows keyed by sequence
        with open(opts.output, "w", newline="") as out_fh:
            writer = csv.writer(out_fh)
            writer.writerow([""] + samples)

            n_written = 0
            n_missing = 0
            for row in reader:
                if not row:
                    continue
                otu_id = row[0]
                seq = centroids.get(otu_id)
                if seq is None:
                    n_missing += 1
                    continue
                writer.writerow([seq] + row[1:])
                n_written += 1

    if n_missing:
        log_error(f"{n_missing} OTU id(s) in the table had no matching centroid sequence")
    log(f"Wrote {n_written} sequence-keyed rows x {len(samples)} samples to: {opts.output}")
    log("\033[0;32m2.2.3-otu-to-seqtable.py completed successfully\033[0m")

################################################################################
# 4. Execute
################################################################################

if __name__ == "__main__":
    main()
