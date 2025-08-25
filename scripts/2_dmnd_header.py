#!/usr/bin/env python3
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
import argparse
import re
import sys


# usage: python headerize.py queries.faa diamond.tsv > diamond_outfmt7.txt

# --- PARSE ARGUMENTS ---
parser = argparse.ArgumentParser(
    description="Python headerize.py queries.faa diamond.tsv > diamond_outfmt7.txt.")
parser.add_argument(
    "-f", "--fasta", 
    required=True, help="Path to protein multifasta file (e.g.: path/to/my_file.fasta)")
parser.add_argument(
    "-d", "--dmnd",
    required=True, help="Path to diamond output file (e.g.:path/to/my__Dblastp.txt)")
parser.add_argument(
    "-o", "--output",
    required=True, help="Name tag for output files (e.g.: strain_name, strain_code, test_id)")
args = parser.parse_args()

# --- FILE PATHS ---
fasta = Path(args.fasta)
dmnd = Path(args.dmnd)
prefix = Path(args.output)
output_file = Path(args.output).with_name(Path(args.output).stem + "_outfmt7.txt")


# read all query IDs from FASTA

all_qids = [rec.id for rec in SeqIO.parse(fasta, "fasta")]

# read DIAMOND tsv
hits = defaultdict(list)
with open(dmnd) as f:
    for line in f:
        if not line.strip(): continue
        cols = line.rstrip("\n").split("\t")
        # sseqid qstart qend sstart send qseq evalue bitscore stitle
        qid = cols[0]
        hits[qid].append(cols)

# write outfmt7-like
with open(output_file, "w") as out:
    out.write("# DIAMOND blastp (outfmt7-like)\n")
    out.write("# Fields: subject id, q. start, q. end, s. start, s. end, query seq, evalue, bit score, subject title\n")
    for qid in all_qids:
        out.write(f"# Query: {qid}\n")
        if qid not in hits or len(hits[qid]) == 0:
            out.write("# 0 hits found\n")
            continue
        out.write(f"# {len(hits[qid])} hits found\n")
        for cols in hits[qid]:
            sseqid   = cols[1]
            qstart    = cols[2]
            qend      = cols[3]
            sstart    = cols[4]
            send      = cols[5]
            qseq      = cols[6]
            evalue    = cols[7]
            bitscore  = cols[8]
            stitle    = cols[9]
            out.write("\t".join([sseqid, qstart, qend, sstart, send, qseq, evalue, bitscore, stitle]) + "\n")
print(f"Results written to {output_file}! :D")

