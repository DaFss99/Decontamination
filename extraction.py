from pathlib import Path
from Bio import SeqIO
import argparse
import re

# --- ARGUMENT PARSING ---
parser = argparse.ArgumentParser(description="Extract sequences with no BLAST hits.")
parser.add_argument("-b", "--blast", required=True, help="Path to BLASTP output file")
parser.add_argument("-p", "--prot", required=True, help="Path to protein multifasta file")
parser.add_argument("-o", "--out", default="noHitsSeq.faa", help="Output fasta file")
parser.add_argument("-l", "--log", default="counting.txt", help="Log output file")

args = parser.parse_args()

# --- FILE PATHS ---
blastp_file = Path(args.blast)
fasta_file = Path(args.prot)
output_file = Path(args.out)
log_file = Path(args.log)

checkAll = 0
checkZero = 0
checkHit = 0

# --- STEP 1: Get query IDs with 0 hits ---
no_hit_ids = set()
with blastp_file.open() as blast:
    current_query = None
    for line in blast:
        if line.startswith("# Query:"):
            checkAll += 1
            current_query = line.strip().split()[-1]
        elif "# 0 hits found" in line and current_query:
            checkZero += 1
            no_hit_ids.add(current_query)

    checkHit = checkAll - checkZero

# --- STEP 2: Extract matching sequences from FASTA ---
with output_file.open("w") as out_fasta:
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in no_hit_ids:
            SeqIO.write(record, out_fasta, "fasta")

# --- STEP 3: Write log ---
log_text = f"{'Total number of queries: '} {checkAll}\n\
{'Number of queries with some hit: '} {checkHit}\n\
{'Number of queries with zero hits'} {checkZero}\n"
counting = "counting.txt"

with open(log_file, "w") as file:
    file.write(log_text)


print(f"Done! Extracted {len(no_hit_ids)} sequences to '{output_file.name}'")
print(f"Log written to '{log_file.name}'")