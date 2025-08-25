#!/usr/bin/env python3
import argparse
import subprocess
from pathlib import Path

# --- PARSE ARGUMENTS ---
parser = argparse.ArgumentParser(
    description="Python script to run Diamond BLASTp search"
)
parser.add_argument(
    "-f", "--fasta",
    required=True, help="Path to protein multifasta file (e.g.: path/to/my_file.fasta)"
)
parser.add_argument(
    "-d", "--database",
    required=True, help="Path to reference database (e.g.: NR, SwissProt)"
)
parser.add_argument(
    "-t", "--threads",
    required=True, type=int, help="Number of threads for running the script."
)
parser.add_argument(
    "-o", "--output",
    required=True, help="Prefix for output files (e.g.: strain_name_database, strain_code_database)"
)
args = parser.parse_args()

# --- FILE PATHS ---
fasta = Path(args.fasta)
db = Path(args.database)
threads = args.threads  # This should be an integer, not a Path object
output_prefix = args.output
output_file = f"{output_prefix}_Dblastp.txt"

# --- DIAMOND COMMAND ---
# Format the outfmt as comma-separated without spaces
outfmt_fields = "sseqid qstart qend sstart send qseq evalue bitscore stitle"

diamond_command = [
    "diamond", "blastp",
    "--query", str(fasta),
    "--db", str(db),
    "--out", output_file,
    "--outfmt", "6", "sseqid", "qstart", "qend", "sstart", "send", "qseq", "evalue", "bitscore", "stitle",
    "--evalue", "1e-5",
    "--threads", str(threads),
    "--max-target-seqs", "5"
]

# --- RUN DIAMOND ---
try:
    subprocess.run(diamond_command, check=True)
    print(f"BLASTp search completed. Results saved to {output_file}")
except subprocess.CalledProcessError as e:
    print(f"Error running Diamond: {e}")
    exit(1)
