from pathlib import Path
from Bio import SeqIO
import argparse
import re
import csv

# --- ARGUMENT PARSING ---
parser = argparse.ArgumentParser(description="Compare KEGG and SwissProt annotations and separate hits.")
parser.add_argument("-b", "--blast", required=True, help="Path to BLASTP output file (outfmt 7)")
parser.add_argument("-k", "--kegg", required=True, help="Path to KEGG .top result file")
parser.add_argument("-f", "--fasta", required=True, help="Path to original multifasta file (from BRAKER3)")
parser.add_argument("-o", "--out", default="strain", help="Output name prefix")
args = parser.parse_args()

# --- FILE PATHS ---
blastp_file = Path(args.blast)
kegg_file = Path(args.kegg)
fasta_file = Path(args.fasta)
prefix = Path(args.out).stem

# --- STEP 1: Parse KEGG ---
kegg_hits = {}
with kegg_file.open() as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 7:
            continue
        query = parts[0].replace("user:", "")
        ko = parts[1]
        taxon = parts[2:6]
        score = float(parts[6])
        if ko and score >= 50:
            kegg_hits[query] = {
                "ko": ko,
                "taxon": " | ".join(taxon),
                "score": score
            }

# --- STEP 2: Parse BLASTP (outfmt 7) ---
blast_hits = {}
with blastp_file.open() as f:
    current_query = None
    for line in f:
        line = line.strip()
        if line.startswith("# Query:"):
            current_query = line.split()[-1]
        elif line.startswith("#") or not line:
            continue
        elif current_query and current_query not in blast_hits:
            cols = line.split('\t')
            subject = cols[0]
            evalue = cols[6]
            bitscore = cols[7]
            subject_title = cols[8] if len(cols) > 8 else "NA"
            match = re.search(r"OS=([^=]+?) OX=", subject_title)
            taxon = match.group(1).strip() if match else "Unknown"
            blast_hits[current_query] = {
                "subject": subject,
                "taxon": taxon,
                "bitscore": float(bitscore),
                "evalue": evalue
            }

# --- STEP 3: Get Intersections ---
queries_kegg = set(kegg_hits.keys())
queries_blast = set(blast_hits.keys())
intersection = queries_kegg & queries_blast
only_kegg = queries_kegg - intersection
only_blast = queries_blast - intersection

# All queries from the fasta
all_queries = set(rec.id for rec in SeqIO.parse(fasta_file, "fasta"))
no_hits = all_queries - (queries_kegg | queries_blast)

# --- STEP 4: Load full sequences from multifasta ---
seq_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

# --- OUTPUT FUNCTION ---
def write_csv(file_path, header, rows):
    with open(file_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(header)
        writer.writerows(rows)

def write_fasta(file_path, ids):
    with open(file_path, "w") as f:
        SeqIO.write((seq_dict[qid] for qid in ids if qid in seq_dict), f, "fasta")

# --- STEP 5: Write KEGG Hits ---
kegg_csv = f"{prefix}_kegg.csv"
kegg_faa = f"{prefix}_kegg.faa"
write_csv(kegg_csv, ["query", "KO", "taxon", "score"], [
    [qid, kegg_hits[qid]["ko"], kegg_hits[qid]["taxon"], kegg_hits[qid]["score"]] for qid in only_kegg
])
write_fasta(kegg_faa, only_kegg)

# --- STEP 6: Write SwissProt Hits ---
sprot_csv = f"{prefix}_sprot.csv"
sprot_faa = f"{prefix}_sprot.faa"
write_csv(sprot_csv, ["query", "subject", "taxon", "bitscore", "evalue"], [
    [qid, blast_hits[qid]["subject"], blast_hits[qid]["taxon"], blast_hits[qid]["bitscore"], blast_hits[qid]["evalue"]] for qid in only_blast
])
write_fasta(sprot_faa, only_blast)

# --- STEP 7: Write Intersection ---
inter_csv = f"{prefix}_intersection.csv"
inter_faa = f"{prefix}_intersection.faa"
write_csv(inter_csv, ["query", "KO", "KEGG_taxon", "KEGG_score", "SwissProt_subject", "SwissProt_taxon", "bitscore", "evalue"], [
    [
        qid,
        kegg_hits[qid]["ko"],
        kegg_hits[qid]["taxon"],
        kegg_hits[qid]["score"],
        blast_hits[qid]["subject"],
        blast_hits[qid]["taxon"],
        blast_hits[qid]["bitscore"],
        blast_hits[qid]["evalue"]
    ] for qid in intersection
])
write_fasta(inter_faa, intersection)

# --- STEP 8: Write No Hits ---
nohit_txt = f"{prefix}_nohits.txt"
nohit_faa = f"{prefix}_nohits.faa"
with open(nohit_txt, "w") as f:
    f.write("\n".join(sorted(no_hits)))
write_fasta(nohit_faa, no_hits)

print("All files written! :D")
