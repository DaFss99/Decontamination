from pathlib import Path
from Bio import SeqIO

# --- FILE PATHS ---
blastp_file = Path("blastTest.txt")
fasta_file = Path("braker.aa")
output_file = Path("noHitsSeq.faa")

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

print(f"✔ Done! Extracted {len(no_hit_ids)} sequences to '{output_file.name}'")

print(checkAll, checkHit, checkZero)