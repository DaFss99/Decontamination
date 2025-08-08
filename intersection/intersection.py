from pathlib import Path
from Bio import SeqIO
import argparse
import re
import csv

# --- ARGUMENT PARSING ---
parser = argparse.ArgumentParser(description="Compare KEGG and SwissProt annotations and separate hits.")
parser.add_argument("-s", "--swissprot", required=True, help="Path to BLASTP output file of nohits against SwissProt (outfmt 7)")
parser.add_argument("-n", "--nr", required=True, help="Path to BLASTP output file of nohits against NR (outfmt 7)")
parser.add_argument("-k", "--kegg", required=True, help="Path to KEGG .top result file")
parser.add_argument("-ks", "--kegg_score", required=True, help="Chosed filtering kegg score.")
parser.add_argument("-f", "--fasta", required=True, help="Path to nohits multifasta file (after blastp agains type organism proteome)")
parser.add_argument("-o", "--out", default="strain", help="Output name prefix")
args = parser.parse_args()

# --- FILE PATHS ---
log_file = Path(args.out).with_name(Path(args.out).stem + ".log")
swissprot_file = Path(args.swissprot)   
nr_file = Path(args.nr)
kegg_file = Path(args.kegg)
kegg_filter = float(args.kegg_score)
fasta_file = Path(args.fasta)
prefix = Path(args.out).stem

# --- CHECKING VARIABLES ---
check_kegg = 0
check_blast = 0
check_nr = 0
check_intersec = 0
check_onlykegg = 0
check_onlyblast = 0
check_onlynr = 0
check_noHits = 0
check_filter_kegg = 0
check_filter_blast = 0
check_filter_nr = 0

# --- STEP 1: Parse KEGG ---
kegg_hits = {}
with kegg_file.open() as f:
    for line in f:
        check_kegg += 1     # for log file
        parts = line.strip().split('\t')
        query = parts[0].replace("user:", "").strip()
        ko = parts[1].strip() if len(parts) > 1 else ""
        taxon = parts[2:6]
        score = float(parts[6]) if len(parts) > 6 else 0.0
        if score >= kegg_filter:
            check_filter_kegg += 1
            kegg_hits[query] = {
                "ko": ko if ko else "None",
                "taxon": " | ".join(taxon),
                "score": score
            }

# --- STEP 2: Parse SwissProt (outfmt 7) ---
swissprot_hits = {}
with swissprot_file.open() as f:
    current_query = None
    for line in f:
        line = line.strip()
        if line.startswith("# Query:"):
            check_blast += 1     # for log file
            current_query = line.split()[-1]
        elif line.startswith("#") or not line:
            continue
        elif current_query and current_query not in swissprot_hits:
            check_filter_blast += 1
            cols = line.split('\t')
            subject = cols[0]
            evalue = cols[6]
            bitscore = cols[7]
            subject_title = cols[8] if len(cols) > 8 else "NA"
            match = re.search(r"OS=([^=]+?) OX=", subject_title)
            taxon = match.group(1).strip() if match else "Unknown"
            swissprot_hits[current_query] = {
                "subject": subject,
                "taxon": taxon,
                "bitscore": float(bitscore),
                "evalue": evalue
            }

# --- STEP 3: Parse NR (outfmt 7) ---
nr_hits = {}
with nr_file.open() as f:
        current_query = None
        for line in f:
            line = line.strip()
            if line.startswith("# Query:"):
                check_nr += 1
                current_query = line.split()[-1]
            elif line.startswith("#") or not line:
                continue
            elif current_query and current_query not in nr_hits:
                check_filter_nr += 1
                cols = line.split('\t')
                subject = cols[0]
                evalue = cols[6]
                bitscore = cols[7]
                subject_title = cols[8] if len(cols) > 8 else "NA"
                match = re.search(r"\[([^\[\]]+)\]", subject_title)
                taxon = match.group(1).strip() if match else "Unknown"
                nr_hits[current_query] = {
                    "subject": subject,
                    "taxon": taxon,
                    "bitscore": float(bitscore),
                    "evalue": evalue
                }

# --- STEP 4: Get Intersections ---
queries_kegg = set(kegg_hits.keys())
queries_sprot = set(swissprot_hits.keys())
queries_nr = set(nr_hits.keys())
intersection = queries_kegg & queries_sprot & queries_nr
only_kegg = queries_kegg - intersection
only_blast = queries_sprot - intersection
only_nr = queries_nr - intersection

check_intersec = len(intersection)
check_onlyblast = len(only_blast)
check_onlykegg = len(only_kegg)
check_nr = len(only_nr)


# All queries from the fasta
all_queries = set(rec.id for rec in SeqIO.parse(fasta_file, "fasta"))
no_hits = all_queries - (queries_kegg | queries_sprot | queries_nr)

# --- STEP 5: Load full sequences from multifasta ---
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

# --- STEP 6: Write KEGG Hits ---
kegg_csv = f"{prefix}_kegg.csv"
kegg_all = f"{prefix}_filtered_kegg.csv"
kegg_faa = f"{prefix}_kegg.faa"
write_csv(kegg_csv, ["query", "KO", "taxon", "score"], [
    [qid, kegg_hits[qid]["ko"], kegg_hits[qid]["taxon"], kegg_hits[qid]["score"]] for qid in only_kegg
])
write_csv(kegg_all, ["query", "KO", "taxon", "score"], [
    [qid, kegg_hits[qid]["ko"], kegg_hits[qid]["taxon"], kegg_hits[qid]["score"]] for qid in queries_kegg
])
write_fasta(kegg_faa, only_kegg)

# --- STEP 7: Write SwissProt Hits ---
sprot_csv = f"{prefix}_sprot.csv"
sprot_faa = f"{prefix}_sprot.faa"
sprot_all = f"{prefix}_filtered_blast.csv"
write_csv(sprot_csv, ["SwissProt_query", "SwissProt_subject", "SwissProt_taxon", "SwissProt_bitscore", "SwissProt_evalue"], [
    [qid, swissprot_hits[qid]["subject"], swissprot_hits[qid]["taxon"], swissprot_hits[qid]["bitscore"], swissprot_hits[qid]["evalue"]] for qid in only_blast
])
write_csv(sprot_all, ["SwissProt_query", "SwissProt_subject", "SwissProt_taxon", "SwissProt_bitscore", "SwissProt_evalue"], [
    [qid, swissprot_hits[qid]["subject"], swissprot_hits[qid]["taxon"], swissprot_hits[qid]["bitscore"], swissprot_hits[qid]["evalue"]] for qid in queries_sprot
])
write_fasta(sprot_faa, only_blast)

# --- STEP 8: Write NR Hits ---
nr_csv = f"{prefix}_nr.csv"
nr_faa = f"{prefix}_nr.faa"
nr_all = f"{prefix}_filtered_nr.csv"
write_csv(nr_csv, ["NR_query", "NR_subject", "NR_taxon", "NR_bitscore", "NR_evalue"], [
    [qid, nr_hits[qid]["subject"], nr_hits[qid]["taxon"], nr_hits[qid]["bitscore"], nr_hits[qid]["evalue"]] for qid in only_nr
])
write_csv(nr_all, ["NR_query", "NR_subject", "NR_taxon", "NR_bitscore", "evalue"], [
    [qid, nr_hits[qid]["subject"], nr_hits[qid]["taxon"], nr_hits[qid]["bitscore"], nr_hits[qid]["evalue"]] for qid in queries_nr
])
write_fasta(nr_faa, only_nr)

# --- STEP : Write Intersection ---
inter_csv = f"{prefix}_intersection.csv"
inter_faa = f"{prefix}_intersection.faa"
write_csv(inter_csv, [
    "query", "KO", "KEGG_taxon", "KEGG_score",
    "SwissProt_subject", "SwissProt_taxon", "SwissProt_bitscore", "SwissProt_evalue",
    "NR_subject", "NR_taxon", "NR_bitscore", "NR_evalue"], [
        [
            qid,
            kegg_hits[qid]["ko"],
            kegg_hits[qid]["taxon"],
            kegg_hits[qid]["score"],
            swissprot_hits[qid]["subject"],
            swissprot_hits[qid]["taxon"],
            swissprot_hits[qid]["bitscore"],
            swissprot_hits[qid]["evalue"],
            nr_hits[qid]["subject"],
            nr_hits[qid]["taxon"],
            nr_hits[qid]["bitscore"],
            nr_hits[qid]["evalue"]
        ] for qid in intersection
    ])
write_fasta(inter_faa, intersection)

# --- STEP : Write No Hits ---
nohit_txt = f"{prefix}_nohits.txt"
nohit_faa = f"{prefix}_nohits.faa"
with open(nohit_txt, "w") as f:
    check_noHits = len(no_hits)
    f.write("\n".join(sorted(no_hits)))
write_fasta(nohit_faa, no_hits)

# --- STEP 9 = Write log file ---
log_text = f"{'Command line: '} {'blastp-kegg_intersec.py'} {"--blast"} \
{swissprot_file} {'--kegg'} {kegg_file} {'--fasta'} {fasta_file} \
{'--out'} {prefix}\n\
\n\
{'=== Annotation Summary ==='}\n\
{'Total KEGG queries .......................: '} {check_kegg}\n\
{'KEGG queries (score ≥ threshold) .........: '} {check_filter_kegg}\n\
{'Queries only in KEGG .....................: '} {check_onlykegg}\n\
{'Total SwissProt queries ..................: '} {check_blast}\n\
{'SwissProt queries (with hit) .............: '} {check_filter_blast}\n\
{'Queries only in SwissProt ................: '} {check_onlyblast}\n\
{'Total NR queries .........................: '} {check_nr}\n\
{'NR queries (with hit) ....................: '} {check_filter_nr}\n\
{'Queries only in NR .......................: '} {check_onlynr}\n\
{'Queries in both (shared) .................: '} {check_intersec}\n\
{'Queries with no identification ...........: '} {check_noHits}\n\
{'=========================================== '}\n"
with open(log_file, "w", encoding="utf-8") as file:
    file.write(log_text)


print("All files written! :D")
