from upsetplot import plot
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
import pandas as pd
import argparse
import textwrap
import csv
import re


# --- ARGUMENT PARSING ---
parser = argparse.ArgumentParser(
    description="Compare KEGG and SwissProt annotations and separate hits.")
parser.add_argument(
    "-s", "--swissprot", 
    required=True, 
    help="Path to BLASTP output file of nohits against SwissProt (outfmt 7)")
parser.add_argument(
    "-n", "--nr", 
    required=True, 
    help="Path to BLASTP output file of nohits against NR (outfmt 7)")
parser.add_argument(
    "-k", "--kegg", 
    required=True, 
    help="Path to KEGG .top result file")
parser.add_argument(
    "-ks", "--kegg_score", 
    required=True, 
    help="Chosed filtering kegg score.")
parser.add_argument(
    "-f", "--fasta", 
    required=True, 
    help="Path to nohits multifasta file (after blastp agains type organism proteome)")
parser.add_argument(
    "-o", "--out", 
    default="strain", 
    help="Output name prefix")
parser.add_argument(
    "--show_all_outputs", 
    action="store_true", 
    help="Show all intermediate outputs. Default is to hide them.")
args = parser.parse_args()

# --- FILE PATHS ---
log_file = Path(args.out).with_name(Path(args.out).stem + ".log")
sprot_file = Path(args.swissprot)   
nr_file = Path(args.nr)
kegg_file = Path(args.kegg)
kegg_filter = float(args.kegg_score)
fasta_file = Path(args.fasta)
prefix = Path(args.out).stem

# --- CHECKING VARIABLES ---
check_kegg = 0
check_sprot = 0
check_nr = 0
check_intersec = 0
check_onlykegg = 0
check_onlysprot = 0
check_onlynr = 0
check_noHits = 0
check_filter_kegg = 0
check_filter_sprot = 0
check_filter_nr = 0

# --- STEP 1: Parse KEGG ---
kegg_hits = {}
with kegg_file.open() as f:
    for line in f:
        # Thats a checkpoint to our log file
        check_kegg += 1
        parts = line.strip().split('\t')
        query = parts[0].replace("user:", "").strip()
        ko = parts[1].strip() if len(parts) > 1 else ""
        taxon = parts[2:6]
        score = float(parts[6]) if len(parts) > 6 else 0.0
        if score >= kegg_filter:
            # Thats a checkpoint to our log file
            check_filter_kegg += 1
            # Saving only the queries with scores bigger
            # than the choosed one
            kegg_hits[query] = {
                "ko": ko if ko else "None",
                "taxon": " | ".join(taxon),
                "score": score
            }

# --- STEP 2: Parse SwissProt (outfmt 7) ---
sprot_hits = {}
with sprot_file.open() as f:
    current_query = None
    for line in f:
        line = line.strip()
        if line.startswith("# Query:"):
            # Thats a checkpoint to our log file
            check_sprot += 1
            current_query = line.split()[-1]
        elif line.startswith("#") or not line:
            continue
        # For now on we are working with the queries with some hit
        elif current_query and current_query not in sprot_hits:
            # Thats a checkpoint to our log file
            check_filter_sprot += 1
            cols = line.split('\t')
            subject = cols[0]
            evalue = cols[6]
            bitscore = cols[7]
            subject_title = cols[8] if len(cols) > 8 else "NA"
            match = re.search(r"OS=([^=]+?) OX=", subject_title)
            taxon = match.group(1).strip() if match else "Unknown"
            sprot_hits[current_query] = {
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
                # Thats a checkpoint to our log file
                check_nr += 1
                current_query = line.split()[-1]
            elif line.startswith("#") or not line:
                continue
            # For now on we are working with the queries with some hit
            elif current_query and current_query not in nr_hits:
                # Thats a checkpoint to our log file
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
# Save the respective queries in a new 
queries_kegg = set(kegg_hits.keys())
queries_sprot = set(sprot_hits.keys())
queries_nr = set(nr_hits.keys())

intersection = queries_kegg & queries_sprot & queries_nr

only_kegg = queries_kegg - queries_sprot - queries_nr
only_sprot = queries_sprot - queries_kegg - queries_nr
only_nr = queries_nr - queries_kegg - queries_sprot

check_intersec = len(intersection)
check_onlysprot = len(only_sprot)
check_onlykegg = len(only_kegg)
check_onlynr = len(only_nr)

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
        SeqIO.write(
            (seq_dict[qid] for qid in ids if qid in seq_dict), 
            f, "fasta")

# --- STEP 6: Write KEGG Hits ---
kegg_csv = f"{prefix}_kegg.csv"
kegg_all = f"{prefix}_filtered_kegg.csv"
kegg_faa = f"{prefix}_kegg.faa"

if args.show_all_outputs:
    write_csv(kegg_csv, [
        "query", 
        "KO", 
        "taxon", 
        "score"], 
        [
        [qid, 
         kegg_hits[qid]["ko"], 
         kegg_hits[qid]["taxon"], 
         kegg_hits[qid]["score"]] 
         for qid in only_kegg
    ])
    write_csv(kegg_all, [
        "query", 
        "KO", 
        "taxon", 
        "score"
        ], [
        [qid, 
         kegg_hits[qid]["ko"], 
         kegg_hits[qid]["taxon"], 
         kegg_hits[qid]["score"]] 
         for qid in queries_kegg
    ])
    write_fasta(kegg_faa, only_kegg)

# --- STEP 7: Write SwissProt Hits ---
sprot_csv = f"{prefix}_sprot.csv"
sprot_faa = f"{prefix}_sprot.faa"
sprot_all = f"{prefix}_filtered_sprot.csv"

if args.show_all_outputs:
    write_csv(sprot_csv, [
        "SwissProt_query", 
        "SwissProt_subject", 
        "SwissProt_taxon", 
        "SwissProt_bitscore", 
        "SwissProt_evalue"
        ], [
        [qid, 
         sprot_hits[qid]["subject"], 
         sprot_hits[qid]["taxon"], 
         sprot_hits[qid]["bitscore"], 
         sprot_hits[qid]["evalue"]] 
         for qid in only_sprot
    ])
    write_csv(sprot_all, [
        "SwissProt_query", 
        "SwissProt_subject", 
        "SwissProt_taxon", 
        "SwissProt_bitscore", 
        "SwissProt_evalue"
        ], [
        [qid, 
         sprot_hits[qid]["subject"], 
         sprot_hits[qid]["taxon"], 
         sprot_hits[qid]["bitscore"], 
         sprot_hits[qid]["evalue"]] 
         for qid in queries_sprot
    ])
    write_fasta(sprot_faa, only_sprot)

# --- STEP 8: Write NR Hits ---
nr_csv = f"{prefix}_nr.csv"
nr_faa = f"{prefix}_nr.faa"
nr_all = f"{prefix}_filtered_nr.csv"

if args.show_all_outputs:
    write_csv(nr_csv, [
        "NR_query", 
        "NR_subject", 
        "NR_taxon", 
        "NR_bitscore", 
        "NR_evalue"
        ], [
        [qid, 
         nr_hits[qid]["subject"], 
         nr_hits[qid]["taxon"], 
         nr_hits[qid]["bitscore"], 
         nr_hits[qid]["evalue"]] 
         for qid in only_nr
    ])
    write_csv(nr_all, [
        "NR_query", 
        "NR_subject", 
        "NR_taxon", 
        "NR_bitscore", 
        "evalue"
        ], [
        [qid, 
         nr_hits[qid]["subject"], 
         nr_hits[qid]["taxon"], 
         nr_hits[qid]["bitscore"], 
         nr_hits[qid]["evalue"]] 
         for qid in queries_nr
    ])
    write_fasta(nr_faa, only_nr)

# --- STEP 9: Write Intersection ---
inter_csv = f"{prefix}_intersection.csv"
inter_faa = f"{prefix}_intersection.faa"
write_csv(inter_csv, [
    "query", "KO", "KEGG_taxon", 
    "KEGG_score", "SwissProt_subject", "SwissProt_taxon", 
    "SwissProt_bitscore", "SwissProt_evalue", "NR_subject", 
    "NR_taxon", "NR_bitscore", "NR_evalue"
    ], [
        [
            qid,
            kegg_hits[qid]["ko"],
            kegg_hits[qid]["taxon"],
            kegg_hits[qid]["score"],
            sprot_hits[qid]["subject"],
            sprot_hits[qid]["taxon"],
            sprot_hits[qid]["bitscore"],
            sprot_hits[qid]["evalue"],
            nr_hits[qid]["subject"],
            nr_hits[qid]["taxon"],
            nr_hits[qid]["bitscore"],
            nr_hits[qid]["evalue"]
        ] for qid in intersection
    ])
write_fasta(inter_faa, intersection)

# --- STEP 10: Write No Hits ---
nohit_txt = f"{prefix}_nohits.txt"
nohit_faa = f"{prefix}_nohits.faa"
with open(nohit_txt, "w") as f:
    check_noHits = len(no_hits)
    f.write("\n".join(sorted(no_hits)))
write_fasta(nohit_faa, no_hits)

# --- STEP 11: Write log file ---
log_text = f"{'Command line: '} {'blastp-kegg_intersec.py'} {"--blast"} \
{sprot_file} {'--kegg'} {kegg_file} {'--fasta'} {fasta_file} \
{'--out'} {prefix}\n\
\n\
{'=== Annotation Summary ==='}\n\
{'Total KEGG queries .......................: '} {check_kegg}\n\
{'KEGG queries (score ≥ threshold) .........: '} {check_filter_kegg}\n\
{'Queries only in KEGG .....................: '} {check_onlykegg}\n\
{'Total SwissProt queries ..................: '} {check_sprot}\n\
{'SwissProt queries (with hit) .............: '} {check_filter_sprot}\n\
{'Queries only in SwissProt ................: '} {check_onlysprot}\n\
{'Total NR queries .........................: '} {check_nr}\n\
{'NR queries (with hit) ....................: '} {check_filter_nr}\n\
{'Queries only in NR .......................: '} {check_onlynr}\n\
{'Queries shared queries ...................: '} {check_intersec}\n\
{'Queries with no identification ...........: '} {check_noHits}\n\
{'=========================================== '}\n\
\n\
\n\
\n\
{'=== List of queries only found against KEGG ==='}\n\
{queries_kegg}\n\
\n\
\n\
{'=== List of queries only found against SwissProt ==='}\n\
{queries_sprot}\n\
\n\
\n\
{'=== List of queries only found against NR ==='}\n\
{queries_nr}\n"

wrapped_text = "\n".join(textwrap.fill(line, width=150) for line in log_text.splitlines())
with open(log_file, "w", encoding="utf-8") as file:
    file.write(wrapped_text)

# --- STEP 12: UpsetPlot ---
set_names = ['KEGG', 'SwissProt', 'NR']
all_elems = queries_kegg.union(queries_sprot).union(queries_nr)
df = pd.DataFrame([
    [e in queries_kegg, e in queries_sprot, e in queries_nr] 
    for e in all_elems], 
    columns = set_names)
df_up = df.groupby(set_names).size()


fig = plt.figure(figsize=(12, 8))  # set figure size
plot(df_up, orientation='horizontal', show_counts=True)

# set and adjust title
plt.suptitle(
    "Hit distribution between KEGG, SwissProt and NR", 
    fontsize=16, 
    y=1.02)

# save figure
plt.savefig(f"{prefix}_upset.png", dpi=300, bbox_inches='tight')
plt.close()




print("All files written! :D")
