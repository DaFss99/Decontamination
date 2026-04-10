import argparse 
import pandas as pd
from pathlib import Path

# --- Argument parsing ---
parser = argparse.ArgumentParser(
    description= "Run score_kegg.py so you can easily get a Kegg score cutoff" \
    "based on a given e-value." \
    "Requested files: kegg.top and blast.txt result" \
    "Output files:" \
    "dataframe.txt with your kegg and blast scores" \
    "score.txt file containing a sugested K-numer cutoff"
)
parser.add_argument(
    "-k", "--kegg",
    required=True,
    help="Path to KEGG .top result file (e.g.:path/to/my_kegg.top)"
)
parser.add_argument(
    "-b", "--blast",
    required=True,
    help="Path to Blast .top result file (e.g.:path/to/my_blast.txt)"
    )
parser.add_argument(
    "-o", "--output",
    required=True,
    help="Name tag for output files (e.g.: strain_name, strain_code, test_id)"
)

args = parser.parse_args() # This line saves the parsed argumnets given as input

# --- File paths ---
kegg = Path(args.kegg)
blast = Path(args.blast)
prefix = Path(args.output)
output_file = Path(args.output).with_name(Path(args.output).stem + "score.txt")
log_file =  Path(args.output).with_name(Path(args.output).stem + "log.txt")


# --- Getting the needed data and storing it in a dictionary ---
kegg_dic = {}   # Creating a dictionary to save my query + it's related score
with kegg.open() as f:
    for line in f:
        parts = line.strip().split('\t')    # Security measure to standardize the archive
        query = parts[0].replace('user:', '').strip() # Use .strip() to remore whitespaces
        knumber = float(parts[6])
        kegg_dic[query] = knumber

blast_dic = {}
with blast.open() as f:
    for line in f:
        parts = line.strip().split('\t')
        query = parts[0]
        evalue = float(parts[7])
        if query not in blast_dic:
            blast_dic[query] = evalue
        elif evalue < blast_dic[query]:
            blast_dic[query] = evalue


# --- Creating a dataframe ---
## My dictionary as dataframe
kegg_df = pd.DataFrame.from_dict(kegg_dic, orient="index", columns=["kegg"])
blast_df = pd.DataFrame.from_dict(blast_dic, orient="index", columns=["blast"])

## Creating final dataframe
df = kegg_df.join(blast_df, how="inner")    # inner: only the queries in the intersection between files.

# --- Getting my KEGG score ---
df_filtrado = df[df["blast"] <= 1e-10]
valores_kegg = df_filtrado["kegg"]
media = valores_kegg.mean()
mediana = valores_kegg.median()

# print(valores_kegg)
print("Média:", media)
print("Mediana:", mediana)
print("Percentis do KEGG score para E-value <= 1e-10:")
print(valores_kegg.quantile([0.1, 0.25, 0.5, 0.75, 0.9]))

# --- Write outputs ---
## Saving my dataframe in a readably file
df.index.name = "Query"
df.to_csv(output_file, sep="\t", index=True)

## Log file
log_text = f"\n\
{'Command line:'}\n\
{'score_kegg.py'} {'--kegg'} {kegg} {'--blast'} {blast} {'--output'} {prefix} \
\n\
\n\
{'=== KEGG Metrics based on {blast} e-values ==='}\n\
{'Mean K-score: '} {media}\n\
{'Median K-score: '} {mediana}\n\
\n\
{'=== Sugested K-score cut-off ==='}\n\
            {mediana}\n"

with open(log_file, "w") as file:
    file.write(log_text)
