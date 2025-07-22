#!/bin/bash

# This is a script to run blastp to a list (id) of species 


braker_path="/home/dafne/inct_leveduras/wyll/filogenia/tree/braker3"
database="/home/dafne/database/protein_db/reference_passalidarum_index/reference_protein_passalidarum"
id=("y6407_B_-_scaffolds" "y2822_B_-_scaffolds" "y7005_B_-_scaffolds")


for item in "${id[@]}"; do

blastp \
-query $braker_path/$item/braker.aa \
-db $database -out "$item"_blastx.txt \
-outfmt "7 sseqid ssac qstart qend sstart send qseq evalue bitscore" \
-evalue 1e-5 \
-num_threads 30

echo "$item has finish! :D"

done
