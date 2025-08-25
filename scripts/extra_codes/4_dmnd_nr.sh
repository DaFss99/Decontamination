#!/bin/bash

# This is a script to run blastp for the reference database (database)to a list (id) of species 


braker_path="/home/dafne/inct_leveduras/wyll/braker3"
database="/home/dafne/database/protein_db/dmnd_reference_passalidarum.dmnd"
id=("y6407_B_scaffolds_500bp" "y2822_B_scaffolds_500bp" "y7005_B_scaffolds_500bp")


for item in "${id[@]}"; do
        diamond blastp \
                --query $braker_path/$item/braker.aa\
                --db "$database" \
                --out "$item"_passalidarum_Dblastp.txt \
                --outfmt 6 qseqid sseqid qstart qend sstart send qseq evalue bitscore stitle \
                --evalue 1e-5 \
                --threads 25 \
                --max-target-seqs 5

        echo "$item has finish! :D"

done

