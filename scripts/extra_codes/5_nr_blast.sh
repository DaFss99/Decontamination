#!/bin/bash

# This is a script to run blastp to a list (id) of species 


protein_path="/home/dafne/inct_leveduras/wyll/cleaning/blastp/"
database="/home/dafne/database/nr_08-2025/nr"
id=("y2822")


for item in "${id[@]}"; do

	blastp \
		-query $protein_path/"$item"_nohits.faa\
		-db $database \
		-out "$item"_nr_blastp.txt \
		-outfmt "7 sseqid ssac qstart qend sstart send qseq evalue bitscore stitle" \
		-evalue 1e-5 \
		-num_threads 15 \
		-max_target_seqs 5

	echo "$item has finish! :D"

done

começou 07/08/2025 às 13:51
terminou 08/08/2025 às 10:37
tempo da run com 698 proteínas: 1 dia, 20 horas e 46 minutos.