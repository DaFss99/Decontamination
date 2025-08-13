#!/bin/bash
set -e # exit the run if there are any problem 

intersect_path="/home/dafne/inct_leveduras/wyll/cleaning/blastp/afterNoHit"
multifasta_path="/home/dafne/inct_leveduras/wyll/cleaning/blastp"
ks="50"
id=("y2822_500bp" "y6407_500bp" "y7005_500bp")
# id="/path/to/my/accession.txt"

# Instead of using the id list, you can also create a .txt file with all the 
# entries you want to run.
# This is a better option if you have more than 5 species.
# You can just edit the second $id and remove the activated one.

for item in "${id[@]}"; do

cd "$intersect_path"
if [ ! -d "$item" ]; then
        mkdir "$item"
fi

cd "$item"

python3 intersection.py \
-k "$intersect_path"/"$item"_kegg.out.top \
-s "$intersect_path"/"$item"_sprot_blastp.txt \
-ks "$ks" \
-n "$intersect_path"/"$item"_nr_blastp.txt \
-f "$multifasta_path"/"$item"_nohits.faa \
-o "$item"

done
