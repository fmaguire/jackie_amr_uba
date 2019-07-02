#!/bin/bash

filecont=$(ls $1/nr* | grep -o nr.[0-9].[0-9]* | uniq)
mkdir -p out
#echo $filecont
echo -e "$filecont" | while read database; do
  #echo "blasting $database"
  blastp -query $2 -db $1/$database -outfmt '6 staxids sseq evalue bitscore' -qcov_hsp_perc 90 -evalue 1e-180
done 
