#!/bin/bash

filecont=$(ls $1/nr* | grep -o nr.[0-9].[0-9]* | uniq)
#echo $filecont
echo -e "$filecont" | while read database; do
  blastp -query $2 -db $1/$database -outfmt '6 staxids sseq sallgi salltitles evalue bitscore' -qcov_hsp_perc 90 -evalue 1e-180
done 
