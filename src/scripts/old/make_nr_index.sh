#!/bin/bash

mkdir -p out
#echo $filecont
count=1

fasta=$(./max_uniq_blast_results.sh $1)

cat $1 | while read line; do
  #echo "blasting $database"
    taxid=$(echo -e "$line" | cut -f1)
#    echo $taxid
    line2=$(echo -e "$line" | sed 's/'$taxid'/ /g')
    echo -e "$count\t$line2\t$taxid"
    ((count++))
done
