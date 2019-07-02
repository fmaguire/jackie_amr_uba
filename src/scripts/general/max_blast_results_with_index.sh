#!/bin/bash

#take the max size of the blast results for each name to retain max info
#then remove gaps for future alignment

#name list (there are multiple subject qury per name)
grep_list=$(cat $1 | cut -f1 | sort | uniq)
IFS=$'\n'
count=1
for line in $grep_list; do
    results_for_head=$(cat $1 | grep $line)
    max_seq_size=0
    max_result=''
    
    for result in $results_for_head; do
	seq=$(echo -e "$result" | cut -f2 | sed 's/-//g') # also clean gaps
	seq_size=$(echo $seq | wc -c)
	if (( $seq_size > $max_seq_size )); then
	    max_seq_size=$seq_size
	    max_result=$result
	fi
    done
    taxid=$(echo -e "$max_result" | cut -f1)
#    echo $taxid
    line2=$(echo -e "$max_result" | sed 's/'$taxid'\t/ /g')
    echo -e "$count\t$line2\t$taxid"

    ((count++))
done
