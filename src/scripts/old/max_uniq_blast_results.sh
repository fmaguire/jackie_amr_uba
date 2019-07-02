#!/bin/bash

#take the max size of the blast results for each name to retain max info
#then remove gaps for future alignment

#name list (there are multiple subject qury per name)
grep_list=$(cat $1 | cut -f1 | sort | uniq)
IFS=$'\n'
for line in $grep_list; do
    results_for_head=$(cat $1 | grep $line)
    max_seq_size=0
    max_seq=''
    for result in $results_for_head; do
	#echo -e "$result"
	seq=$(echo -e "$result" | cut -f2 | sed 's/-//g') # also clean gaps
	seq_size=$(echo $seq | wc -c)
	if (( $seq_size > $max_seq_size )); then
	    max_seq_size=$seq_size
	    max_seq=$seq
	fi
    done
    echo -e ">$line"
    echo -e "$max_seq"
done
