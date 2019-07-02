#!/bin/bash
#
# create a tabular file with an index as a first field and the header
# information from the fasta in a second feild
#
# input 1: fasta file path (note: this file should have a corresponding fasta
# file with indicies)
#

count=1
cat $1 | while read header && read sequence; do
    info=$(echo $header | sed 's/>//g')
    echo -e "$count\t$info"
    ((count++))
done
