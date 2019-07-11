#!/bin/bash
#
#
# create headers for a fasta file with indicies which correspond to metadata
# generated in another file
#
# input 1: fasta file path (note: this file should have a corresponding index
# file)
#
# input 2: string to label the type of sequence (example: canon for canonical)
#

count=1
cat $1 | while read header && read sequence; do
    echo -e ">lcl|$2|$count\n$sequence"
    ((count++))
done
