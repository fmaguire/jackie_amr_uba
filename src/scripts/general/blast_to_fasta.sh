#!/bin/bash
#
# input 1: blast result modified with indicies as the first feild and sequence
# as the second
#
# input 2: string to label (for example, prev for prevalence)
#

cat $1 | while read line; do
    id=$(echo $line | cut -d' ' -f1)
    seq=$(echo $line | cut -d ' ' -f2)
    echo -e ">lcl|$2|$id"
    echo -e "$seq" | sed 's/-//g'
done
