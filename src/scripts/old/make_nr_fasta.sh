#!/bin/bash
uba_id=''
count=1
cat $1 | while read line; do
    nr_id=$(echo -e "$line"| cut -f1)
    seq=$(echo -e "$line" | cut -f2 | sed 's/ //g'| sed 's/-//g')
    echo -e ">lcl|nrdb|$nr_id\n$seq"
done
