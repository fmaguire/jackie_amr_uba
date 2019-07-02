#!/bin/bash
uba_id=''
count=1
tail -n +2 $1 | while read line; do
    uba_id=$(echo -e "$line"| cut -f1)
    seq=$(echo -e "$line" | cut -f21)
    echo -e ">lcl|uba|$uba_id\n$seq"
done
