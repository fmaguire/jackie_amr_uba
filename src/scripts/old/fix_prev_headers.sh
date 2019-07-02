#!/bin/bash
cat $1 | while read header && read sequence; do
    lcllbl=$(echo $header | cut -d '|' -f1 | cut -d':' -f2)
    echo -e ">lcl|prev|$lcllbl\n$sequence"
done
