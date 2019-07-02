#!/bin/bash

#index for the prev data

cat $1 | while read line && read sequence; do

    number=$(echo $line | cut -d '|' -f1 | cut -d ':' -f2)
    line2=$(echo -e "$line" | sed 's/>//g')
    echo -e "$number\t$line2"
 #    old_name=$(echo -e "$line" | cut -f1)
#    new_line=$(echo -e "$line" | sed 's/'$old_name'\t//g')
    #    new_line=
    
#    echo -e "$count\t$new_line\t$old_name"
done
