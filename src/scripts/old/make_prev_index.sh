#!/bin/bash
uba_id=''
id_count=1
count=1
#column_labels=$(cat $1/UBA*txt | head -n 1)
echo -e "local_id\tuba_grep_result\t"$column_labels
cat $1 | while read header && read sequence; do
    aro=$(echo $header | cut -d '|' -f3)
    id=$(echo $header | cut -d '|' -f1)
    echo $id'\t'$aro
done
