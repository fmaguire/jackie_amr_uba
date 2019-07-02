#!/bin/bash
uba_id=''
id_count=1
count=1
column_labels=$(cat $1/UBA*txt | head -n 1)
#echo -e "local_id\tuba_grep_result\t"$column_labels
tail -n +1 $1/UBA*txt | grep -v 'ORF_ID' | sed '/^$/d' | grep "==>.*<==\|$2" | while read line; do
    if [[ $line == "==>"*"<==" ]]; then
	uba_id=$(echo "$line" | egrep -o 'UBA[0-9]*\.txt'| sed 's/.txt//g')
	count=1
     else
	echo -e "$id_count\t$uba_id-hit-$count\t$line"
	count=$((count+1))
	id_count=$((id_count+1))
    fi
   
done
