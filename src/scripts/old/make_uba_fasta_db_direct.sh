#!/bin/bash
uba_id=''
#echo -e "local_id\tuba_grep_result\t"$column_labels
tail -n +1 $1/UBA*txt |grep -v 'ORF_ID'|  sed '/^$/d' | grep "==>.*<==\|$2" | while read line; do
    if [[ $line == "==>"*"<==" ]]; then
	uba_id=$(echo "$line" | egrep -o 'UBA[0-9]*\.txt'| sed 's/.txt//g')
	count=1
    else
	seq=$(echo -e "$line" | cut -f19)
	line1=$(echo -e "$line" | cut --complement -f 18,19,20)
	line2=$(echo -e "$line1" | sed 's/\t/|/g')
	echo -e ">$uba_id|$line2"
	echo -e "$seq"
	count=$((count+1))
	id_count=$((id_count+1))
    fi
   
done
