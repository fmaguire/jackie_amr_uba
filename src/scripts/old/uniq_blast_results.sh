#!/bin/bash
grep_list=$(cat $1 | cut -f1 | sort | uniq | cut -d'|' -f3 | sed 's/^/\||/g' | sed 's/$/\$/g' | tr '\n' '\\|' | rev | cut -c 2- | rev | cut -c 2-)

cat $2 | grep -A1 $grep_list | sed 's/--//g' | sed '/^[[:space:]]*$/d'
#echo -e "'$grep_list'"
