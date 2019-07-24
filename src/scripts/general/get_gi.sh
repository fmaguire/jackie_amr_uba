#!/bin/bash
#
#TODO check for installed esearch and eformat
gi_num=$1
esearch -db nuccore -query $gi_num | efetch -format fasta
