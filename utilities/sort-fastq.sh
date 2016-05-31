#!/bin/bash
# ARGS: inputfile outputfile min_length min_qual
cat $1 | make-linear-fastq.pl | filt-fastq.pl $3 $4 | sort --stable --reverse -n -k6.9,6 -k7.12,7 | split-linear-fastq.pl > $2

