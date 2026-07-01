#!/bin/bash

# Convert primer_pooler.csv to FASTA format for PrimerPooler

sed 's/^/>/' ../data/derived/primer_pooler_filt.csv \
| tr ',' '\n' \
| sed '1,2d' \
| sed 's/"//g' \
> primer_pooler_filt.fasta

# Append GT-seq tag sequences
cat <<EOF >> primer_pooler_filt.fasta
>tagF
ACACTCTTTCCCTACACGACGCTCTTCCGATCT
>tagR
GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
EOF

echo "Created primer_pooler_filt.fasta"