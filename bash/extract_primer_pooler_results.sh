#!/bin/bash

# Extract primer pooler results

sed -n '1~8p' ../data/derived/primer_pooler_4thresh.txt > ../data/derived/primer_pooler_4thres_primers.txt

echo "Extracted primers below threshold"