#!/bin/bash

# List of chromosomes to include
include_chromosomes="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY chrM"

# Get the unique chromosome names from the bed file and filter them
chromosomes=$(cut -f1 intermediate.coordinates.bed | sort | uniq | grep -w -E "chr[1-9]|chr1[0-9]|chrX|chrY|chrM")

# Loop through each chromosome name
for chromo in $chromosomes; do
  # Check if the chromosome is in the include list
  if echo "$include_chromosomes" | grep -w -q "$chromo"; then
    # Construct the pairwise alignment file path
    alignment_file="B6_FVBNJ_pairwisealignment_${chromo}.txt"

    # Run the Python script with the current chromosome name and alignment file
    python3 CoordinateMapping.py \
      -b intermediate.coordinates.bed \
      -p $alignment_file \
      -t g2 \
      -c $chromo \
      -a B6.anchor_info.csv FVBNJ.anchor_info.csv \
      -i 1 \
      -r 1 \
      -cpg 0 \
      -gn B6 FVBNJ \
      -o intermediate.coordinates.mapped.to.FVB.coordinates
  fi
done
