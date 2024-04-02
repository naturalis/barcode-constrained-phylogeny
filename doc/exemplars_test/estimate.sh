#!/usr/bin/env bash

# Activate environment
# mamba env create -f env.yml
# mamba activate exemplars_test

# Define locations of assets
DB=../databases/BOLD_COI-5P_barcodes.db
HMM=../../resources/hmm/COI-5P.hmm

# Iterate over tree files
for TREEFILE in $(ls *.tre); do
  BASE=${TREEFILE%.*}

  # Get sequences from database
  python seqs4tre.py -d $DB -i $TREEFILE > ${BASE}.fa

  # HMM align to Stockholm file
  hmmalign -o ${BASE}.sto --trim $HMM ${BASE}.fa

  # Convert Stockholm to FASTA
  python sto2fa.py ${BASE}.sto ${BASE}-aln.fa

  # Estimate the branch lengths
  raxml -m GTRGAMMA -n $BASE -s ${BASE}-aln.fa -f e -t $TREEFILE
done
