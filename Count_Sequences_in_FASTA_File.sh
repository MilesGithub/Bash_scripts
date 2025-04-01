#!/bin/bash

# Simple script to count sequences in a FASTA file

# Input: Specify the FASTA file
fasta_file="input_sequences.fasta"

# Check if the file exists
if [ ! -f "$fasta_file" ]; then
  echo "Error: File '$fasta_file' not found!"
  exit 1 # Exit with a non-zero status indicating an error
fi

# Use grep to find lines starting with '>' and count them (-c flag)
sequence_count=$(grep -c "^>" "$fasta_file")

echo "Number of sequences found in '$fasta_file': $sequence_count"

exit 0 # Exit successfully