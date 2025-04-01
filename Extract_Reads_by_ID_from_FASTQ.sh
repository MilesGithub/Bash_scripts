#!/bin/bash

# Extracts specific reads from a FASTQ file based on a list of IDs

# --- Configuration ---
id_list_file="read_ids_to_extract.txt" # File with one read ID per line
input_fastq="raw_reads/large_dataset.fastq.gz" # Input FASTQ (can be gzipped)
output_fastq="extracted_reads.fastq" # Output FASTQ file

# --- Pre-run Checks ---
if ! command -v seqtk &> /dev/null; then
    echo "Error: seqtk is not installed or not in PATH. Please install it."
    exit 1
fi
if [ ! -f "$id_list_file" ]; then echo "Error: ID list file '$id_list_file' not found."; exit 1; fi
if [ ! -f "$input_fastq" ]; then echo "Error: Input FASTQ file '$input_fastq' not found."; exit 1; fi

echo "Extracting reads listed in '$id_list_file' from '$input_fastq'..."

# Use seqtk subseq: takes name list and fastq file
# Handles gzipped input automatically
seqtk subseq "$input_fastq" "$id_list_file" > "$output_fastq"

if [ $? -ne 0 ]; then
    echo "Error during seqtk execution."
    exit 1
fi

# Optional: Count extracted reads
output_count=$(grep -c '^@' "$output_fastq") # Count lines starting with @ (FASTQ header)
echo "Extraction complete. Saved $output_count reads to '$output_fastq'."

exit 0