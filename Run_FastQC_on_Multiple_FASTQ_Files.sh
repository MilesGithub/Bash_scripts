#!/bin/bash

# Runs FastQC on all FASTQ files in a given directory

# Configuration
fastq_directory="raw_reads/"      # Directory containing your FASTQ files
output_directory="fastqc_reports/" # Where to save FastQC output

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Check if input directory exists
if [ ! -d "$fastq_directory" ]; then
    echo "Error: Input directory '$fastq_directory' not found!"
    exit 1
fi

echo "Starting FastQC analysis..."

# Loop through all files ending with .fastq or .fastq.gz in the directory
# Using 'find' is often more robust than a simple 'for file in *.fastq*' loop
find "$fastq_directory" -maxdepth 1 \( -name "*.fastq" -o -name "*.fastq.gz" \) -print0 | while IFS= read -r -d $'\0' fastq_file; do
    echo " -> Processing: $fastq_file"
    # Run fastqc, specifying the output directory
    fastqc "$fastq_file" --outdir "$output_directory"

    # Basic check if FastQC ran without immediate error
    if [ $? -ne 0 ]; then
        echo "   Warning: FastQC returned an error for $fastq_file"
    fi
done

echo "FastQC processing complete. Reports are in '$output_directory'."
exit 0