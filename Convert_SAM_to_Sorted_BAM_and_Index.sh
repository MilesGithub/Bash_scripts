#!/bin/bash

# Converts SAM to sorted, indexed BAM using samtools

# --- Configuration ---
input_sam="aligned_reads/sampleA.sam"
output_bam_prefix="aligned_reads/sampleA.sorted" # Output will be sampleA.sorted.bam
threads=4 # Threads for sorting

# --- Pre-run Checks ---
if [ ! -f "$input_sam" ]; then
    echo "Error: Input SAM file '$input_sam' not found!"
    exit 1
fi

output_bam="${output_bam_prefix}.bam"
echo "Converting '$input_sam' to sorted BAM '$output_bam'..."

# --- Processing Steps using Pipes ---
# 1. samtools view: Convert SAM to BAM (-b), include header (-h)
# 2. samtools sort: Sort the BAM data, use specified threads (-@), output to file (-o)
#    The '-' indicates that sort should read from standard input (the pipe)
samtools view -bh "$input_sam" | samtools sort -@ "$threads" -o "$output_bam" -

# Check if view | sort pipeline was successful
if [ $? -ne 0 ]; then
    echo "Error during SAM -> Sorted BAM conversion."
    exit 1
fi
echo "Sorted BAM file created: $output_bam"

# --- Indexing ---
echo "Indexing '$output_bam'..."
samtools index "$output_bam"

# Check if indexing was successful
if [ $? -ne 0 ]; then
    echo "Error during BAM indexing."
    exit 1
fi
echo "BAM index created: ${output_bam}.bai"

echo "Process completed successfully."
exit 0