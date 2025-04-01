#!/bin/bash

# Basic BWA MEM alignment script for paired-end reads

# --- Configuration ---
reference_genome="reference/hg38.fasta"
read1_fastq="raw_reads/sampleA_R1.fastq.gz"
read2_fastq="raw_reads/sampleA_R2.fastq.gz"
output_sam="aligned_reads/sampleA.sam"
threads=8 # Number of CPU threads to use

# --- Pre-run Checks ---
# Ensure output directory exists
output_dir=$(dirname "$output_sam")
mkdir -p "$output_dir"

# Check if required files exist
if [ ! -f "$reference_genome" ]; then echo "Error: Reference genome '$reference_genome' not found."; exit 1; fi
if [ ! -f "$read1_fastq" ]; then echo "Error: Read 1 FASTQ '$read1_fastq' not found."; exit 1; fi
if [ ! -f "$read2_fastq" ]; then echo "Error: Read 2 FASTQ '$read2_fastq' not found."; exit 1; fi

# --- Step 1: Index Reference Genome (if needed) ---
# Check for a common index file (.bwt). If missing, run bwa index.
if [ ! -f "${reference_genome}.bwt" ]; then
    echo "Reference index not found. Indexing '$reference_genome'..."
    bwa index "$reference_genome"
    if [ $? -ne 0 ]; then echo "Error during BWA indexing. Exiting."; exit 1; fi
    echo "Indexing complete."
else
    echo "Reference index found."
fi

# --- Step 2: Alignment ---
echo "Starting BWA MEM alignment..."
echo "Reference: $reference_genome"
echo "Read 1: $read1_fastq"
echo "Read 2: $read2_fastq"
echo "Output: $output_sam"
echo "Threads: $threads"

# Run BWA MEM -t threads reference read1 read2 > output.sam
bwa mem -t "$threads" "$reference_genome" "$read1_fastq" "$read2_fastq" > "$output_sam"

# --- Post-run Check ---
if [ $? -ne 0 ]; then
    echo "Error during BWA MEM alignment. Check BWA output/logs."
    # Optional: remove potentially incomplete SAM file
    # rm -f "$output_sam"
    exit 1
fi

echo "Alignment finished successfully. Output saved to '$output_sam'."
exit 0