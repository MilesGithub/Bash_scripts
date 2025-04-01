#!/bin/bash

# Reads a file named 'input.txt' line by line

input_file="input.txt"

if [ ! -f "$input_file" ]; then
    echo "Error: Input file '$input_file' not found."
    exit 1
fi

echo "Reading file '$input_file':"
line_num=0
while IFS= read -r line; do
    line_num=$((line_num + 1))
    echo "Line $line_num: $line"
    # You could add more processing here, e.g., grep for patterns
done < "$input_file" # Redirect the file content into the while loop

echo "Finished reading '$input_file'."