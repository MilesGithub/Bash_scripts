#!/bin/bash

# Counts lines in all *.log files in the current directory

shopt -s nullglob # Important: Ensures loop doesn't run if no .log files exist

echo "Processing log files..."

for logfile in *.log; do
  line_count=$(wc -l < "$logfile") # Use input redirection for wc
  # Remove leading whitespace from wc output
  line_count=$(echo $line_count | xargs) 
  echo " -> File '$logfile' has $line_count lines."
done

shopt -u nullglob # Turn off nullglob if desired (good practice)

echo "Finished processing."