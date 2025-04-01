#!/bin/bash

# Creates a timestamped backup of a given file

# Check if a filename was provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <filename_to_backup>"
  exit 1 # Exit with an error code
fi

original_file="$1"
timestamp=$(date +%Y%m%d_%H%M%S) # Get current date and time
backup_file="${original_file}_backup_${timestamp}"

# Check if the original file exists
if [ ! -f "$original_file" ]; then
  echo "Error: File '$original_file' not found!"
  exit 1
fi

echo "Backing up '$original_file' to '$backup_file'..."
cp "$original_file" "$backup_file"

# Check if the copy was successful
if [ $? -eq 0 ]; then
  echo "Backup successful!"
else
  echo "Backup failed!"
  exit 1
fi

exit 0 # Exit successfully