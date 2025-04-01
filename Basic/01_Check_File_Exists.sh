#!/bin/bash

#Checks if 'file.txt' exists

filename="file.txt"

if [ -f "$filename" ]; then
  echo "File '$filename' exists."
else
  echo "File '$filename' does not exist."
fi