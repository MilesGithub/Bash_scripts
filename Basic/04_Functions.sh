#!/bin/bash

# Demonstrates using a function to print system info

# Define the function
print_system_info() {
  echo "--- System Information ---"
  echo "Hostname: $(hostname)"
  echo "Kernel: $(uname -r)"
  echo "Uptime: $(uptime -p)"
  echo "--------------------------"
}

# Call the function
echo "Gathering system info..."
print_system_info
echo "Done."