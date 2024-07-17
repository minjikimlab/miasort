#!/bin/bash

# Stop on errors
set -Eeuo pipefail
set -x

# Define variables for paths and options
path1="GM12878-cohesin-pooled_comp_FDR_0.1_ALL_motifext4kbboth.region.PEanno"
path2="test-july-1.bedte"
type="abc"
graphs="BtoA;BtoC;AtoB;AtoC;CtoA;CtoB;AandC;Bcentered"
numfrag_min=2
anchor_options="no"
out_dir="test_folder"

# Check if the script is already running
if pgrep -f "python main.py"; then
    echo "The script is already running. Exiting."
    exit 1
fi

# Run the Python script in the background
python main.py --path1 "$path1" --path2 "$path2" --type "$type" --graphs "$graphs" --numfrag_min "$numfrag_min" --anchor_options "$anchor_options" --out_dir "$out_dir" &

# Get the PID of the Python script
script_pid=$!

# Initialize max memory usage
max_mem_usage=0

# Function to check memory usage
check_memory_usage() {
    mem_usage_kb=$(ps -o rss= -p "$script_pid")
    mem_usage_mb=$((mem_usage_kb / 1024))
    echo "Current Memory Usage: ${mem_usage_mb} MB"
    if (( mem_usage_mb > max_mem_usage )); then
        max_mem_usage=$mem_usage_mb
    fi
}

# Wait a bit to let the process start and gather some memory usage data
sleep 5

# Check memory usage periodically
while kill -0 "$script_pid" 2> /dev/null; do
    check_memory_usage
    sleep 10 # Adjust the sleep time as needed
done

# Wait for the script to finish
wait "$script_pid"

# Print the maximum memory usage
echo "Maximum Memory Usage: ${max_mem_usage} MB"
