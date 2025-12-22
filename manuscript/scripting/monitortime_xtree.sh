#!/bin/bash

# Validate input files exist
if [[ ! -f "TRIAL_SAMPLES_RUNTIME" || ! -f "dblocs" ]]; then
    echo "ERROR: TRIAL_SAMPLES_RUNTIME or dblocs file not found!"
    exit 1
fi

# Define output CSV file
output_csv="parallel_runtime_results.csv"

# Initialize CSV with headers
echo "database,user_time_sec,system_time_sec,cpu_percent,elapsed_time,max_memory_kb" > "$output_csv"

# Read each database from `dblocs`
while IFS= read -r db; do
    # Extract clean database name for logging
    db_name=$(basename "$db" | sed 's/[^a-zA-Z0-9]/_/g')

    # Create a temporary `foo` file
    foo_file=foobar

    # Read each sample name from TRIAL_SAMPLES_RUNTIME and populate the foo file
    while IFS= read -r name; do
        echo "zcat ${name}*.fastq.gz | xtree --seqs - --threads 12 --db $db --ref-out prof/${name}.ref --cov-out prof/${name}.cov --redistribute --doforage" >> "$foo_file"
    done < TRIAL_SAMPLES_RUNTIME

    echo "Running parallel for DB: $db_name"

    # Run parallel and collect timing/memory stats
    /usr/bin/time -v bash -c "cat $foo_file | parallel -j3" 2>&1 | tee parallel_memory_runtime.log

    # Extract relevant stats
    user_time=$(awk -F': ' '/User time \(seconds\)/ {print $2}' parallel_memory_runtime.log)
    system_time=$(awk -F': ' '/System time \(seconds\)/ {print $2}' parallel_memory_runtime.log)
    cpu_percent=$(awk -F': ' '/Percent of CPU this job got/ {print $2}' parallel_memory_runtime.log)
    elapsed_time=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' parallel_memory_runtime.log)
    max_memory=$(awk -F': ' '/Maximum resident set size/ {print $2}' parallel_memory_runtime.log)

    # Append results to CSV
    echo "$db,$user_time,$system_time,$cpu_percent,$elapsed_time,$max_memory" >> "$output_csv"

    # Clean up temp file
    rm -f "$foo_file"

done < dblocs

echo "All runs completed. Results saved in $output_csv"
