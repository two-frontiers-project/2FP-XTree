#!/bin/bash

# Validate input files exist
if [[ ! -f "TRIAL_SAMPLES_RUNTIME" || ! -f "dblocs" ]]; then
    echo "ERROR: TRIAL_SAMPLES_RUNTIME or dblocs file not found!"
    exit 1
fi

# Define output CSV file
output_csv="parallel_runtime_results_xtree.csv"

# Initialize CSV with headers
echo "database,num_processes,command,user_time_sec,system_time_sec,cpu_percent,elapsed_time,max_memory_kb" > "$output_csv"

# Loop through each database
while IFS= read -r db; do
    # Extract clean database name for logging
    db_name=$(basename "$db" | sed 's/[^a-zA-Z0-9]/_/g')

    # Loop through different `j` values (1, 2, 3)
    for j in {1..3}; do
        # Create a temporary `foo` file
        foo_file=$(mktemp)

        # Read each sample name from TRIAL_SAMPLES_RUNTIME and populate the foo file
        while IFS= read -r name; do
            echo "zcat ${name}*.fastq.gz | xtree --seqs - --threads 12 --db $db --ref-out prof/${name}.ref --cov-out prof/${name}.cov --redistribute --doforage" >> "$foo_file"
        done < TRIAL_SAMPLES_RUNTIME

        echo "Running parallel -j$j for DB: $db_name"

        # Run parallel with the current `j` value and collect timing/memory stats
        /usr/bin/time -v bash -c "cat $foo_file | parallel -j$j" 2>&1 | tee parallel_memory_runtime.log

        # Extract relevant stats
        user_time=$(awk -F': ' '/User time \(seconds\)/ {print $2}' parallel_memory_runtime.log)
        system_time=$(awk -F': ' '/System time \(seconds\)/ {print $2}' parallel_memory_runtime.log)
        cpu_percent=$(awk -F': ' '/Percent of CPU this job got/ {print $2}' parallel_memory_runtime.log)
        elapsed_time=$(awk -F': ' '/Elapsed \(wall clock\) time/ {print $2}' parallel_memory_runtime.log)
        max_memory=$(awk -F': ' '/Maximum resident set size/ {print $2}' parallel_memory_runtime.log)

        # Append results to CSV
        echo "$db,$j,$(cat $foo_file),$user_time,$system_time,$cpu_percent,$elapsed_time,$max_memory" >> "$output_csv"

        # Clean up temp file
        rm -f "$foo_file"

    done
done < dblocs

echo "All runs completed. Results saved in $output_csv"
