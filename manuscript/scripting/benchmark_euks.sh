#!/bin/bash

# --- Core Simulation Parameters ---
MIN_COVERAGE=0.01
MAX_COVERAGE=0.1
MIN_GENOMES=10
MAX_GENOMES=25
# Define the final log file name
FINAL_COVERAGE_LOG="synthetic_metagenome_coverage_log.tsv"

# --- Function Definitions ---

random_coverage() {
    echo "scale=3; $MIN_COVERAGE + ($MAX_COVERAGE - $MIN_COVERAGE) * $(awk -v seed=$RANDOM 'BEGIN { srand(seed); print rand() }')" | bc
}

calculate_min_coverage() {
    genome_length=$(grep -v ">" $1 | tr -d '\n' | wc -c)
    echo "scale=3; (151 / $genome_length)" | bc
}


# --- Core Function to Create a Single Metagenome (Safe for Parallel) ---
# Argument 1: Multi-FASTA file path
# Argument 2: Output base prefix (e.g., 'community')
# Argument 3: Unique Metagenome ID (e.g., '1', '2', '3')
# Argument 4: Threads for dwgsim (max available threads)

create_single_metagenome() {
    local multi_fasta_file=$1
    local output_base_prefix=$2
    local metagenome_id=$3
    local threads=$4
    
    # Define the unique temporary log file for this parallel job
    local TEMP_LOG_FILE="coverage_log_${metagenome_id}_$$.tmp"
    
    # Unique directory for this specific run, using the Metagenome ID and PID for isolation
    local TMP_DIR="temp_sim_${metagenome_id}_$$"
    mkdir -p "$TMP_DIR"
    
    local OUTPUT_PREFIX="${output_base_prefix}_${metagenome_id}"
    
    echo "--- [M$metagenome_id] Starting Metagenome: $OUTPUT_PREFIX ---"
    
    # --- Step 1: Select N Genomes ---
    local num_genomes=$(( $MIN_GENOMES + $RANDOM % ($MAX_GENOMES - $MIN_GENOMES + 1) ))

    # Get the list of all headers (passed to function for speed)
    local all_headers=$(grep '^>' "$multi_fasta_file")
    local selected_headers=$(echo "$all_headers" | shuf -n "$num_genomes")
    
    # --- Step 2: Split and Simulate (Isolating files by ID) ---
    
    # Use AWK to filter the multi-FASTA and split the selected genomes into the isolated TMP_DIR
    echo "$selected_headers" | awk -v multi_fasta="$multi_fasta_file" -v out_dir="$TMP_DIR" '
        BEGIN {
            while ((getline line < "-") > 0) {
                sub(/^>/, "", line);
                selected[line] = 1;
            }
            close("-");
            
            while ((getline line < multi_fasta) > 0) {
                if (line ~ /^>/) {
                    header = line;
                    sub(/^>/, "", header);
                    if (selected[header]) {
                        is_selected = 1;
                        filename = header;
                        gsub(/[ \/\\:;=,]/, "_", filename);
                        filename = out_dir "/" filename ".fasta";
                        print line > filename;
                    } else {
                        is_selected = 0;
                    }
                } else if (is_selected) {
                    print line >> filename;
                }
            }
        }
    ' -

    local individual_fasta_files=$(ls "$TMP_DIR"/*.fasta)
    local sim_commands=""
    
    for fasta_file in $individual_fasta_files; do
        
        local coverage=$(random_coverage)
        local min_cov=$(calculate_min_coverage "$fasta_file")

        if (( $(echo "$coverage < $min_cov" | bc -l) )); then
            coverage=$min_cov
        fi
        
        # Log the critical information to the unique temporary file
        local genome_base=$(basename "$fasta_file" .fasta)
        echo -e "${metagenome_id}\t${OUTPUT_PREFIX}\t${genome_base}\t${coverage}" >> "$TEMP_LOG_FILE"
        
        # Unique read prefix includes Metagenome ID and Genome name
        local dwgsim_read_prefix="${OUTPUT_PREFIX}_${genome_base}"
        
        sim_commands+=$"dwgsim -e 0.001 -E 0.001 -1 150 -2 150 -C $coverage -r 0 -X 0 -o 1 $fasta_file $dwgsim_read_prefix"$'\n'
    done
    
    # Execute all simulation commands for THIS single metagenome
    echo "$sim_commands" | parallel -j "$threads" --halt-on-error 1
    
    # --- Step 3: Concatenate Reads ---
    
    local FINAL_R1="${OUTPUT_PREFIX}_R1.fastq.gz"
    local FINAL_R2="${OUTPUT_PREFIX}_R2.fastq.gz"

    : > "$FINAL_R1"
    : > "$FINAL_R2"

    for fasta_file in $individual_fasta_files; do
        local genome_base=$(basename "$fasta_file" .fasta)
        local dwgsim_read_prefix="${OUTPUT_PREFIX}_${genome_base}"
        
        cat "${dwgsim_read_prefix}.bwa.read1.fastq.gz" >> "$FINAL_R1"
        cat "${dwgsim_read_prefix}.bwa.read2.fastq.gz" >> "$FINAL_R2"
        
        # Clean up the individual dwgsim files
        rm "${dwgsim_read_prefix}".bwa.read*.fastq.gz
        rm "${dwgsim_read_prefix}".mutations.* 2>/dev/null
    done
    
    # --- Step 4: Final Cleanup ---
    rm -rf "$TMP_DIR"
    
    echo "--- [M$metagenome_id] Complete. Logged to $TEMP_LOG_FILE. ---"
}

# Export the function and variables so GNU Parallel can access them
export -f create_single_metagenome
export MIN_COVERAGE MAX_COVERAGE MIN_GENOMES MAX_GENOMES FINAL_COVERAGE_LOG
export -f random_coverage calculate_min_coverage

# --- Main Script Execution ---

# 1. Input Multi-FASTA file
multi_fasta_file=$1
# 2. Output Metagenome Base Prefix (e.g., 'community')
output_base_prefix=$2
# 3. Number of synthetic metagenomes (N)
n_metagenomes=$3
# 4. Total available threads
total_threads=$4

# --- Validation ---
if [ -z "$multi_fasta_file" ] || [ -z "$output_base_prefix" ] || [ -z "$n_metagenomes" ] || [ -z "$total_threads" ]; then
  echo "USAGE: $0 <multi_fasta_file> <output_base_prefix> <N_metagenomes> <total_threads>"
  echo "Example: $0 fungal_lingenome_subsampled.fasta community 10 8"
  exit 1
fi

echo "--- Starting Parallel Metagenome Creation ---"
echo "Creating $n_metagenomes synthetic metagenomes."
echo "Final coverage log will be saved to: $FINAL_COVERAGE_LOG"
echo "----------------------------------------------"

# --- Run Parallel Jobs ---
seq 1 "$n_metagenomes" | parallel -j "$total_threads" "create_single_metagenome $multi_fasta_file $output_base_prefix {} $total_threads"

# --- Final Log Aggregation ---
echo "--- Aggregating Final Log File ---"

# 1. Write the header
echo -e "Metagenome_ID\tMetagenome_Prefix\tGenome_ID\tCoverage" > "$FINAL_COVERAGE_LOG"

# 2. Concatenate all temporary log files (safely)
cat coverage_log_*.tmp >> "$FINAL_COVERAGE_LOG"

# 3. Clean up the temporary log files
rm coverage_log_*.tmp

echo "SUCCESS! All $n_metagenomes metagenomes created and logged."
echo "Results and log saved in the current directory."
