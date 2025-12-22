#!/bin/bash

random_coverage() {
    echo "scale=2; 0.01 + (10 - 0.01) * $(awk -v seed=$RANDOM 'BEGIN { srand(seed); print rand() }')" | bc
}

random_coverage_enforce_halfx_above() {
    echo "scale=2; 0.25 + (10 - 0.25) * $(awk -v seed=$RANDOM 'BEGIN { srand(seed); print rand() }')" | bc
}

calculate_min_coverage() {
    genome_length=$(grep -v ">" $1 | tr -d '\n' | wc -c)
    echo "scale=2; (151 / $genome_length)" | bc
}

# input folder of genomes to simulate from
genomes_to_simulate_from=$1

ghash=$(echo $genomes_to_simulate_from | sed 's/\///g')

# input location of xtree database
xtree_db=$2

# input number of threads
threads=$3

# input whether to actually do spike-in (true) or just background (false)
do_spikein=$4

# whether to use XTree 0.93+'s new "forage" mode
do_forage=$5

# an integer which, if supplied, will govern the number of backgrounds to simulate (default omitted)
num_backgrounds=$6

# set up coverages (all nonzero backgrounds must be present in spike-in coverages)
# Also, coverages include all coverages the background can simulate from (except 0)
coverages="0.01 0.02 0.05 0.125 0.25 0.5 1 2 4 8 16 32"
coveragesBG="0 0 0 0 0 0 0 0 0.01 0.02 0.05 0.125 0.25 0.5 1 2 4 8 16"
#coveragesX="0.12 0.25 0.35 0.5 0.75 1 2 4 8 12 16"

# Introduce the program
echo "Welcome to the benchmarking script for the xtree package. Version 0.2.0"
echo "This script will simulate mixed genome communities (+spike-in) for benchmarking."

# if input folder not provided, use current directory
if [ -z "$genomes_to_simulate_from" ]; then
  genomes_to_simulate_from="."
fi

# if input xtree database not provided, use /dev/shm/simdb.xtr
if [ -z "$xtree_db" ]; then
  xtree_db="/dev/shm/simdb.xtr"
  # later we will test if this exists; but proceed anyway until then
fi

# if input number of threads not provided, use all available
if [ -z "$threads" ]; then
  threads=$(nproc)
fi

# if input for spike-in not provided, use true
if [ -z "$do_spikein" ]; then
  do_spikein=true
  echo "Using spike-in mode -- to disable, set the fourth argument to false"
fi

# if we're not doing spike-in, just use the background coverages
if [ "$do_spikein" == "false" ]; then
    coverages=$coveragesBG
fi

# if input forage flag not provided, use true
if [ -z "$do_forage" ]; then
  do_forage=true
  echo "Using forage mode (REQUIRES XTree 0.93+) -- to disable, set the fifth argument to false"
fi

echo "Using xtree database at $xtree_db"
echo "Using $threads threads"
echo "Using spike-in mode: $do_spikein"
echo "Using forage mode: $do_forage"

if [ ! -z "$num_backgrounds" ]; then
  do_spikein=false
  genomes=$(ls $genomes_to_simulate_from/*.fasta | sed 's/.*\///')
  echo "Generating $num_backgrounds simulated backgrounds"
  echo "Using coverage levels ranging from 0.01x to 10x, sampled uniformly"
  backgrounds=$genomes
  echo -e "Genome\tCoverage\tnum_background" > Background."$ghash".log
  for i in $(seq 1 $num_backgrounds); do
    for genome in $backgrounds; do
        coverage=$(random_coverage)
        min_coverage=$(calculate_min_coverage $genomes_to_simulate_from/$genome)

        if [ "$coverage" == "0" ]; then
            continue
        fi

        if [ "$coverage" == "" ]; then
            echo "Random coverage generation failed, trying again"
            coverage=$(random_coverage)
            echo $coverage
        fi

        if (( $(echo "$coverage < $min_coverage" | bc -l) )); then
            echo "Coverage of "$coverage" is below "$min_coverage", setting to 0.25x or above..."
            #coverage=$(random_coverage_enforce_halfx_above)
            coverage=$min_coverage
            echo "dwgsim -e 0.001 -E 0.001 -1 150 -2 150 -C "$coverage" -r 0 -X 0 -o 1 "$genomes_to_simulate_from"/"$genome" TEST"
        fi

        echo -e "${genome%.*}\t$coverage\t$i" >> Background."$ghash".log

        # Run dwgsim command
        dwgsim -e 0.001 -E 0.001 -1 150 -2 150 -C $coverage -r 0 -X 0 -o 1 $genomes_to_simulate_from/$genome $genomes_to_simulate_from/${genome%.*}_${i}_TEMP 2> /dev/null &
        pid_dwgsim=$!
        wait $pid_dwgsim
        # Concatenate the simulated reads
        #cat $genomes_to_simulate_from/${genome%.*}_${i}_TEMP.bwa.read1.fastq.gz >> "${i}"_"$ghash"_Background.R1.fastq.gz
        #cat $genomes_to_simulate_from/${genome%.*}_${i}_TEMP.bwa.read2.fastq.gz >> "${i}"_"$ghash"_Background.R2.fastq.gz

    done
    echo "Merging individual genomes"
    for genome in $backgrounds; do
        cat $genomes_to_simulate_from/${genome%.*}_${i}_TEMP.bwa.read1.fastq.gz >> "${i}"_"$ghash"_Background.R1.fastq.gz 
        cat $genomes_to_simulate_from/${genome%.*}_${i}_TEMP.bwa.read2.fastq.gz >> "${i}"_"$ghash"_Background.R2.fastq.gz 
    done
  done
  genomes_to_simulate_from="SKIP"
fi

if [ "$genomes_to_simulate_from" != "SKIP" ]; then
    echo "Simulating from genomes in $genomes_to_simulate_from"
    echo "Coverage levels: $coverages"
    echo "Background coverages: $coveragesBG"

    # Step 1: Create an array of all genome files' base names
    genomes=$(ls $genomes_to_simulate_from/*.fasta | sed 's/.*\///')

    # Step 2: Loop over each genome file to simulate reads
    echo "Generating simulated reads for each genome at each coverage level"
    for genome in $genomes; do
        min_coverage=$(calculate_min_coverage $genomes_to_simulate_from/$genome)
        for coverage in $coverages; do
            # Ensure coverage is not below the minimum threshold
            if (( $(echo "$coverage < $min_coverage" | bc -l) )); then
                echo "Coverage below floor, setting to 1 read..."
                coverage=$min_coverage
            fi
            # Create a unique prefix for each simulation
            prefix="${genome%.*}_$coverage"
            # Run dwgsim command
            dwgsim -e 0.001 -E 0.001 -1 150 -2 150 -C $coverage -r 0 -X 0 -o 1 $genomes_to_simulate_from/$genome $genomes_to_simulate_from/${genome%.*}_${i}_TEMP 2> /dev/null &
            pid_dwgsim=$!
            wait $pid_dwgsim
        done
    done

    # Step 3: Create backgrounds and spike-ins
    echo "Creating backgrounds"
    echo -e "Genome\tCoverage\tToSpike" > Background."$ghash".log
    for omit in $genomes; do
        echo -e "${omit%.*}\t0\t${omit%.*}" >> Background."$ghash".log
        backgrounds=$(echo "$genomes" | tr ' ' '\n' | grep -v $omit | tr '\n' ' ')
        
        for genome in $backgrounds; do
            coverage=$(echo $coveragesBG | tr ' ' '\n' | shuf -n 1)
            echo -e "${genome%.*}\t$coverage\t${omit%.*}" >> Background."$ghash".log

            if [ "$coverage" == "0" ]; then
                continue
            fi
            cat "${genome%.*}_$coverage.bwa.read1.fastq.gz" >> "${omit%.*}_"$ghash"_Background.R1.fastq.gz"
            cat "${genome%.*}_$coverage.bwa.read2.fastq.gz" >> "${omit%.*}_"$ghash"_Background.R2.fastq.gz"
        done

        if [ "$do_spikein" == "true" ]; then
            for coverage in $coverages; do
                prefix="${omit%.*}_"$ghash"_Spikein_Cov_$coverage"
                cp "${omit%.*}_"$ghash"_$coverage.bwa.read1.fastq.gz" "$prefix."$ghash".R1.fastq.gz"
                cp "${omit%.*}_"$ghash"_$coverage.bwa.read2.fastq.gz" "$prefix."$ghash".R2.fastq.gz"
            done
        fi
    done

    rm *.mutations.*
    rm *.read[12].fastq.gz

else
    genomes=$(ls *"$ghash"*_Background.R1.fastq.gz | sed 's/_Background.R1.fastq.gz//')
    echo "Backgrounds: $(echo $genomes | tr '\n' ' ')"
fi

# Step 4: Align reads from each genome's spike-in to the reference metagenomic background
echo "Aligning reference backgrounds (and spike-ins, if applicable) to xtree database"

# Test if the xtree database exists and fail if it doesn't
if [ ! -f "$xtree_db" ]; then
  echo "ERROR: xtree database not found at $xtree_db"
    exit 1
fi

# Copy xtree database to /dev/shm/ for faster access, with dd
# But first, check if it's already there (and record if it isn't so we can delete it later)
copied_to_devshm=false
if [ ! -f "/dev/shm/simdb.xtr" ]; then
  echo "Copying xtree database to /dev/shm/ for faster access"
  numactl --interleave=all dd if=$xtree_db of=/dev/shm/simdb.xtr bs=1M
    # Record that we copied the database to /dev/shm/ so we can delete it later
    copied_to_devshm=true
    # otherwise warn the user that we're using the existing database
else
    echo "NOTE: Using existing xtree database at /dev/shm/simdb.xtr"
fi

# Make output directory if it doesn't exist
mkdir -p xtree_output

# An option for the user, if they want to add the new forage mode into the xtree invocation as "--doforage"
addforage=""
if [ "$do_forage" == "true" ]; then
    addforage="--doforage"
fi

commands=""
# The same list of genomes applies to both backgrounds and spike-ins
for genome in $genomes; do
  # Align spike-in to its dedicated background. Pipe the combined reads from both the spike-in and background to xtree
  if [ "$do_spikein" == "true" ]; then
    for coverage in $coverages; do
      # To be executed later in parallel
      commands+=$"echo ${genome%.*}_COV_$coverage ; "
      commands+=$"paste -dN <(zcat ${genome%.*}_Spikein_Cov_$coverage.R1.fastq.gz ${genome%.*}_Background.R1.fastq.gz) "
      commands+=$"<(zcat ${genome%.*}_Spikein_Cov_$coverage.R2.fastq.gz ${genome%.*}_Background.R2.fastq.gz) | "
      commands+=$"xtree --db /dev/shm/simdb.xtr --seqs - --threads "$threads" --ref-out xtree_output/${genome%.*}_COV_$coverage.ref "
      commands+=$"--cov-out xtree_output/${genome%.*}_COV_$coverage.cov --redistribute $addforage"$'\n'
    done
  fi
  # Align the background by itself (guaranteed omission of spike-in genome)

  commands+=$"echo ${genome%.*}_COV_0 ; "
  commands+=$"paste -dN <(zcat ${genome%.*}_Background.R1.fastq.gz) <(zcat ${genome%.*}_Background.R2.fastq.gz) | "
  commands+=$"xtree --db /dev/shm/simdb.xtr --seqs - --threads "$threads" --ref-out xtree_output/${genome%.*}_COV_0.ref "
  commands+=$" --cov-out xtree_output/${genome%.*}_COV_0.cov --redistribute $addforage"$'\n'
done

# Save the commands to a file
commands=$(echo "$commands" | sed 's/  */ /g')
echo "$commands" > xtree_commands_"$ghash".txt

# Run xtree commands in parallel (8 threads per command)
#threads=$((threads/8))
if [ $threads -lt 1 ]; then
  threads=1
fi
echo "$commands" | parallel -j $threads | grep -v "Processed .* queries" > xtree_output/xaligns_"${ghash}".log

mv  *Background*gz xtree_output
mv xtree_commands_genomeset_*txt xtree_output
mv Background.genomeset_*log xtree_output
