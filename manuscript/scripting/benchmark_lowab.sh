#!/bin/bash

random_coverage() {
    echo "scale=4; 0.0001 + (1 - 0.0001) * $(awk -v seed=$RANDOM 'BEGIN { srand(seed); print rand() }')" | bc
}

random_coverage_enforce_halfx_above() {
    echo "scale=4; "$1" + (1 - "$1") * $(awk -v seed=$RANDOM 'BEGIN { srand(seed); print rand() }')" | bc
}

calculate_min_coverage() {
    genome_length=$(grep -v ">" $1 | tr -d '\n' | wc -c)
    echo "scale=5; (151 / $genome_length)" | bc
}

# input folder of genomes to simulate from
genomes_to_simulate_from=$1

ghash=$(echo $genomes_to_simulate_from | sed 's/\///g')

# input location of xtree database
xtree_db=/mnt/f/braden2/gtdbr214_k29.xtr2

# input number of threads
threads=124
# whether to use XTree 0.93+'s new "forage" mode
do_forage=true

# an integer which, if supplied, will govern the number of backgrounds to simulate (default omitted)
num_backgrounds=1

echo "Kicking of synthetic metagenome generation..."

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


# if input forage flag not provided, use true
if [ -z "$do_forage" ]; then
  do_forage=true
  echo "Using forage mode (REQUIRES XTree 0.93+) -- to disable, set the fifth argument to false"
fi

echo "Using xtree database at $xtree_db"
echo "Using $threads threads"
echo "Using forage mode: $do_forage"

if [ ! -z "$num_backgrounds" ]; then
  genomes=$(ls $genomes_to_simulate_from/*.fasta | sed 's/.*\///')
  echo "Generating $num_backgrounds simulated backgrounds"
  echo "Using coverage levels ranging from 0.0001x to 1x, sampled uniformly"
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
            echo "Coverage of "$coverage" is below "$min_coverage", setting floor at min possible coverage of "$min_coverage" to attain one read..."
            coverage=$(random_coverage_enforce_halfx_above $min_coverage)
            #coverage=$min_coverage
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
fi

genomes=$(ls *"$ghash"*_Background.R1.fastq.gz | sed 's/_Background.R1.fastq.gz//')
echo "Backgrounds: $(echo $genomes | tr '\n' ' ')"

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

