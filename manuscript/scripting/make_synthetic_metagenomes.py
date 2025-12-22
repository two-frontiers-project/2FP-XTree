import sys
import os
from random import sample
import hashlib
import random

def read_fasta(filename, filetype):
    print(f"Loading fasta file: {filename}...")
    genomes = []
    with open(filename, 'r') as f:
        genome = []
        identifier = ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if genome:
                    genomes.append((identifier, "".join(genome)))
                    genome = []
                identifier = line.split()[0][1:]
            else:
                genome.append(line)
        if genome:
            genomes.append((identifier, "".join(genome)))
    print(f"Finished loading {filename}. Total genomes: {len(genomes)}")
    return [(filetype, identifier, genome) for identifier, genome in genomes]

def generate_file_name(viral_only):
    """Generates a unique file name using a hash."""
    prefix = "viralonly_" if viral_only else ""
    random_str = ''.join(random.choices('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', k=5))
    hashed_str = hashlib.md5(random_str.encode()).hexdigest()[:10]  # take only the first 10 characters for brevity
    return prefix + hashed_str + ".fasta"

def write_fasta(filetype, identifier, genome, filename):
    with open(filename, 'w') as f:
        f.write(">" + identifier + "\n")  # write the identifier with '>' prefix
        f.write(genome + "\n")  # write the genome sequence

def generate_directory_name(viral_only=False):
    """Generates a unique directory name using a hash."""
    random_str = 'genomes_' + ''.join(random.choices('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789', k=10))
    hashed_str = hashlib.md5(random_str.encode()).hexdigest()
    prefix = "genomeset_viralonly_" if viral_only else "genomeset_"
    return prefix + hashed_str


def main(genomedirnumber, genomemin, genomemax, viralproportionmin, viralproportionmax):
    # Load genomes
    genomes1 = read_fasta("complete_high_medium_low_quality.fa", "viral")
    genomes2 = []
    if not (viralproportionmin == 1 and viralproportionmax == 1):
        genomes2 = read_fasta("GTDB_phage_contigs_removed.fasta", "bacterial")
        viral_only = False
    else:
        viral_only = True

    metadata = []
    for _ in range(genomedirnumber):
        # Randomly select total genomes and proportions
        totalgenomes = random.randint(genomemin, genomemax)
        proportion = random.uniform(viralproportionmin, viralproportionmax)

        out_dir = generate_directory_name(viral_only)
        os.makedirs(out_dir)

        n_genomes1 = int(totalgenomes * proportion)
        n_genomes2 = totalgenomes - n_genomes1

        # Subsample genomes
        subsampled_genomes1 = sample(genomes1, n_genomes1)
        subsampled_genomes2 = sample(genomes2, n_genomes2) if genomes2 else []

        # Write genomes to separate files in the specified directory
        for _, (filetype, identifier, genome) in enumerate(subsampled_genomes1):
            file_name = generate_file_name(viral_only)
            write_fasta(filetype, identifier, genome, os.path.join(out_dir, file_name))

        if genomes2:  # Check if there are bacterial genomes
            for _, (filetype, identifier, genome) in enumerate(subsampled_genomes2):
                file_name = generate_file_name(False)  # Never "viralonly" for bacterial genomes
                write_fasta(filetype, identifier, genome, os.path.join(out_dir, file_name))


        metadata.append((out_dir, totalgenomes, proportion))

    # Write metadata to TSV file
    with open("genomes_metadata.tsv", "w") as f:
        f.write("Directory\tTotalGenomes\tProportionViral\n")
        for dir_name, genomes, prop in metadata:
            f.write(f"{dir_name}\t{genomes}\t{prop}\n")

    print("Script completed successfully!")

if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: script_name.py genomedirnumber genomemin genomemax viralproportionmin viralproportionmax")
        sys.exit(1)

    genomedirnumber = int(sys.argv[1])
    genomemin = int(sys.argv[2])
    genomemax = int(sys.argv[3])
    viralproportionmin = float(sys.argv[4])
    viralproportionmax = float(sys.argv[5])

    main(genomedirnumber, genomemin, genomemax, viralproportionmin, viralproportionmax)





    