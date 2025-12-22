# full bacterial pipeline

# 1. identify synthetic metagenome sets for REP and for NONREP

python make_synthetic_metagenomes.py 50 75 100 0 0 gtdb_r214_linearized.fa bacteria-gtdbrep
python make_synthetic_metagenomes.py 50 75 100 0 0 gtdb_nonreps_linearized.fna bacteria-gtdbnonrep

ls | grep genomeset | grep -v NONREP > rep_config
ls | grep genomeset | grep NONREP > nonrep_config

# 2. generate actual metagenomes from sets made in step 1 for lowab and highab -- this will run xtree at k29

rm -rf /dev/shm/simb.xtr

while read p; do ./benchmark.sh $p  /mnt/f/braden2/gtdbr214_k29.xtr2 124 false true 1; done<rep_config

mv xtree_output xtree_output_bac_29_highab_REP

### REDOING
while read p; do ./benchmark.sh $p  /mnt/f/braden2/gtdbr214_k29.xtr2 124 false true 1; done<nonrep_config

mv xtree_output xtree_output_bac_29_highab_NONREP

### REDOING

while read p; do ./benchmark_lowab.sh $p /mnt/f/braden2/gtdbr214_k29.xtr2 124 false true 1; done<rep_config

mv xtree_output xtree_output_bac_29_lowab_REP

while read p; do ./benchmark_lowab.sh $p /mnt/f/braden2/gtdbr214_k29.xtr2 124 false true 1; done<nonrep_config

mv xtree_output xtree_output_bac_29_lowab_NONREP

# 3. generate config files for the k29 fastqs, run at k24 and k21 and k17

ls xtree_output_bac_29_highab_REP/*R1*fastq.gz > r1s
ls xtree_output_bac_29_highab_REP/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_bac_29_highab_REP/*R1*fastq.gz | cut -f2 -d/ | cut -f1 -d. > names
paste names r1s r2s > xtree_output_bac_29_highab_REP_CONFIG

rm -rf names r1s r2s

ls xtree_output_bac_29_lowab_REP/*R1*fastq.gz > r1s
ls xtree_output_bac_29_lowab_REP/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_bac_29_lowab_REP/*R1*fastq.gz | cut -f2 -d/ |  cut -f1 -d.  > names
paste names r1s r2s > xtree_output_bac_29_lowab_REP_CONFIG

rm -rf names r1s r2s

ls xtree_output_bac_29_highab_NONREP/*R1*fastq.gz > r1s
ls xtree_output_bac_29_highab_NONREP/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_bac_29_highab_NONREP/*R1*fastq.gz | cut -f2 -d/ |  cut -f1 -d.  > names
paste names r1s r2s > xtree_output_bac_29_highab_NONREP_CONFIG

rm -rf names r1s r2s

ls xtree_output_bac_29_lowab_NONREP/*R1*fastq.gz > r1s
ls xtree_output_bac_29_lowab_NONREP/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_bac_29_lowab_NONREP/*R1*fastq.gz | cut -f2 -d/ |  cut -f1 -d.  > names
paste names r1s r2s > xtree_output_bac_29_lowab_NONREP_CONFIG

rm -rf names r1s r2s

cat *CONFIG > all_synbac_metagenome_config

for kmer in 17 ; do

rm -rf /dev/shm/simb.xtr

rsync gtdbr214_k"$kmer".xtr2 /dev/shm/simb.xtr

while read p; do

name=$(echo $p | cut -d ' ' -f1)
name2=$(echo $name | sed 's/_COV_0//g')
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
ab=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/ | cut -f5 -d_)
rep=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/ | cut -f6 -d_)
output=xtree_output_bac_"$kmer"_aligned_to_29_"$ab"_"$rep"
mkdir -p $output
echo $output
echo $name

paste -dN <(zcat $r1) <(zcat $r2) | xtree --db /dev/shm/simb.xtr --seqs - --threads 58 --ref-out $output/$name.ref --cov-out $output/$name.cov --redistribute --doforage > "$output"/xaligns_"$name2".log

done<all_synbac_metagenome_config

rm -rf /dev/shm/simb.xtr

done

# 4. run kraken2 on everything

conda activate

while read p; do

name=$(echo $p | cut -d ' ' -f1)
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
output=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/)
OUTPUT_PREFIX="$output"/"$name"
#mkdir -p $output

DB_PATH=/mnt/f/braden2/db/xtreer214/gtdb_r214_kraken2

echo "kraken2 --confidence 0.1 --report-minimizer-data --db $DB_PATH --paired $r1 $r2 --output ${OUTPUT_PREFIX}_kraken2_output.txt --report ${OUTPUT_PREFIX}_kraken2_report.txt --threads 4 &> ${OUTPUT_PREFIX}_kraken_bracken.log; bracken -d $DB_PATH -i ${OUTPUT_PREFIX}_kraken2_report.txt -o ${OUTPUT_PREFIX}_bracken_output_genus.txt -r 150 -l G &>> ${OUTPUT_PREFIX}_kraken_bracken.log; bracken -d $DB_PATH -i ${OUTPUT_PREFIX}_kraken2_report.txt -o ${OUTPUT_PREFIX}_bracken_output_species.txt -r 150 -l S &>> ${OUTPUT_PREFIX}_kraken_bracken.log"

done<tmpconfig > kraken_full_config

# 5. run metaphlan4 on everything

conda activate mpa

while read p; do

name=$(echo $p | cut -d ' ' -f1)
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
output=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/)
#mkdir -p $output

metaphlan "$r1","$r2"  --input_type fastq --force --bowtie2out $output/"$name".bt2 --bowtie2db /mnt/f/braden2/db/metaphlan/ -o $output/"$name"_metaphlan.tsv --nproc 124 

done<all_synbac_metagenome_config

done

# 6. parse output 

Rscript bac_only_benchmark.R



