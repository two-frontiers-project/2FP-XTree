# full viral pipeline

# 1. identify synthetic metagenome sets for REP and for NONREP

python make_synthetic_metagenomes.py 50 100 250 1 1 pvc_nonrep_seqs.fa viral-nonrep-genomes
python make_synthetic_metagenomes.py 50 100 250 1 1 complete_high_quality.fa viral-ch-genomes
python make_synthetic_metagenomes.py 50 100 250 1 1 complete_high_medium_quality.fa viral-chm-genomes
python make_synthetic_metagenomes.py 50 100 250 1 1 complete_high_medium_low_quality.fa viral-chml-genomes

ls | grep genomeset | grep pvc_nonrep_seqs > viral_nonrep_config
ls | grep genomeset | grep complete_high_quality > viral_ch_config
ls | grep genomeset | grep complete_high_medium_quality > viral_chm_config
ls | grep genomeset | grep complete_high_medium_low_quality > viral_chml_config

# 2. generate actual metagenomes from sets made in step 1 for lowab and highab -- this will run xtree at k29

rm -rf /dev/shm/simb.xtr

while read p; do ./benchmark.sh $p  /mnt/f/braden2/complete_high_quality.fa.lin.fa.29.xtr 124 false true 1; done<viral_ch_config

mv xtree_output xtree_output_viral_29_ch

while read p; do ./benchmark.sh $p  /mnt/f/braden2/complete_high_medium_quality.fa.lin.fa.29.xtr 124 false true 1; done<viral_chm_config

mv xtree_output xtree_output_viral_29_chm

while read p; do ./benchmark.sh $p  /mnt/f/braden2/complete_high_medium_low_quality.fa.lin.fa.29.xtr 124 false true 1; done<viral_chml_config

mv xtree_output xtree_output_viral_29_chml

while read p; do ./benchmark.sh $p  /mnt/f/braden2/complete_high_medium_low_quality.fa.lin.fa.29.xtr 124 false true 1; done<viral_nonrep_config

mv xtree_output xtree_output_viral_29_nonrep_V2

# 3. generate config files for the k29 fastqs, run at k24 and k21 and k17

ls xtree_output_viral_29_ch/*R1*fastq.gz > r1s
ls xtree_output_viral_29_ch/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_viral_29_ch/*R1*fastq.gz | cut -f2 -d/ | cut -f1 -d. > names
paste names r1s r2s > xtree_output_viral_29_complete_high_CONFIG

rm -rf names r1s r2s

ls xtree_output_viral_29_chm/*R1*fastq.gz > r1s
ls xtree_output_viral_29_chm/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_viral_29_chm/*R1*fastq.gz | cut -f2 -d/ | cut -f1 -d. > names
paste names r1s r2s > xtree_output_viral_29_complete_high_medium_CONFIG

rm -rf names r1s r2s

ls xtree_output_viral_29_chml/*R1*fastq.gz > r1s
ls xtree_output_viral_29_chml/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_viral_29_chml/*R1*fastq.gz | cut -f2 -d/ |  cut -f1 -d.  > names
paste names r1s r2s > xtree_output_viral_29_complete_high_medium_low_CONFIG

rm -rf names r1s r2s

ls xtree_output_viral_29_nonrep_V2/*R1*fastq.gz > r1s
ls xtree_output_viral_29_nonrep_V2/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_viral_29_nonrep_V2/*R1*fastq.gz | cut -f2 -d/ |  cut -f1 -d.  > names
paste names r1s r2s > xtree_output_viral_29_nonrep_CONFIG

rm -rf names r1s r2s

for kmer in 17 21 24 ; do

for c in complete_high complete_high_medium complete_high_medium_low; do

rm -rf /dev/shm/simb.xtr

rsync "$c"_quality.fa.lin.fa."$kmer".xtr /dev/shm/simb.xtr

while read p; do

name=$(echo $p | cut -d ' ' -f1)
name2=$(echo $name | sed 's/_COV_0//g')
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
ab=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/ | cut -f5 -d_)
rep=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/ | cut -f6 -d_)
output=xtree_output_viral_"$kmer"_aligned_to_29_"$c"
mkdir -p $output
echo $output
echo $name

paste -dN <(zcat $r1) <(zcat $r2) | xtree --db /dev/shm/simb.xtr --seqs - --threads 124 --ref-out $output/$name.ref --cov-out $output/$name.cov --redistribute --doforage > "$output"/xaligns_"$name2".log

done<xtree_output_viral_29_"$c"_CONFIG

rm -rf /dev/shm/simb.xtr

done

done

# one more time for nonrep

for kmer in 17 21 24 ; do

rm -rf /dev/shm/simb.xtr

rsync complete_high_medium_low_quality.fa.lin.fa."$kmer".xtr /dev/shm/simb.xtr
echo 'it moved'
while read p; do

name=$(echo $p | cut -d ' ' -f1)
name2=$(echo $name | sed 's/_COV_0//g')
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
ab=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/ | cut -f5 -d_)
rep=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/ | cut -f6 -d_)
output=xtree_output_viral_"$kmer"_aligned_to_29_nonrep_V2
mkdir -p $output
echo $output
echo $name

paste -dN <(zcat $r1) <(zcat $r2) | xtree --db /dev/shm/simb.xtr --seqs - --threads 56 --ref-out $output/$name.ref --cov-out $output/$name.cov --redistribute --doforage > "$output"/xaligns_"$name2".log

done<xtree_output_viral_29_nonrep_CONFIG

rm -rf /dev/shm/simb.xtr

done

# 3. generate config files for the k29 fastqs, run at k24 and k21

rm -rf names r1s r2s

ls xtree_output_viral_29_nonrep_V2/*R1*fastq.gz > r1s
ls xtree_output_viral_29_nonrep_V2/*R1*fastq.gz | sed 's/.R1./.R2./g' > r2s
ls xtree_output_viral_29_nonrep_V2/*R1*fastq.gz | cut -f2 -d/ |  cut -f1 -d.  > names
paste names r1s r2s > xtree_output_vir_NONREP_CONFIG

# 4. run kraken2 on everything

conda activate

while read p; do

name=$(echo $p | cut -d ' ' -f1)
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
output=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/)
OUTPUT_PREFIX="$output"/"$name"
#mkdir -p $output

DB_PATH=/mnt/f/braden2/db/k2_vir

echo "kraken2 --confidence 0.1 --report-minimizer-data --db $DB_PATH --paired $r1 $r2 --output ${OUTPUT_PREFIX}_viralnonrep_kraken2_output.txt --report ${OUTPUT_PREFIX}_viralnonrep_kraken2_report.txt --threads 4 &> ${OUTPUT_PREFIX}_viralnonrep_kraken_bracken.log; bracken -d $DB_PATH -i ${OUTPUT_PREFIX}_viralnonrep_kraken2_report.txt -o ${OUTPUT_PREFIX}_bracken_output_Phylum.txt -r 150 -l P &>> ${OUTPUT_PREFIX}_viralnonrep_kraken_bracken.log; bracken -d $DB_PATH -i ${OUTPUT_PREFIX}_viralnonrep_kraken2_report.txt -o ${OUTPUT_PREFIX}_bracken_output_Class.txt -r 150 -l C &>> ${OUTPUT_PREFIX}_viralnonrep_kraken_bracken.log; bracken -d $DB_PATH -i ${OUTPUT_PREFIX}_viralnonrep_kraken2_report.txt -o ${OUTPUT_PREFIX}_bracken_output_Order.txt -r 150 -l O &>> ${OUTPUT_PREFIX}_viralnonrep_kraken_bracken.log; bracken -d $DB_PATH -i ${OUTPUT_PREFIX}_viralnonrep_kraken2_report.txt -o ${OUTPUT_PREFIX}_bracken_output_Family.txt -r 150 -l F &>> ${OUTPUT_PREFIX}_viralnonrep_kraken_bracken.log"

done<xtree_output_vir_NONREP_CONFIG > kraken_full_config

# 5. run metaphlan4 on everything

conda activate mpa

while read p; do

name=$(echo $p | cut -d ' ' -f1)
r1=$(echo $p | cut -d ' '  -f2)
r2=$(echo $p | cut -d ' '  -f3)
output=$(echo $p | cut -d ' '  -f2 | cut -f1 -d/)
#mkdir -p $output

metaphlan "$r1","$r2"  --input_type fastq --force --bowtie2out $output/"$name".bt2 --bowtie2db /mnt/f/braden2/db/metaphlan/ -o $output/"$name"_metaphlan.tsv --nproc 124 

done<all_synvir_metagenome_config

done

# 6. parse output 

Rscript bac_only_benchmark.R



