#!/bin/bash

while read fastafile;

do

while read kmers;

do

echo sbatch -p panda -c 24 -n 1 --mem=96G -t 4-11:11 ./buildtree.sh $fastafile $kmers

./xtree BUILD --seqs $fastafile --log-out "$fastafile".xlog --db-out $fastafile.xtr --comp 1 --k $kmers --threads 24

done<kmerlist

done<fastalist


