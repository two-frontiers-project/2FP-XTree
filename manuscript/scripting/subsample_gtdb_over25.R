#subsample gtdb randomly

library(tidyverse)
library(stringr)

str_reverse <- function(x) {
  sapply(lapply(strsplit(x, NULL), rev), paste, collapse = "")
}

d = read.delim('bac120_metadata_r214.tsv')
d2 = d %>% filter(gtdb_representative == 'f') %>% mutate(accession = gsub('RS_','',accession)) %>% mutate(accession = gsub('GB_','',accession)) 

genbank = read.delim('assembly_summary_genbank.txt',sep='\t') %>% dplyr::rename(accession = X.assembly_accession)
refseq = read.delim('assembly_summary_refseq.txt',sep='\t') %>% dplyr::rename(accession = X.assembly_accession)

ncbi = bind_rows(genbank %>% select(accession,ftp_path),refseq %>% select(accession,ftp_path)) %>% distinct

d3 = inner_join(d2,ncbi)

d4 = d3 %>% select(gtdb_taxonomy) %>% table %>% data.frame %>% filter(Freq>25) %>% select(gtdb_taxonomy)
d5 = inner_join(d4,d3)
d5 = d5 %>% group_by(gtdb_taxonomy) %>% sample_n(25) %>% separate(gtdb_taxonomy, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';', fill = 'right')

out = d5 %>% mutate(name = sapply(strsplit(ftp_path, "/"), function(x) str_reverse(str_reverse(x[length(x)])))) %>% mutate(out = paste(ftp_path, '/', name, '_genomic.fna.gz', sep = '')) %>% select(out)

write.table(out,'validation_to_download.csv',quote=F ,row.names=FALSE,col.names=FALSE)







