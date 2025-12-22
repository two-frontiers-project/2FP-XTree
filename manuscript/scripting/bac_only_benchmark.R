library(tidyverse)
library(ggplot2)
library(reshape2)
library(stringi)
### post process viral only benchmarking

# functions

loaddata <- function(inputdir){
	merged= list()
	files = dir(inputdir,pattern = "*.cov")
	for (file in files) {
	  t = read.delim(paste(inputdir,file,sep='/'),T,row=1,as.is=T)
	  t = t %>% rownames_to_column('genome_id') %>% select(genome_id,Proportion_covered,Unique_proportion_covered,Adamantium_covered,Adamantium_prop,Coverage_est,XCov_adamantium) %>% mutate(observed_ra = pmin(Coverage_est,XCov_adamantium,na.rm=T))%>% mutate(observed_ra = observed_ra/sum(observed_ra,na.rm=T))
	  sname = gsub(".cov.*","",file)
	  t = t %>% mutate(sample = sname) %>% mutate(abset = inputdir%>% gsub('-V2','',.)%>% gsub('-v2','',.)  %>% stri_reverse() %>% strsplit('_') %>% map_chr(2) %>% stri_reverse()) %>% mutate(genomeset = inputdir %>% gsub('-V2','',.)%>% gsub('-v2','',.) %>% stri_reverse() %>% strsplit('_') %>% map_chr(1) %>% stri_reverse()) 
	  merged[[sname]] = t 
	}
	merged = bind_rows(merged)
	otus = list()
	files = dir(inputdir,pattern = "*.ref")
	for (file in files) {
	  t = read.delim(paste(inputdir,file,sep='/'),T,row=1,as.is=T)
	  colnames(t)[1] = 'otus'
	  t = t %>% rownames_to_column('genome_id') 
	  sname = gsub(".ref.*","",file)
	  t = t %>% mutate(sample = sname) 
	  otus[[sname]] = t
	}
	otus = bind_rows(otus)	
	merged = inner_join(merged,otus)
	return(merged)
}

# directories with files

dirs=c("xtree_output_bac_21_aligned_to_29_highab_NONREP","xtree_output_bac_21_aligned_to_29_highab_REP","xtree_output_bac_21_aligned_to_29_lowab_NONREP","xtree_output_bac_21_aligned_to_29_lowab_REP","xtree_output_bac_24_aligned_to_29_highab_NONREP","xtree_output_bac_24_aligned_to_29_highab_REP","xtree_output_bac_24_aligned_to_29_lowab_NONREP","xtree_output_bac_24_aligned_to_29_lowab_REP","xtree_output_bac_29_highab_NONREP","xtree_output_bac_29_highab_REP","xtree_output_bac_29_lowab_NONREP","xtree_output_bac_29_lowab_REP")

# load .cov and .ref files, parse into long form, compute relative abundance
parseddat = list()
for(d in dirs){
	loaded=loaddata(d) %>% mutate(kmer = d %>% strsplit('_') %>% map_chr(4)) 
	parseddat[[d]] = loaded
}

parsed = parseddat %>% bind_rows %>% mutate(mock_commmunity_id = strsplit(sample,'genomeset_','') %>% map_chr(2) %>% gsub('_COV_0','',.))
parsed[is.na(parsed)] = 0
parsed$sample = gsub('1_genomeset_','',parsed$sample)
parsed$sample = gsub('bacterialonly_gtdb214-NONREP','',parsed$sample)
parsed$sample = gsub('bacterialonly_gtdb214-REP','',parsed$sample)
parsed$sample = gsub('bacterialonly_gtdb214','',parsed$sample)
parsed$sample = gsub('_Background','',parsed$sample)
parsed$sample = gsub('_COV_0','',parsed$sample)
parsed$observed_ra[is.na(parsed$observed_ra)] = 0

truedata = list()
# load true coverage data of genomes used in community
#dirs2=c("xtree_output_bac_29_highab_NONREP","xtree_output_bac_29_highab_REP","xtree_output_bac_29_lowab_NONREP","xtree_output_bac_29_lowab_REP")
for(d in dirs){
	print(d)
	d2=gsub('_21_aligned_to','',d)
	d2=gsub('_24_aligned_to','',d2)
	d2=gsub('_17_aligned_to','',d2)
	files = dir(d2,pattern = "Background\\.genomeset.*log")
	kval=strsplit(d,'_') %>% map_chr(4)
	for(f in files){
		samplesetname = gsub(d,'',f)
		samplesetname = gsub('/','',f)
		samplesetname = gsub('Background\\.','',f)
		samplesetname = gsub('\\.log','',samplesetname)
		dat = read.delim(paste0(d2,'/',f),sep='\t') %>% mutate(sample = samplesetname)
		fasta_files <- list.files(paste0('.','/',samplesetname), pattern = "\\.fasta$", full.names = TRUE)
		identifiers <- c()
		file_names <- c()
		for(file in fasta_files) {
		    first_line <- readLines(file, n = 1)
		    identifier <- substring(first_line, 2)
		    identifiers <- c(identifiers, identifier)
		    file_names <- c(file_names, tools::file_path_sans_ext(basename(file)))
		}
		result_df <- data.frame(genome_id = identifiers, mock_commmunity_id = file_names)
		dat = inner_join(dat %>% select(-num_background) %>% rename(mock_commmunity_id = Genome),result_df) %>% mutate(kmer = kval) %>% mutate(abset = d %>% gsub('-V2','',.)%>% gsub('-v2','',.) %>% stri_reverse() %>% strsplit('_') %>% map_chr(2) %>% stri_reverse()) %>% mutate(genomeset = d %>% gsub('-V2','',.)%>% gsub('-v2','',.)%>% stri_reverse() %>% strsplit('_') %>% map_chr(1) %>% stri_reverse()) 
		truedata[[paste(f,d,kval)]] = dat 
	}
}

truedata = bind_rows(truedata)
#truedata2 = truedata %>% mutate(kmer = "24")
#truedata3 = truedata %>% mutate(kmer = "21")

#truedata = bind_rows(truedata,truedata2,truedata3)

truedata$sample = gsub('bacterialonly_gtdb214-NONREP','',truedata$sample)
truedata$sample = gsub('bacterialonly_gtdb214-REP','',truedata$sample)
truedata$sample = gsub('bacterialonly_gtdb214','',truedata$sample)
truedata$sample = gsub('genomeset_','',truedata$sample)

# load in full database info
gtdb = read.delim('bac120_metadata_r214.tsv',header=T) %>% dplyr::rename(genome_id = accession)
gtdb$genome_id = gsub('RS_','',gtdb$genome_id)
gtdb$genome_id = gsub('GB_','',gtdb$genome_id)
gtdb = gtdb %>% select(genome_id,gtdb_taxonomy)
gtdb = gtdb %>% mutate(genome_id = genome_id %>% strsplit('\\.') %>% map_chr(1))
colnames(gtdb) = c('genome_id','lineage')

tndat = dim(gtdb)[[1]]

# merge true positive data
parsed = parsed %>% mutate(genome_id = genome_id %>% strsplit('\\.') %>% map_chr(1)) 
parsed = left_join(parsed,gtdb) %>% separate(lineage, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';', fill = 'right') %>% select(-Kingdom,-Phylum,-Class,-Order,-Family)
parsed = parsed %>% filter(!is.na(Species),!is.na((Genus)))

parsed$Species = as.factor(parsed$Species)
parsed$Genus = as.factor(parsed$Genus)
parsed$sample = as.factor(parsed$sample)
parsed$genomeset = as.factor(parsed$genomeset)
parsed$abset = as.factor(parsed$abset)

library(data.table)
parsed_dt = as.data.table(parsed)
result = parsed_dt[, lapply(.SD, mean, na.rm = TRUE), by = .(Species, sample, abset, genomeset, kmer), .SDcols = sapply(parsed_dt, is.numeric)]

parsed_s = as.data.frame(result) %>% mutate(rank = 'Species')
colnames(parsed_s)[1] = 'taxonomy'

#result = parsed_dt[, lapply(.SD, mean, na.rm = TRUE), by = .(Genus, sample, abset, genomeset, kmer), .SDcols = sapply(parsed_dt, is.numeric)]

#parsed_g = as.data.frame(result) %>% mutate(rank = 'Genus')
#colnames(parsed_g)[1] = 'taxonomy'

parsed_rank = bind_rows(parsed_s)

truedata = truedata %>% mutate(genome_id = genome_id %>% strsplit('\\.') %>% map_chr(1))%>% select(-mock_commmunity_id)
truedata = left_join(truedata,gtdb) %>% separate(lineage, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';', fill = 'right') %>% select(-Kingdom,-Phylum,-Class,-Order,-Family)
###
truedata$Species = as.factor(truedata$Species)
truedata$Genus = as.factor(truedata$Genus)
truedata$sample = as.factor(truedata$sample)
truedata$genomeset = as.factor(truedata$genomeset)
truedata$abset = as.factor(truedata$abset)

truedata_dt = as.data.table(truedata)
result = truedata_dt[, lapply(.SD, mean, na.rm = TRUE), by = .(Species, sample, abset, genomeset, kmer), .SDcols = sapply(truedata_dt, is.numeric)]

truedata_s = as.data.frame(result) %>% mutate(rank = 'Species')
colnames(truedata_s)[1] = 'taxonomy'

#result = truedata_dt[, lapply(.SD, mean, na.rm = TRUE), by = .(Genus, sample, abset, genomeset, kmer), .SDcols = sapply(truedata_dt, is.numeric)]

#truedata_g = as.data.frame(result) %>% mutate(rank = 'Genus')
#colnames(truedata_g)[1] = 'taxonomy'

truedata_rank = bind_rows(truedata_s)

###

parsed_tp = full_join(parsed_rank,truedata_rank) %>% dplyr::rename(expected_ra = Coverage)

parsed_tp = parsed_tp %>% mutate(expected_ra = if_else(is.na(expected_ra), 0,expected_ra)) %>% mutate(observed_ra = if_else(is.na(observed_ra), 0,observed_ra)) %>% filter(!(expected_ra == 0 & observed_ra == 0))

parsed_tp = parsed_tp %>% mutate(tp_fp = 'FALSE POSITIVE',tp_fp = if_else(expected_ra>0 & observed_ra>0,'TRUE POSITIVE',if_else(expected_ra==0 & observed_ra>0,'FALSE POSITIVE',if_else(expected_ra>0 &observed_ra == 0,'FALSE NEGATIVE',tp_fp))))


parsed_tp = parsed_tp %>% mutate(dbsize = tndat)
#parsed_tp = left_join(parsed_tp,gtdb)

### CHECK THAT YOU HAVE 200 UNIQUE SAMPLES ACROSS SAMPLEID, ABSET,GENOMESET 
# parsed %>% select(sample,abset,genomeset) %>% distinct
# parsed_tp %>% select(sample,abset,genomeset) %>% distinct
# truedata %>% select(sample,abset,genomeset) %>% distinct

saveRDS(parsed_tp,'data_packet/comparison_data_processed/gtdbr214_benchmarking_xtree_full_dataset_auc.rds')
























