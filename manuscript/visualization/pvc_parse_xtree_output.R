library(tidyverse)

### ARGUMENTS
# input directory
# output directory
# database (bacterial/viral)
# thres
# uthres

args = commandArgs(trailingOnly=TRUE)

inputdir = "./metagenomics/prof/"#args[[1]]
outputdir = "xtree_parsed/"#args[[2]]
database = "pvc-viral"#args[[3]]
thres = .1#args[[4]]
uthres = 0.05#args[[5]]
namemap = "pvc_taxonomic_data_20240123.csv"#args[[6]]
dtype = "metagenomics"#args[[7]]

covs = data.frame()
covsU = data.frame()
otus = data.frame()
# covsR = data.frame()
files = dir(inputdir,pattern = "*.cov")
for (file in files) {
  t = read.delim(paste(inputdir,file,sep='/'),T,row=1,as.is=T)
  sname = gsub(".cov.*","",file)
  covs[rownames(t),sname] = t[,"Proportion_covered"]
  covsU[rownames(t),sname] = t[,"Unique_proportion_covered"]
  otus[rownames(t),sname] = pmin(t[,"Coverage_est"],t[,'XCov_adamantium'])
  #covsR[rownames(t),sname] = t[,"Reads_covered"]
}
covs[is.na(covs)]=0
covsU[is.na(covsU)]=0
otus[is.na(otus)]=0

# covsR[is.na(covsR)] = 0

#otus = data.frame()
#files = dir(inputdir,pattern = "*.ref")
#for (file in files) {
#  tryCatch({
#  t = read.delim(paste(inputdir,file,sep='/'),F,row=1,as.is=T)
#  },
#  error=function(cond){
#  print('No lines in file')
#  })
#  otus[rownames(t),gsub(".ref.*","",file)] = t
#}
#otus[is.na(otus)]=0

# Prepare the coverage + ref cross-maps with taxonomy
if(namemap != 'None'){
  tkey= read.delim(namemap,header=T,sep=',')
  tkey = tkey %>% select(query,species_rep_90perc,genomad_lineage,ncbi_organism_name,ncbi_taxid) %>% distinct %>% column_to_rownames('query')
  covs = covs[rownames(covs) %in% rownames(tkey), ,drop=F ]
  covsU = covsU[rownames(covsU) %in% rownames(tkey),,drop=F ]
  taxaconv = tkey[rownames(covs),1] 
  otus = otus[rownames(otus) %in% rownames(tkey),,drop=F ]
  taxaconvR = tkey[rownames(otus),1]
}
if(namemap == 'None'){
  taxaconv = rownames(covs)
  taxaconvR = rownames(otus)
}

orig.abund = colSums(otus)

mask = (covsU <= as.numeric(uthres) | covs <= as.numeric(thres)) #& covs <= hthres
unmasked = rowsum(0+!mask,taxaconv)

tax.tm = rowsum(otus,taxaconvR)
if(ncol(tax.tm) == 1){
  tax.tm[!unmasked[rownames(tax.tm),],]=0
}
if(ncol(tax.tm)>1){
  tax.tm[!unmasked[rownames(tax.tm),]]=0
}
tax.tm = tax.tm[rowSums(tax.tm) > 0,,drop=F]
tax.tm["Unknown",]=orig.abund - colSums(tax.tm)
tax.cm = rowsum(otus,taxaconvR)[rownames(tax.tm)[-nrow(tax.tm)],,drop=F]
tax.tmr=sweep(tax.tm,2,colSums(tax.tm),'/')

tax.tm = tax.tm %>% rownames_to_column('query')
tax.cm = tax.cm %>% rownames_to_column('query')
tax.tmr = tax.tmr %>% rownames_to_column('query')

tax.tm = inner_join(tax.tm,tkey %>% rownames_to_column('query'))
tax.cm = inner_join(tax.cm,tkey %>% rownames_to_column('query'))
tax.tmr = inner_join(tax.tmr,tkey %>% rownames_to_column('query'))

write.table(tax.tm,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',dtype,'_','species_counts.tsv',sep=''),F,F,'\t')
write.table(tax.tmr,paste(outputdir,'/',database,'_',as.character(thres),'_',as.character(uthres),'_',dtype,'_','species_ra.tsv',sep=''),F,F,'\t')
write.table(tax.cm,paste(outputdir,'/',database,'_',dtype,'_','species_rawcounts.tsv',sep=''),F,F,'\t')





