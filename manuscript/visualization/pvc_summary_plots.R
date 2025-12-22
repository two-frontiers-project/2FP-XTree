library(tidyverse)
library(cowplot)
library(RColorBrewer)

### TARGET OUTPUT PLOTS -- 
# counts of viruses by quality
# counts of phyla colored by quality
# counts of size of gene catalog barplot
# kegg overlap annotation

setwd('/Users/bradentemp/Library/CloudStorage/GoogleDrive-braden@twofrontiers.org/Shared drives/The Two Frontiers Project/STUDIES/PUBLICATION_2024-11-12_xtree_pvc')

# load pvc taxonomic and quality data
qual = read.delim('old/data_packet/taxonomy/pvc_checkv.tsv') %>% dplyr::rename(seq_name = contig_id)
tax = read.delim('old/data_packet/taxonomy/genomad_taxonomy_20231120_complete_high_medium_low_quality_taxonomy.tsv') 
mapping = read.delim('./old/data_packet/taxonomy/unique_terms_with_ranks_and_genome_type.csv',sep=',')

mapping_list <- split(mapping$Term, mapping$Rank)
mapping_list2 <- split(mapping$Term, mapping$Genome_Type)

# Function to extract the term for a given rank from a lineage string
extract_term <- function(lineage, rank) {
  terms <- strsplit(lineage, ";")[[1]]
  matched_terms <- terms[terms %in% mapping_list[[rank]]]
  if (length(matched_terms) > 0) return(matched_terms[1])
  return(NA)
}

extract_genome_composition <- function(lineage, mapping) {
  terms <- strsplit(lineage, ";")[[1]]
  genome_types <- mapping$Genome_Type[match(terms, mapping$Term)]
  genome_types <- na.omit(genome_types)  # Remove NA values
  
  if (any(genome_types == "RNA")) {
    return("RNA")
  } else if (any(genome_types == "DNA")) {
    return("DNA")
  } else if (any(genome_types == "Both")) {
    return("Both")
  } else {
    return("Unknown")
  }
}

# Create a dataframe mapping lineage to each of the ranks and genome composition
unique_lineages <- unique(tax$lineage)
lineage_rank_df <- data.frame(lineage = unique_lineages)

ranks <- unique(mapping$Rank)
for (rank in ranks) {
  lineage_rank_df[[rank]] <- sapply(unique_lineages, function(l) extract_term(l, rank))
}

# Add genome composition column
lineage_rank_df$Genome_Composition <- sapply(unique_lineages, function(l) extract_genome_composition(l, mapping))

tax = left_join(tax,lineage_rank_df)

qualtax = left_join(qual,tax)
qualtax$checkv_quality = factor(qualtax$checkv_quality,levels = c('Complete','High-quality','Medium-quality','Low-quality','Not-determined'))

#p# viruses by quality
ggplot(data = qualtax,aes(x = checkv_quality,fill=checkv_quality)) + geom_bar(stat = 'count') + theme_cowplot() + ylab('# of viruses') + xlab('') + scale_fill_brewer(palette='Set1') + theme(axis.text.x = element_text(angle= 45,hjust = 1))+ ggtitle('Viral counts by quality') + scale_y_log10() + theme(legend.position = 'none')
ggsave('plots/pvc_quality.pdf',width=3,height=4)

#p# taxonomy by quality
qualtax <- qualtax %>%group_by(Class) %>%mutate(Count = n()) %>%ungroup() %>%mutate(Class = reorder(Class, desc(Count)))

# Create a color mapping based on genome composition
color_mapping <- setNames(c("red", "blue", "green"), c("RNA", "DNA", "Both"))

# Create the plot
qualtaxclass = qualtax %>% select(Class,checkv_quality) %>% group_by(Class,checkv_quality) %>% count() %>% ungroup%>%mutate(Class = reorder(Class, desc(n)))
 
ggplot(data = qualtaxclass, aes(x = Class,y = log10(n), fill = checkv_quality)) +theme(legend.position = 'none')+geom_bar(stat = 'identity') + theme_cowplot() + ylab('# of viruses') + xlab('') + scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle('PVC taxonomies by quality') 
ggsave('plots/pvc_taxonomy_qual_genomad.pdf',width=6,height=3)

qualtaxcomp = qualtax %>% select(Genome_Composition,checkv_quality) %>% group_by(Genome_Composition,checkv_quality) %>% count() %>% ungroup%>%mutate(Genome_Composition = reorder(Genome_Composition, desc(n)))
ggplot(data = qualtaxcomp %>% filter(Genome_Composition == 'DNA' | Genome_Composition == 'RNA'), aes(x = Genome_Composition,y = log10(n), fill = checkv_quality)) +geom_bar(stat = 'identity') + theme_cowplot() + ylab('# of viruses') + xlab('') + scale_fill_brewer(palette = 'Set1') + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +ggtitle('PVC genome composition by quality') + theme(legend.position = 'none')
ggsave('plots/pvc_taxonomy_genome_comp.pdf',width=2,height=4)

### will need to repeat for ncbi and consensus taxonomy
  
# load functional data
library(tidyverse)

pgap = read.delim('old/data_packet/annotation/pvc_pgap_annotations_long.tsv')
pgapgenes = pgap %>% select(gene) %>% mutate(pgap = 1) %>% distinct

kegg = read.delim('old/data_packet/annotation/pvc_kegg_annotations_long.tsv')
kegggenes = kegg %>% select(gene)%>% mutate(kegg = 1)%>% distinct

pk = full_join(pgapgenes,kegggenes)

genecat = read.delim('old/data_packet/genecat/gene_cluster_analysis_results.csv',sep=',') 
genecat= genecat %>% melt
ggplot(data = genecat, aes(x = fct_reorder(PID, value), y = value, fill = factor(variable, levels = c('ClustersWithPgap', 'ClustersWithKegg', 'UniqueGenes')))) + geom_bar(stat = 'identity', position = 'dodge') + geom_text(aes(label = scales::comma(value)), position = position_dodge(width = 0.9), vjust = .1, hjust = 1,size = 3, check_overlap = TRUE) + scale_fill_brewer(palette = 'Set2') + scale_y_continuous(labels = scales::number_format()) + scale_y_log10() + theme_cowplot() + labs(y = "Value") + theme(legend.position = 'bottom',legend.title = element_blank()) + xlab('Clustering Identity') + ylab('') + coord_flip() + ggtitle('PVC gene catalog size and annotation data')
ggsave('plots/gene_catalog_stats.pdf',width=5,height=4)






