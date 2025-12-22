# standardized plots for true false positive rate PVC-X

library(tidyverse)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(PRROC)
library(fuzzyjoin)

compute_pr_auc <- function(data) {
  # Initialize dataframes to store results for total and unique coverage
  total_results <- tibble(cutoff = double(), Precision = double(), Recall = double(), AUC_PR = double())
  unique_results <- tibble(cutoff = double(), Precision = double(), Recall = double(), AUC_PR = double())
  
  # Iterate over total coverage cutoffs in increments of 0.01
  for (cutoff in seq(0, 1, by = 0.01)) {
    
    # Filter and calculate for total coverage
    filtered_data_total <- data |>
      filter(!is.na(sample)) |>
      mutate(tp_fp = if_else(Proportion_covered < cutoff & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) |>
      mutate(tp_fp = if_else(Proportion_covered < cutoff & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) |> 
      filter(tp_fp != 'TRUE NEGATIVE') |>
      filter(Proportion_covered >= cutoff | tp_fp %in% c('FALSE NEGATIVE')) |>
      group_by(sample) |>
      summarise(TP = sum(tp_fp == "TRUE POSITIVE"), 
                FP = sum(tp_fp == "FALSE POSITIVE"),
                FN = sum(tp_fp == "FALSE NEGATIVE")) |>
      mutate(Precision = TP / (TP + FP),
             Recall = TP / (TP + FN))
    
    # Aggregate and store total coverage results
    aggregate_data_total <- filtered_data_total |>
      summarise(Avg_Precision = mean(Precision, na.rm = TRUE), 
                Avg_Recall = mean(Recall, na.rm = TRUE))
    
    pr_curve_total <- pr.curve(scores.class0 = filtered_data_total$FP, scores.class1 = filtered_data_total$TP)
    auc_pr_total <- pr_curve_total$auc.integral
    
    total_results <- total_results |>
      add_row(cutoff = cutoff, Precision = aggregate_data_total$Avg_Precision, Recall = aggregate_data_total$Avg_Recall, AUC_PR = auc_pr_total)
  }
  
  # Iterate over unique coverage cutoffs in increments of 0.01
  for (cutoff in seq(0, 1, by = 0.01)) {
    
    # Filter and calculate for unique coverage
    filtered_data_unique <- data |>
      filter(!is.na(sample)) |>
      mutate(tp_fp = if_else(unique_coverage_union < cutoff & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) |>
      mutate(tp_fp = if_else(unique_coverage_union < cutoff & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) |> 
      filter(tp_fp != 'TRUE NEGATIVE') |>
      filter(unique_coverage_union >= cutoff | tp_fp %in% c('FALSE NEGATIVE')) |>
      group_by(sample) |>
      summarise(TP = sum(tp_fp == "TRUE POSITIVE"), 
                FP = sum(tp_fp == "FALSE POSITIVE"),
                FN = sum(tp_fp == "FALSE NEGATIVE")) |>
      mutate(Precision = TP / (TP + FP),
             Recall = TP / (TP + FN))
    
    # Aggregate and store unique coverage results
    aggregate_data_unique <- filtered_data_unique |>
      summarise(Avg_Precision = mean(Precision, na.rm = TRUE), 
                Avg_Recall = mean(Recall, na.rm = TRUE))
    
    pr_curve_unique <- pr.curve(scores.class0 = filtered_data_unique$FP, scores.class1 = filtered_data_unique$TP)
    auc_pr_unique <- pr_curve_unique$auc.integral
    
    unique_results <- unique_results |>
      add_row(cutoff = cutoff, Precision = aggregate_data_unique$Avg_Precision, Recall = aggregate_data_unique$Avg_Recall, AUC_PR = auc_pr_unique)
  }
  
  return(list(total_results = total_results, unique_results = unique_results))
}
setwd('/Users/bradentemp/Library/CloudStorage/GoogleDrive-braden@twofrontiers.org/Shared drives/The Two Frontiers Project/STUDIES/PUBLICATION_2024-11-12_xtree_pvc')

tax = read.delim('data_packet/pvc_tax/pvc_taxonomic_data_20240123.csv',sep=',') 
tax = tax %>% select(-X) %>% dplyr::rename(genome_id = query)
tax2 = read.delim('data_packet/pvc_tax/unique_terms_with_ranks_and_genome_type.csv',sep=',')
colnames(tax2) = c('Term','genomad Rank','genomad Genome_Type')

viral = readRDS('data_packet/comparison_data_processed/viral_benchmarking_xtree_full_dataset.rds') %>% mutate(type = 'All Genomes')
viral = inner_join(viral,tax)
parsed_tp = viral%>% mutate(type = 'All Genomes') %>% mutate(unique_coverage_union = pmin(Adamantium_prop,Unique_proportion_covered,na.rm=T))
parsed_tp$genomeset = parsed_tp$dbqual
parsed_tp = parsed_tp %>% mutate(genomeset = if_else(grepl('nonrep',sample),'PVC Non-Representative Genomes',genomeset))
parsed_tp$genomeset = gsub('complete_high_medium_low','CHML',parsed_tp$genomeset)
parsed_tp$genomeset = gsub('complete_high_medium','CHM',parsed_tp$genomeset)
parsed_tp$genomeset = gsub('complete_high','CH',parsed_tp$genomeset)
parsed_tp$dbqual[parsed_tp$dbqual == 'nonrep_V2'] = 'complete_high_medium_low'

viral_phy <- regex_inner_join(parsed_tp,tax2 %>% filter(`genomad Rank` == "Phylum"), by = c("lineage" = "Term")) %>% mutate(type = 'Phylum') %>% select(-lineage) %>% dplyr::rename(lineage = Term) %>% group_by(lineage,genomeset,sample,dbqual,kmer,type) %>% summarise(Proportion_covered = mean(Proportion_covered,na.rm=T),unique_coverage_union = mean(unique_coverage_union,na.rm=T),observed_ra = sum(observed_ra,na.rm=T),expected_ra = sum(expected_ra,na.rm=T)) %>% mutate(tp_fp = 'FALSE POSITIVE',tp_fp = if_else(expected_ra>0 & observed_ra>0,'TRUE POSITIVE',if_else(expected_ra==0 & observed_ra>0,'FALSE POSITIVE',if_else(expected_ra>0 &observed_ra == 0,'FALSE NEGATIVE',tp_fp))))
viral_cla <- regex_inner_join(parsed_tp,tax2 %>% filter(`genomad Rank` == "Class"), by = c("lineage" = "Term")) %>% mutate(type = 'Class')  %>% select(-lineage) %>% dplyr::rename(lineage = Term) %>% group_by(lineage,genomeset,sample,dbqual,kmer,type) %>% summarise(Proportion_covered = mean(Proportion_covered,na.rm=T),unique_coverage_union = mean(unique_coverage_union,na.rm=T),observed_ra = sum(observed_ra,na.rm=T),expected_ra = sum(expected_ra,na.rm=T)) %>% mutate(tp_fp = 'FALSE POSITIVE',tp_fp = if_else(expected_ra>0 & observed_ra>0,'TRUE POSITIVE',if_else(expected_ra==0 & observed_ra>0,'FALSE POSITIVE',if_else(expected_ra>0 &observed_ra == 0,'FALSE NEGATIVE',tp_fp))))
viral_ord <- regex_inner_join(parsed_tp,tax2 %>% filter(`genomad Rank` == "Order"), by = c("lineage" = "Term")) %>% mutate(type = 'Order')  %>% select(-lineage) %>% dplyr::rename(lineage = Term) %>% group_by(lineage,genomeset,sample,dbqual,kmer,type) %>% summarise(Proportion_covered = mean(Proportion_covered,na.rm=T),unique_coverage_union = mean(unique_coverage_union,na.rm=T),observed_ra = sum(observed_ra,na.rm=T),expected_ra = sum(expected_ra,na.rm=T)) %>% mutate(tp_fp = 'FALSE POSITIVE',tp_fp = if_else(expected_ra>0 & observed_ra>0,'TRUE POSITIVE',if_else(expected_ra==0 & observed_ra>0,'FALSE POSITIVE',if_else(expected_ra>0 &observed_ra == 0,'FALSE NEGATIVE',tp_fp))))
viral_fam <- regex_inner_join(parsed_tp,tax2 %>% filter(`genomad Rank` == "Family"), by = c("lineage" = "Term")) %>% mutate(type = 'Family')  %>% select(-lineage) %>% dplyr::rename(lineage = Term) %>% group_by(lineage,genomeset,sample,dbqual,kmer,type) %>% summarise(Proportion_covered = mean(Proportion_covered,na.rm=T),unique_coverage_union = mean(unique_coverage_union,na.rm=T),observed_ra = sum(observed_ra,na.rm=T),expected_ra = sum(expected_ra,na.rm=T)) %>% mutate(tp_fp = 'FALSE POSITIVE',tp_fp = if_else(expected_ra>0 & observed_ra>0,'TRUE POSITIVE',if_else(expected_ra==0 & observed_ra>0,'FALSE POSITIVE',if_else(expected_ra>0 &observed_ra == 0,'FALSE NEGATIVE',tp_fp))))
parsed_tp = bind_rows(parsed_tp,viral_phy,viral_cla,viral_ord,viral_fam)

parsed_tp$tp_fp <- factor(parsed_tp$tp_fp, levels = c("FALSE POSITIVE", "TRUE POSITIVE","FALSE NEGATIVE","TRUE NEGATIVE"))

parsed_tp$type = factor(parsed_tp$type,levels = c('All Genomes','Phylum','Class','Order','Family','NCBI Species'))

kraken2 = readRDS('./data_packet/comparison_data_processed/kraken2_benchmarking_vir_full_dataset_auc_withkmers.rds')
kraken2$genomeset[kraken2$genomeset == 'nonrep'] = 'PVC Non-Representative Genomes'

viral_f1_early <- parsed_tp %>% filter(kmer == 21) %>% mutate(tp_fp = if_else(unique_coverage_union<0.01 & tp_fp == 'TRUE POSITIVE','FALSE NEGATIVE',tp_fp)) %>% mutate(tp_fp = if_else(unique_coverage_union<0.01 & tp_fp == 'FALSE POSITIVE','TRUE NEGATIVE',tp_fp)) %>% mutate(plot_type = case_when(type == 'All Genomes' & genomeset %in% c('CH','CHM','CHML') ~ 'Rep', type == 'All Genomes' & genomeset == 'PVC Non-Representative Genomes' ~ 'Non-Rep', type %in% c('Phylum','Class','Order','Family') & genomeset == 'PVC Non-Representative Genomes' ~ 'Taxonomy', TRUE ~ NA_character_)) %>% filter(!is.na(plot_type)) %>% group_by(plot_type, sample) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE"), .groups = 'drop') %>% mutate(Precision = TP / (TP + FP), Recall = TP / (TP + FN), F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% filter(!is.na(TP), !is.na(Precision), !is.na(Recall)) %>% group_by(plot_type) %>% summarise(mean_precision = mean(Precision, na.rm=TRUE), sd_precision = sd(Precision, na.rm=TRUE), mean_recall = mean(Recall, na.rm=TRUE), sd_recall = sd(Recall, na.rm=TRUE), mean_f1 = mean(F1, na.rm=TRUE), sd_f1 = sd(F1, na.rm=TRUE), .groups = 'drop')
print("Viral Xtree F1 scores:")
print(viral_f1_early)

kraken_f1_early <- kraken2 %>% filter(!is.na(rank), genomeset == 'PVC Non-Representative Genomes') %>% group_by(sample,rank) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE"), .groups = 'drop') %>% mutate(Precision = TP / (TP + FP), Recall = TP / (TP + FN), F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% dplyr::rename(type = rank) %>% mutate(plot_type = case_when(type == 'All Genomes' ~ 'Non-Rep', type %in% c('Phylum','Class','Order','Family') ~ 'Taxonomy', TRUE ~ NA_character_)) %>% filter(!is.na(plot_type), plot_type %in% c('Non-Rep','Taxonomy'), !is.na(TP), !is.na(Precision), !is.na(Recall)) %>% group_by(plot_type) %>% summarise(mean_precision = mean(Precision, na.rm=TRUE), sd_precision = sd(Precision, na.rm=TRUE), mean_recall = mean(Recall, na.rm=TRUE), sd_recall = sd(Recall, na.rm=TRUE), mean_f1 = mean(F1, na.rm=TRUE), sd_f1 = sd(F1, na.rm=TRUE), .groups = 'drop')
print("Viral Kraken2 F1 scores:")
print(kraken_f1_early)

# get precision recall curves
curvedatalist1 = list()
curvedatalist2 = list()
for(t in 'All Genomes'){
for(num in c(17,21,24,29)){
    for(gset in unique(parsed_tp$genomeset)){
        a = compute_pr_auc(parsed_tp %>% filter(kmer == num,genomeset == gset,type ==t)) 
        a1 = a[[1]] %>% mutate(curvetype = 'P-R',kmer = num,covtype = 'Total Coverage',genomeset = gset,type = t)
        a2 = a[[2]] %>% mutate(curvetype = 'P-R',kmer = num,covtype = 'Unique Coverage',genomeset = gset,type = t)
        # b = compute_tpr_fpr(parsed_tp %>% filter(kmer == num,rank == r,abset == ab,genomeset == gset))
        #  b1 = b[[1]] %>% mutate(curvetype = 'TPR-FPR',kmer = num,rank = r,covtype = 'Total Coverage',abset = ab,genomeset = gset)
        # b2 = b[[2]] %>% mutate(curvetype = 'TPR-FPR',kmer = num,rank = r,covtype = 'Unique Coverage',abset = ab,genomeset = gset)
        curvedatalist1[[paste(num,gset,t)]] = bind_rows(a1,a2)
        #  curvedatalist2[[paste(num,'b',gset,r,ab)]] = bind_rows(b1,b2)     
    }
  }
}

curvedata = bind_rows(curvedatalist1)
#curvedata2 = bind_rows(curvedatalist2) %>% filter(rank=='Species')

# plot out the curve data and pick your cutoffs

# pick cutoffs

optimize_cutoff <- function(data) {
  data %>%
    mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>%
    arrange(desc(F1)) %>%
    slice(1)
}

optimized_df <- curvedata %>% group_by(curvetype, kmer, covtype, genomeset, type) %>% do(optimize_cutoff(.))%>% select(curvetype,kmer,covtype,cutoff, type,genomeset) %>% mutate(cutoffselect=TRUE)

curvedata = left_join(curvedata,optimized_df) %>% mutate(cutoffselect = if_else(cutoffselect == TRUE,cutoffselect,FALSE))

#### overall f1

curvdatatoplot = curvedata %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall))
ggplot(curvdatatoplot %>% filter(type == 'All Genomes'),aes(x = cutoff,y = F1, shape = type,color = factor(kmer))) + geom_point(size=2.5,alpha=.5) + theme_cowplot() + facet_grid(covtype ~ genomeset) + scale_color_brewer(palette = 'Set3') + theme(legend.position = 'bottom')
ggsave('plots/viral_f1_plot.pdf',width=8,height=4)

# full 
cd2 = curvedata %>% filter(type == 'All Genomes')
ggplot() + theme_cowplot() + geom_point(alpha=0.5,data = cd2 %>% ungroup %>% filter(genomeset == 'PVC Non-Representative Genomes',covtype == 'Total Coverage',kmer == 21), aes(x = Precision, y = Recall, color = cutoff)) + facet_grid(. ~ type)  + scale_color_viridis_c(option = 'E') + geom_text(data = filter(cd2 %>% filter(genomeset == 'PVC Non-Representative Genomes',covtype == 'Total Coverage',kmer == 21), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.5), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(cd2%>% filter(genomeset == 'PVC Non-Representative Genomes',covtype == 'Total Coverage',kmer == 21), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.5), colour = "grey50")+ ggnewscale::new_scale_color()+ geom_point(alpha=0.5,data = cd2 %>% ungroup %>% filter(genomeset == 'PVC Non-Representative Genomes',covtype == 'Unique Coverage',kmer == 21), aes(x = Precision, y = Recall, color = cutoff))   + scale_color_viridis_c(option = 'G') + geom_text(data = filter(cd2 %>% filter(genomeset == 'PVC Non-Representative Genomes',covtype == 'Unique Coverage',kmer == 21), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(cd2%>% filter(genomeset == 'PVC Non-Representative Genomes',covtype == 'Unique Coverage',kmer == 21), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50") + theme(legend.position = 'bottom')
ggsave('plots/viral_k21_both_coverages_nonrep.pdf',width=4,height=4)


ggplot() + theme_cowplot() + geom_point(alpha=0.5,data = curvedata %>% ungroup %>% filter(covtype == 'Total Coverage',kmer == 21), aes(x = Precision, y = Recall, color = cutoff)) + facet_grid(type ~ genomeset)  + scale_color_viridis_c(option = 'E') + geom_text(data = filter(curvedata %>% filter(covtype == 'Total Coverage',kmer == 21), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.5), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(covtype == 'Total Coverage',kmer == 21), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.5), colour = "grey50")+ ggnewscale::new_scale_color()+ geom_point(alpha=0.5,data = curvedata %>% ungroup %>% filter(covtype == 'Unique Coverage',kmer == 21), aes(x = Precision, y = Recall, color = cutoff))   + scale_color_viridis_c(option = 'G') + geom_text(data = filter(curvedata %>% filter(covtype == 'Unique Coverage',kmer == 21), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(covtype == 'Unique Coverage',kmer == 21), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50") + theme(legend.position = 'bottom')
ggsave('plots/SUPP_viral_both_coverages_nonrep.pdf',width=12,height=12)


# merge in optimized cutoff data 
parsed_tp_op = left_join(parsed_tp,optimized_df %>% reshape2::dcast(curvetype + type + kmer +genomeset  ~ covtype,value.var = 'cutoff' ) %>% filter(type == 'All Genomes') %>% select(-type) %>% filter(genomeset == 'PVC Non-Representative Genomes') %>% select(-genomeset) %>% mutate(kmer = as.factor(kmer)))

# plot raw data on coverage cutoffs
#ggplot(parsed_tp_op |> filter(tp_fp != 'FALSE NEGATIVE',kmer!=21)|> filter(!is.na(kmer)) |> group_by(kmer) ,aes(x = Proportion_covered, y = unique_coverage_union)) + geom_point(aes(color= tp_fp),alpha=.5,size=.6)  + theme_cowplot() + geom_hline(aes(yintercept = `Unique Coverage`),color = 'red',linetype = 'dashed')+geom_vline(aes(xintercept = `Total Coverage`),color = 'red',linetype = 'dashed') + facet_grid(kmer + type ~ genomeset,scales = 'free')  + scale_fill_brewer(palette = 'Dark2')+ scale_color_brewer(palette = 'Dark2') + ylab('Unique coverage') + xlab('Total coverage') + xlim(0,1) + ylim(0,1)
#ggsave('plots/SUPP_viral_fp_rate_by_coverage_cutoff_sp.pdf',width=12,height=9)

parsed_tp_op$genomeset[parsed_tp_op$genomeset == 'PVC Non-Representative Genomes'] = 'Non-Representative'  
ggplot(parsed_tp_op |> filter(tp_fp != 'FALSE NEGATIVE',kmer==21) |> filter(!is.na(kmer)) |> group_by(kmer) ,aes(x = Proportion_covered, y = unique_coverage_union)) + geom_point(aes(color= tp_fp),alpha=.9,size=.6)  + theme_cowplot() + geom_hline(aes(yintercept = `Unique Coverage`),color = 'red',linetype = 'dashed')+geom_vline(aes(xintercept = `Total Coverage`),color = 'red',linetype = 'dashed')  + scale_fill_manual(values = c('grey','black'))+ scale_color_manual(values = c('grey','black'))+ ylab('Unique coverage') + xlab('Total coverage') + theme(legend.position = 'bottom') + facet_grid(type ~ genomeset)+ xlim(0,1) + ylim(0,1)
ggsave('plots/viral_k21_fp_rate_by_coverage_cutoff_sp_nonrep.png',width=10,height=8)

# compute FP/TP rates etc by points - CREATE plot_type BEFORE grouping to match expected vs observed structure
precision_recall_viral_opt <-  parsed_tp_op %>% mutate(tp_fp = if_else(unique_coverage_union<0.01 & tp_fp == 'TRUE POSITIVE','FALSE NEGATIVE',tp_fp))%>% mutate(tp_fp = if_else(unique_coverage_union<0.01 & tp_fp == 'FALSE POSITIVE','TRUE NEGATIVE',tp_fp)) %>%
  mutate(plot_type = case_when(
    type == 'All Genomes' & genomeset %in% c('CH','CHM','CHML') ~ 'All Genomes',
    type == 'All Genomes' & genomeset == 'Non-Representative' ~ 'All Genomes',
    type %in% c('Phylum','Class','Order','Family') ~ type,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(plot_type)) %>%
  group_by(kmer,genomeset,plot_type, sample) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE"), .groups = 'drop') %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))

precision_recall_kraken <- kraken2 %>% filter(!is.na(rank))%>% group_by(sample,rank) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN)) %>% dplyr::rename(type = rank) %>%
  mutate(plot_type = case_when(
    type == 'All Genomes' ~ 'Non-Rep',
    type %in% c('Phylum','Class','Order','Family') ~ type,
    TRUE ~ NA_character_
  )) %>%
  filter(!is.na(plot_type)) %>%
  # Kraken2 only has non-representative data, so filter out any that might not match
  filter(plot_type %in% c('Non-Rep','Phylum','Class','Order','Family'))

longkrak <- tidyr::pivot_longer(precision_recall_kraken, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "Kraken2")
long_df_opt <- tidyr::pivot_longer(precision_recall_viral_opt, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "Xtree (1% coverage cutoff)")

combined_data = bind_rows(longkrak,long_df_opt)
combined_data$aligner = factor(combined_data$aligner,levels = c('Xtree (1% coverage cutoff)','Kraken2'))
# Factor levels - Representative only appears in Xtree data, not Kraken2
combined_data$plot_type = factor(combined_data$plot_type, levels = c('Rep','Non-Rep','Phylum','Class','Order','Family'))

ggplot(data =  combined_data, aes(x = metric, y = value, color = aligner)) + theme_cowplot() + geom_violin(alpha = 0.2, scale = "width", position = position_dodge(width = 0.8)) + ggbeeswarm::geom_quasirandom(width = 0.2, dodge.width = 0.8, alpha=0.5) + stat_summary(fun = median, geom = "crossbar", width = 1, position = position_dodge(width = 0.8),aes(color = aligner)) + facet_wrap(~ plot_type, nrow = 1) + ylim(0, 1)+ scale_fill_manual(values = c('#00008B','#8B0000'))  + scale_color_manual(values = c('#00008B','#8B0000')) + theme(legend.position = 'bottom') 
ggsave('../../plots/viral_precision_recall_by_sample.pdf',width=8,height=3.5)

### FOR PUBLICATION
aligner_summary_table = bind_rows(precision_recall_viral_opt %>% mutate(aligner = "xtree (1% coverage cutoff)") ,precision_recall_kraken %>% mutate(aligner = "kraken2", kmer = NA, genomeset = 'Unknown')) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% filter(!is.na(TP)) %>% group_by(kmer,genomeset,plot_type,aligner)%>% summarise(mean_precision = mean(Precision),sd_precision = sd(Precision),mean_recall = mean(Recall),sd_recall = sd(Recall),mean_f1 = mean(F1),sd_f1 = sd(F1), .groups = 'drop')
write.csv(aligner_summary_table,'tables/aligner_summary_table_vir_AVG.csv')

aligner_summary_table = bind_rows(precision_recall_viral_opt %>% mutate(aligner = "xtree (1% coverage cutoff)") ,precision_recall_kraken %>% mutate(aligner = "kraken2", kmer = NA, genomeset = 'Unknown')) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% filter(!is.na(TP))
write.csv(aligner_summary_table,'tables/aligner_summary_table_vir.csv')

#### by clade and composition
ranks = c('Phylum','Class','Order','Family')

curvedatalist1 = list()
curvedatalist2 = list()
for(num in c(21)){
  for(r in ranks){
    print(r)
    ofinterest = bind_rows(viral_phy,viral_cla,viral_ord,viral_fam) %>% filter(type == r) %>% ungroup %>% select(lineage) %>% unlist %>% unname %>% unique
    for(f in ofinterest){
      print(f)
      a = compute_pr_auc(bind_rows(viral_phy,viral_cla,viral_ord,viral_fam)  %>% filter(kmer == num, lineage == f)) 
      a1 = a[[1]] %>% mutate(curvetype = 'P-R',rank = r,taxonomy = f,kmer = num,covtype = 'Total Coverage')
      a2 = a[[2]] %>% mutate(curvetype = 'P-R',rank = r,taxonomy = f,kmer = num,covtype = 'Unique Coverage')
      curvedatalist1[[paste(num,f,r)]] = bind_rows(a1,a2)
    }
  }
}
curvedata_fam = bind_rows(curvedatalist1)

optimize_cutoff2 <- function(data) {
  data  %>% group_by(curvetype, kmer, covtype,rank,taxonomy) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% arrange(desc(F1)) %>% slice(1)
}

optimized_cutoffs = optimize_cutoff2(curvedata_fam) %>% mutate(cutoffselect=TRUE) %>% filter(!is.na(F1))

curvedata_fam = inner_join(curvedata_fam,optimized_cutoffs) %>% mutate(cutoffselect = if_else(cutoffselect == TRUE,cutoffselect,FALSE))

ggplot(curvedata_fam %>% filter(kmer == 21,rank == 'Order') %>% mutate(taxonomy = fct_reorder(taxonomy, cutoff * (covtype == "Total Coverage"), .desc = FALSE)), aes(y = cutoff, x = taxonomy, color = covtype, alpha = F1)) + geom_point(size=3) + theme_cowplot() + theme(axis.text.x = element_text(angle = 60,hjust=1))
ggplot(curvedata_fam %>% filter(kmer == 21,rank == 'Class') %>% mutate(taxonomy = fct_reorder(taxonomy, cutoff * (covtype == "Total Coverage"), .desc = FALSE)), aes(y = cutoff, x = taxonomy, color = covtype, alpha = F1)) + geom_point(size=3) + theme_cowplot() + theme(axis.text.x = element_text(angle = 60,hjust=1))
ggplot(curvedata_fam %>% filter(kmer == 21,rank == 'Phylum') %>% mutate(taxonomy = fct_reorder(taxonomy, cutoff * (covtype == "Total Coverage"), .desc = FALSE)), aes(y = cutoff, x = taxonomy, color = covtype, alpha = F1)) + geom_point(size=3) + theme_cowplot() + theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave('plots/clade_specific_phyla.pdf',width=6,height=4)

#### actual vs expected

foo_rep = bind_rows(parsed_tp %>% filter(kmer == 21, type == 'All Genomes', genomeset %in% c('CH','CHM','CHML')) %>% mutate(kmer = paste('k =',kmer),aligner = 'Xtree')) %>% 
  select(sample,type,genomeset,expected_ra,observed_ra,tp_fp,kmer) %>% 
  filter(tp_fp == "TRUE POSITIVE") %>% 
  mutate(plot_type = 'Rep') %>%
  select(sample,plot_type,expected_ra,observed_ra,kmer)

foo_nonrep = bind_rows(parsed_tp %>% filter(kmer == 21, type == 'All Genomes', genomeset == 'PVC Non-Representative Genomes') %>% mutate(kmer = paste('k =',kmer),aligner = 'Xtree')) %>% 
  select(sample,type,genomeset,expected_ra,observed_ra,tp_fp,kmer) %>% 
  filter(tp_fp == "TRUE POSITIVE") %>% 
  mutate(plot_type = 'Non-Rep') %>%
  select(sample,plot_type,expected_ra,observed_ra,kmer)

foo_combined = bind_rows(foo_rep %>% mutate(plot_type = 'Rep'), foo_nonrep %>% mutate(plot_type = 'Non-Rep')) %>% mutate(plot_type = factor(plot_type, levels = c('Rep', 'Non-Rep'))) %>% group_by(sample, plot_type) %>% mutate(expected_ra = expected_ra/sum(expected_ra, na.rm=TRUE), observed_ra = observed_ra/sum(observed_ra, na.rm=TRUE)) %>% ungroup()
corr_labels = foo_combined %>% group_by(plot_type) %>% summarise(corr = cor(expected_ra, observed_ra, use = "pairwise.complete.obs"), .groups = 'drop') %>% mutate(label = paste0("Pearson = ", round(corr, 3)))
ggplot(foo_combined, aes(x = expected_ra, y = observed_ra)) + geom_point(alpha=0.5, color = "black") + geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") + facet_wrap(~ plot_type, nrow = 1, scales = "free") + geom_text(data = corr_labels, aes(x = 0, y = Inf, label = label), hjust = 0, vjust = 1.5, inherit.aes = FALSE) + xlim(0,0.04) + ylim(0,0.04) + labs(x = "Expected Relative Abundance", y = "Observed Relative Abundance") + theme_cowplot(font_size = 14) + theme(strip.background = element_rect(fill = "gray90", color = "gray70"), strip.text = element_text(face = "bold", size = 16), panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2), axis.text = element_text(size = 14), axis.title = element_text(size = 16))
ggsave('plots/vir_expected_observed_rep_nonrep.pdf',width=8,height=3)
