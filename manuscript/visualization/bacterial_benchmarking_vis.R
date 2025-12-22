# standardized plots for true false positive rate PVC-X

library(tidyverse)
library(reshape2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(PRROC)

compute_tpr_fpr <- function(data) {
  # Initialize dataframes to store results for total and unique coverage
  total_results <- tibble(cutoff = double(), TPR = double(), FPR = double())
  unique_results <- tibble(cutoff = double(), TPR = double(), FPR = double())
  
  # Iterate over total coverage cutoffs in increments of 0.01
  for (cutoff in seq(0, 1, by = 0.01)) {

    # Filter and calculate for total coverage
    filtered_data_total <- data %>%
      filter(!is.na(sample)) %>%
      mutate(tp_fp = if_else(Proportion_covered < cutoff & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>%
      mutate(tp_fp = if_else(Proportion_covered < cutoff & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) %>%
      group_by(sample) %>%
      summarise(TP = sum(tp_fp == "TRUE POSITIVE"), 
                FP = sum(tp_fp == "FALSE POSITIVE"),
                FN = sum(tp_fp == "FALSE NEGATIVE"),
                TN = dbsize - TP - FP - FN) %>%
      mutate(TPR = TP / (TP + FN),
             FPR = FP / (FP + TN)) %>% distinct
    
    # Aggregate and store total coverage results
    aggregate_data_total <- filtered_data_total %>%
      summarise(Avg_TPR = mean(TPR, na.rm = TRUE), 
                Avg_FPR = mean(FPR, na.rm = TRUE))
    
    total_results <- total_results %>%
      add_row(cutoff = cutoff, TPR = aggregate_data_total$Avg_TPR, FPR = aggregate_data_total$Avg_FPR)
  }
  
  # Iterate over unique coverage cutoffs in increments of 0.01
  for (cutoff in seq(0, 1, by = 0.01)) {

    # Filter and calculate for unique coverage
    filtered_data_unique <- data %>%
      filter(!is.na(sample)) %>%
      mutate(tp_fp = if_else(unique_coverage_union < cutoff & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>%
      mutate(tp_fp = if_else(unique_coverage_union < cutoff & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) %>%
      group_by(sample) %>%
      summarise(TP = sum(tp_fp == "TRUE POSITIVE"), 
                FP = sum(tp_fp == "FALSE POSITIVE"),
                FN = sum(tp_fp == "FALSE NEGATIVE"),
                TN = dbsize - TP - FP - FN) %>%
      mutate(TPR = TP / (TP + FN),
             FPR = FP / (FP + TN)) %>% distinct
    
    # Aggregate and store unique coverage results
    aggregate_data_unique <- filtered_data_unique %>%
      summarise(Avg_TPR = mean(TPR, na.rm = TRUE), 
                Avg_FPR = mean(FPR, na.rm = TRUE))
    
    unique_results <- unique_results %>%
      add_row(cutoff = cutoff, TPR = aggregate_data_unique$Avg_TPR, FPR = aggregate_data_unique$Avg_FPR)
  }
  
  return(list(total_results = total_results, unique_results = unique_results))
}

compute_pr_auc <- function(data) {
  # Initialize dataframes to store results for total and unique coverage
  total_results <- tibble(cutoff = double(), Precision = double(), Recall = double(), AUC_PR = double())
  unique_results <- tibble(cutoff = double(), Precision = double(), Recall = double(), AUC_PR = double())
  
  # Iterate over total coverage cutoffs in increments of 0.01
  for (cutoff in seq(0, 1, by = 0.01)) {

    # Filter and calculate for total coverage
    filtered_data_total <- data %>%
      filter(!is.na(sample)) %>%
      mutate(tp_fp = if_else(Proportion_covered < cutoff & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>%
      mutate(tp_fp = if_else(Proportion_covered < cutoff & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) %>% 
      filter(tp_fp != 'TRUE NEGATIVE') %>%
      filter(Proportion_covered >= cutoff | tp_fp %in% c('FALSE NEGATIVE')) %>%
      group_by(sample) %>%
      summarise(TP = sum(tp_fp == "TRUE POSITIVE"), 
                FP = sum(tp_fp == "FALSE POSITIVE"),
                FN = sum(tp_fp == "FALSE NEGATIVE")) %>%
      mutate(Precision = TP / (TP + FP),
             Recall = TP / (TP + FN))
    
    # Aggregate and store total coverage results
    aggregate_data_total <- filtered_data_total %>%
      summarise(Avg_Precision = mean(Precision, na.rm = TRUE), 
                Avg_Recall = mean(Recall, na.rm = TRUE))
    
    pr_curve_total <- pr.curve(scores.class0 = filtered_data_total$FP, scores.class1 = filtered_data_total$TP)
    auc_pr_total <- pr_curve_total$auc.integral
    
    total_results <- total_results %>%
      add_row(cutoff = cutoff, Precision = aggregate_data_total$Avg_Precision, Recall = aggregate_data_total$Avg_Recall, AUC_PR = auc_pr_total)
  }
  
  # Iterate over unique coverage cutoffs in increments of 0.01
  for (cutoff in seq(0, 1, by = 0.01)) {

    # Filter and calculate for unique coverage
    filtered_data_unique <- data %>%
      filter(!is.na(sample)) %>%
      mutate(tp_fp = if_else(unique_coverage_union < cutoff & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>%
      mutate(tp_fp = if_else(unique_coverage_union < cutoff & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) %>% 
      filter(tp_fp != 'TRUE NEGATIVE') %>%
      filter(unique_coverage_union >= cutoff | tp_fp %in% c('FALSE NEGATIVE')) %>%
      group_by(sample) %>%
      summarise(TP = sum(tp_fp == "TRUE POSITIVE"), 
                FP = sum(tp_fp == "FALSE POSITIVE"),
                FN = sum(tp_fp == "FALSE NEGATIVE")) %>%
      mutate(Precision = TP / (TP + FP),
             Recall = TP / (TP + FN))
    
    # Aggregate and store unique coverage results
    aggregate_data_unique <- filtered_data_unique %>%
      summarise(Avg_Precision = mean(Precision, na.rm = TRUE), 
                Avg_Recall = mean(Recall, na.rm = TRUE))
    
    pr_curve_unique <- pr.curve(scores.class0 = filtered_data_unique$FP, scores.class1 = filtered_data_unique$TP)
    auc_pr_unique <- pr_curve_unique$auc.integral
    
    unique_results <- unique_results %>%
      add_row(cutoff = cutoff, Precision = aggregate_data_unique$Avg_Precision, Recall = aggregate_data_unique$Avg_Recall, AUC_PR = auc_pr_unique)
  }
  
  return(list(total_results = total_results, unique_results = unique_results))
}

setwd('/Users/bradentemp/Library/CloudStorage/GoogleDrive-braden@twofrontiers.org/Shared drives/The Two Frontiers Project/STUDIES/PUBLICATION_2024-11-12_xtree_pvc')

gtdb = readRDS('data_packet/comparison_data_processed/gtdbr214_benchmarking_xtree_full_dataset_auc.rds') %>% mutate(type = 'xtreeGTDB')
parsed_tp = gtdb %>% mutate(unique_coverage_union = pmin(Adamantium_prop,Unique_proportion_covered))
parsed_tp$tp_fp <- factor(parsed_tp$tp_fp, levels = c("FALSE POSITIVE", "TRUE POSITIVE","FALSE NEGATIVE","TRUE NEGATIVE"))
parsed_tp = parsed_tp %>% filter(!is.na(taxonomy))
parsed_tp$abset = gsub('highab','High Abundance',parsed_tp$abset)
parsed_tp$abset = gsub('lowab','Low Abundance',parsed_tp$abset)
parsed_tp$genomeset = gsub('NONREP','GTDB Non-Representative Species',parsed_tp$genomeset)
parsed_tp$genomeset = gsub('REP','GTDB Representative Species',parsed_tp$genomeset)

#parsed_tp = parsed_tp %>% separate(lineage, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';', fill = 'right')

metaphlan = readRDS('data_packet/comparison_data_processed/metaphlan4_benchmarking_bac_full_dataset_auc.rds')
metaphlan$abset = gsub('highab','High Abundance',metaphlan$abset)
metaphlan$abset = gsub('lowab','Low Abundance',metaphlan$abset)
metaphlan$genomeset = gsub('NONREP','GTDB Non-Representative Species',metaphlan$genomeset)
metaphlan$genomeset = gsub('REP','GTDB Representative Species',metaphlan$genomeset)

kraken2 = readRDS('data_packet/comparison_data_processed/kraken2_benchmarking_bac_full_dataset_auc_withkmers.rds')%>% filter(!is.na(conf))
kraken2$abset = gsub('highab','High Abundance',kraken2$abset)
kraken2$abset = gsub('lowab','Low Abundance',kraken2$abset)
kraken2$genomeset = gsub('NONREP','GTDB Non-Representative Species',kraken2$genomeset)
kraken2$genomeset = gsub('REP','GTDB Representative Species',kraken2$genomeset)


# cov dist

#propelt = parsed_tp %>% select(taxonomy,genomeset,abset,kmer,Proportion_covered,Unique_proportion_covered,Adamantium_prop) %>% reshape2::melt(id.vars = c('taxonomy','genomeset','abset','kmer'))
#ggplot(propelt %>% sample_n(10000),aes(x = taxonomy, y = value, fill = variable)) + facet_grid(kmer + variable ~ genomeset) + geom_bar(alpha=0.3,stat = 'identity') 

# get precision recall curves
curvedatalist1 = list()
curvedatalist2 = list()
for(num in c(17, 21,24,29)){
  for(ab in c('High Abundance','Low Abundance')){
    for(gset in c('GTDB Representative Species','GTDB Non-Representative Species')){
      for(r in c('Species')){
        a = compute_pr_auc(parsed_tp %>% filter(kmer == num,abset == ab,genomeset == gset)) 
        a1 = a[[1]] %>% mutate(curvetype = 'P-R',rank = r,kmer = num,covtype = 'Total Coverage',abset = ab,genomeset = gset)
        a2 = a[[2]] %>% mutate(curvetype = 'P-R',rank = r,kmer = num,covtype = 'Unique Coverage',abset = ab,genomeset = gset)
       # b = compute_tpr_fpr(parsed_tp %>% filter(kmer == num,rank == r,abset == ab,genomeset == gset))
      #  b1 = b[[1]] %>% mutate(curvetype = 'TPR-FPR',kmer = num,rank = r,covtype = 'Total Coverage',abset = ab,genomeset = gset)
       # b2 = b[[2]] %>% mutate(curvetype = 'TPR-FPR',kmer = num,rank = r,covtype = 'Unique Coverage',abset = ab,genomeset = gset)
        curvedatalist1[[paste(num,'a',gset,r,ab)]] = bind_rows(a1,a2)
      #  curvedatalist2[[paste(num,'b',gset,r,ab)]] = bind_rows(b1,b2)     
      }
    }
  }
}

curvedata = bind_rows(curvedatalist1) %>% filter(rank=='Species')
#curvedata2 = bind_rows(curvedatalist2) %>% filter(rank=='Species')

# plot out the curve data and pick your cutoffs

# pick cutoffs

optimize_cutoff <- function(data) {
  data %>%
    mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>%
    arrange(desc(F1)) %>%
    slice(1)
}

optimized_df <- curvedata %>% group_by(curvetype, rank, kmer, covtype, abset, genomeset) %>% do(optimize_cutoff(.))%>% select(curvetype,kmer,covtype,cutoff,rank, abset, genomeset,F1) %>% mutate(cutoffselect=TRUE)

curvedata = left_join(curvedata,optimized_df) %>% mutate(cutoffselect = if_else(cutoffselect == TRUE,cutoffselect,FALSE))

#### overall f1

curvdatatoplot = curvedata %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall))
ggplot(curvdatatoplot,aes(x = cutoff,y = F1, color = factor(kmer), shape = genomeset)) + geom_point(size=2.5,alpha=.5) + theme_cowplot() + facet_grid(covtype ~ abset) + scale_color_brewer(palette = 'Set1') + theme(legend.position = 'bottom')
ggsave('plots/f1_plot.pdf',width=8,height=4)

# full 
ggplot(data = curvedata %>% filter(covtype == 'Total Coverage'), aes(x = Precision, y = Recall, color = cutoff)) + theme_cowplot() + geom_point() + facet_grid(kmer ~ abset + genomeset)  + scale_color_viridis_c(option = 'D') + geom_text(data = filter(curvedata %>% filter(covtype == 'Total Coverage'), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(covtype == 'Total Coverage'), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50")
#ggsave('~/HMS Dropbox/Braden Tierney/viraldb/plots/total_coverage_cutoff.pdf',width=8,height=8)

ggplot(data = curvedata %>% filter(covtype == 'Total Coverage',kmer == 29), aes(x = Precision, y = Recall, color = cutoff)) + theme_cowplot() + geom_point() + facet_grid(kmer ~ abset + genomeset)  + scale_color_viridis_c(option = 'D') + geom_text(data = filter(curvedata %>% filter(covtype == 'Total Coverage',kmer == 29), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(covtype == 'Total Coverage',kmer == 29), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50")
#ggsave('~/HMS Dropbox/Braden Tierney/viraldb/plots/k29_total_coverage_cutoff.pdf',width=8,height=4)

ggplot(data = curvedata %>% filter(covtype == 'Unique Coverage'), aes(x = Precision, y = Recall, color = cutoff)) + theme_cowplot() + geom_point() + facet_grid(kmer ~ abset + genomeset)  + scale_color_viridis_c(option = 'D') + geom_text(data = filter(curvedata %>% filter(covtype == 'Unique Coverage'), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(covtype == 'Unique Coverage'), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50")
#ggsave('~/HMS Dropbox/Braden Tierney/viraldb/plots/unique_coverage_cutoff.pdf',width=8,height=8)

ggplot(data = curvedata %>% ungroup %>% filter(covtype == 'Unique Coverage',kmer == 29), aes(x = Precision, y = Recall, color = cutoff)) + theme_cowplot() + geom_point() + facet_grid(kmer ~ abset + genomeset)  + scale_color_viridis_c(option = 'D') + geom_text(data = filter(curvedata %>% filter(covtype == 'Unique Coverage',kmer == 29), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(covtype == 'Unique Coverage',kmer == 29), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50")
#ggsave('~/HMS Dropbox/Braden Tierney/viraldb/plots/k29_unique_coverage_cutoff.pdf',width=8,height=4)

ggplot() + theme_cowplot() + geom_point(alpha=0.5,data = curvedata %>% ungroup %>% filter(genomeset == 'GTDB Non-Representative Species',covtype == 'Total Coverage',kmer == 29), aes(x = Precision, y = Recall, color = cutoff)) + facet_grid(kmer ~ abset)  + scale_color_viridis_c(option = 'E') + geom_text(data = filter(curvedata %>% filter(genomeset == 'GTDB Non-Representative Species',covtype == 'Total Coverage',kmer == 29), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.5), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(genomeset == 'GTDB Non-Representative Species',covtype == 'Total Coverage',kmer == 29), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.5), colour = "grey50")+ ggnewscale::new_scale_color()+ geom_point(alpha=0.5,data = curvedata %>% ungroup %>% filter(genomeset == 'GTDB Non-Representative Species',covtype == 'Unique Coverage',kmer == 29), aes(x = Precision, y = Recall, color = cutoff)) + facet_grid(. ~ abset)  + scale_color_viridis_c(option = 'G') + geom_text(data = filter(curvedata %>% filter(genomeset == 'GTDB Non-Representative Species',covtype == 'Unique Coverage',kmer == 29), cutoffselect == TRUE), aes(label = cutoff, x = 0.25, y = 0.25), vjust = 0.75, hjust = 0.5) + geom_segment(data = filter(curvedata%>% filter(genomeset == 'GTDB Non-Representative Species',covtype == 'Unique Coverage',kmer == 29), cutoffselect == TRUE), aes(xend = Precision, yend = Recall, x = 0.25, y = 0.25), colour = "grey50") + theme(legend.position = 'bottom')
ggsave('plots/k29_both_coverages_nonrep.pdf',width=8,height=4)

# merge in optimized cutoff data 
parsed_tp_op = left_join(parsed_tp,optimized_df %>% dcast(curvetype + kmer +genomeset + rank + abset ~ covtype,value.var = 'cutoff' ) %>% mutate(kmer = as.factor(kmer)))

# plot raw data on coverage cutoffs
ggplot(parsed_tp_op %>% filter(tp_fp != 'FALSE NEGATIVE',kmer!=29) %>% filter(rank=='Species')%>% filter(!is.na(kmer)) %>% group_by(kmer) ,aes(x = Proportion_covered, y = unique_coverage_union)) + geom_point(aes(color= tp_fp),alpha=.5,size=.6)  + theme_cowplot() + geom_hline(aes(yintercept = `Unique Coverage`),color = 'red',linetype = 'dashed')+geom_vline(aes(xintercept = `Total Coverage`),color = 'red',linetype = 'dashed') + facet_grid(kmer ~ genomeset + abset,scales = 'free')  + scale_fill_brewer(palette = 'Dark2')+ scale_color_brewer(palette = 'Dark2') + ylab('Unique breadth of coverage') + xlab('Total breadth of coverage') + xlim(0,1) + ylim(0,1)
ggsave('plots/SUPP_fp_rate_by_coverage_cutoff_sp.pdf',width=12,height=9)

ggplot(parsed_tp_op %>% filter(tp_fp != 'FALSE NEGATIVE',kmer==29,genomeset == 'GTDB Non-Representative Species') %>% filter(rank=='Species')%>% filter(!is.na(kmer)) %>% group_by(kmer) ,aes(x = Proportion_covered, y = unique_coverage_union)) + geom_point(aes(color= tp_fp),alpha=.5,size=.6)  + theme_cowplot() + geom_hline(aes(yintercept = `Unique Coverage`),color = 'red',linetype = 'dashed')+geom_vline(aes(xintercept = `Total Coverage`),color = 'red',linetype = 'dashed') + facet_grid(. ~ abset,scales = 'free')  + scale_fill_manual(values = c('gray','black'))+ scale_color_manual(values = c('gray','black'))+ ylab('Unique coverage') + xlab('Total coverage') + theme(legend.position = 'bottom')+ xlim(0,1) + ylim(0,1) 
ggsave('plots/k29_fp_rate_by_coverage_cutoff_sp_nonrep.png',width=8,height=4)

# compute FP/TP rates etc by points
parsed_tp_op_filt = parsed_tp_op%>% mutate(tp_fp = if_else(unique_coverage_union<0.05 & tp_fp == 'TRUE POSITIVE','FALSE NEGATIVE',tp_fp))%>% mutate(tp_fp = if_else(unique_coverage_union<0.05 & tp_fp == 'FALSE POSITIVE','TRUE NEGATIVE',tp_fp))

precision_recall_gtdb <- parsed_tp_op %>% group_by(kmer,abset,genomeset, sample,rank) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))
precision_recall_gtdb_opt <- parsed_tp_op_filt %>% group_by(kmer,abset,genomeset, sample,rank) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))
precision_recall_kraken <- kraken2 %>% filter(rank!='Genus') %>% group_by(genomeset, abset,sample,rank,conf) %>% mutate(tp_fp = if_else(approx_cov < 0.01 & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))
precision_recall_kraken_1perc <- kraken2 %>% filter(rank!='Genus') %>% group_by(genomeset, abset,sample,rank,conf) %>% mutate(tp_fp = if_else(approx_cov < 0.01 & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>%  mutate(tp_fp = if_else(approx_cov < 0.01 & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))
precision_recall_kraken_5perc <- kraken2 %>% filter(rank!='Genus') %>% group_by(genomeset, abset,sample,rank,conf) %>% mutate(tp_fp = if_else(approx_cov < 0.05 & tp_fp == 'TRUE POSITIVE', 'FALSE NEGATIVE', tp_fp)) %>%  mutate(tp_fp = if_else(approx_cov < 0.05 & tp_fp == 'FALSE POSITIVE', 'TRUE NEGATIVE', tp_fp)) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))
precision_recall_metaph <- metaphlan%>% filter(rank!='Genus') %>% group_by(genomeset, abset,sample,rank) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))

#long_df <- tidyr::pivot_longer(precision_recall_gtdb, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "xtree (no cutoffs)")
long_df_opt <- tidyr::pivot_longer(precision_recall_gtdb_opt, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "xtree (5% coverage)")
long_df_krak <- tidyr::pivot_longer(precision_recall_kraken, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "kraken2")
#long_df_krak1 <- tidyr::pivot_longer(precision_recall_kraken_1perc, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "kraken2 (1% approx coverage)")
long_df_krak5 <- tidyr::pivot_longer(precision_recall_kraken_5perc, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "kraken2 (5% coverage)")
long_df_metaph <- tidyr::pivot_longer(precision_recall_metaph, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "MetaPhlAn4")

combineddata= bind_rows(long_df_opt, long_df_krak, long_df_krak5, long_df_metaph) 
combineddata$aligner = factor(combineddata$aligner, levels = c('xtree (5% coverage)','kraken2','kraken2 (5% coverage)','MetaPhlAn4'))

ggplot(data = combineddata %>% filter(rank != 'Genus', genomeset == 'GTDB Non-Representative Species'), aes(x = metric, y = value, fill = aligner, color = aligner, group = paste(metric, aligner))) + theme_cowplot() + geom_violin(alpha = 0.2, scale = "width", position = position_dodge(width = 0.8)) + ggbeeswarm::geom_quasirandom(width = 0.2, dodge.width = 0.8, alpha=0.5) + stat_summary(fun = median, geom = "crossbar", width = 1, position = position_dodge(width = 0.8), aes(color = aligner)) + facet_grid(. ~ abset) + ylim(0, 1) + scale_fill_manual(values = c("xtree (5% coverage)" = "#00008B", "kraken2" = "#8B0000", "kraken2 (5% coverage)" = "#E66100", "MetaPhlAn4" = "#228B22")) + scale_color_manual(values = c("xtree (5% coverage)" = "#00008B", "kraken2" = "#8B0000", "kraken2 (5% coverage)" = "#E66100", "MetaPhlAn4" = "#228B22")) + theme(legend.position = 'bottom')
ggsave('plots/precision_recall_by_sample_nonrep.pdf',width=8,height=4)

### FOR PUBLICATION
aligner_summary_table = bind_rows(precision_recall_gtdb_opt %>% mutate(aligner = "xtree (5% coverage cutoff)") ,precision_recall_kraken %>% mutate(aligner = "kraken2") ,precision_recall_kraken_5perc  %>% mutate(aligner = "kraken2 (5% coverage cutoff)"),precision_recall_metaph %>% mutate(aligner = "MetaPhlAn4")) %>% group_by(aligner,kmer,abset,genomeset,rank) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% filter(!is.na(TP))
write.csv(aligner_summary_table,'tables/aligner_summary_table_bac.csv')

aligner_summary_table = bind_rows(precision_recall_gtdb_opt %>% mutate(aligner = "xtree (5% coverage cutoff)") ,precision_recall_kraken %>% mutate(aligner = "kraken2") ,precision_recall_kraken_5perc  %>% mutate(aligner = "kraken2 (5% coverage cutoff)"),precision_recall_metaph %>% mutate(aligner = "MetaPhlAn4")) %>% group_by(aligner,kmer,abset,genomeset,rank) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% filter(!is.na(TP))%>% summarise(mean_precision = mean(Precision),sd_precision = sd(Precision),mean_recall = mean(Recall),sd_recall = sd(Recall),mean_f1 = mean(F1),sd_f1 = sd(F1))
write.csv(aligner_summary_table,'tables/aligner_summary_table_bac_AVG.csv')

#### clade specific

gtdbtax = read.delim('bac120_taxonomy_rn.tsv',sep='\t',header=F) %>% select(V2) %>% separate(V2, into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = ";", remove = TRUE, extra = "merge") %>% distinct

parsed_tp_op_tax = left_join(parsed_tp_op,gtdbtax %>% dplyr::rename(taxonomy = Species)) %>% select(Phylum,Class,Order,Family,Genus,kmer,genomeset,abset,sample,tp_fp,Proportion_covered,unique_coverage_union,expected_ra,observed_ra) %>% filter(genomeset == 'GTDB Non-Representative Species')#%>% group_by(Family,kmer,genomeset,abset,sample) 
#parsed_tp_op_tax = parsed_tp_op_tax %>% mutate(tp_fp = if_else(any(tp_fp == "TRUE POSITIVE"), "TRUE POSITIVE",if_else(any(tp_fp == "FALSE POSITIVE"), "FALSE POSITIVE",if_else(all(tp_fp %in% c("FALSE NEGATIVE", "TRUE NEGATIVE")), "FALSE NEGATIVE", first(tp_fp))))) %>% ungroup()

#parsed_tp_op_tax = parsed_tp_op_tax %>% group_by(Family,kmer,genomeset,abset,sample,tp_fp) %>% summarise(Proportion_covered = mean(Proportion_covered),unique_coverage_union = mean(unique_coverage_union),expected_ra = mean(expected_ra),expected_ra = mean(expected_ra))

# get precision recall curves

ranks = c('Phylum','Class','Order')

curvedatalist1 = list()
curvedatalist2 = list()
for(num in c(21,24,29)){
  for(r in ranks){
    ofinterest = parsed_tp_op_tax %>% filter(!is.na(r)) %>% select(r) %>% unlist %>% unname %>% unique
      for(f in ofinterest){
        print(f)
        a = compute_pr_auc(parsed_tp_op_tax %>% filter(kmer == num, get(r) == f)) 
        a1 = a[[1]] %>% mutate(curvetype = 'P-R',rank = r,taxonomy = f,kmer = num,covtype = 'Total Coverage')
        a2 = a[[2]] %>% mutate(curvetype = 'P-R',rank = r,taxonomy = f,kmer = num,covtype = 'Unique Coverage')
        curvedatalist1[[paste(num,f,r)]] = bind_rows(a1,a2)
      }
  }
}
saveRDS(curvedatalist1,'./tmp_rankcurvedata.rds')
curvedatalist1 = readRDS('tmp_rankcurvedata.rds')
curvedata_fam = bind_rows(curvedatalist1)

optimize_cutoff2 <- function(data) {
  data  %>% group_by(curvetype, kmer, covtype,rank,taxonomy) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% arrange(desc(F1)) %>% slice(1)
}

optimized_cutoffs = optimize_cutoff2(curvedata_fam) %>% mutate(cutoffselect=TRUE) %>% filter(!is.na(F1))

curvedata_fam = inner_join(curvedata_fam,optimized_cutoffs) %>% mutate(cutoffselect = if_else(cutoffselect == TRUE,cutoffselect,FALSE))

ggplot(curvedata_fam %>% filter(kmer == 29,rank == 'Order') %>% mutate(taxonomy = fct_reorder(gsub('o__','',taxonomy), cutoff * (covtype == "Total Coverage"), .desc = FALSE)), aes(y = cutoff, x = taxonomy, shape = covtype, size = F1)) + geom_point(alpha=.8) + theme_cowplot() + theme(axis.text.x = element_text(angle = 60,hjust=1)) + scale_size_continuous(range = c(0.05,5)) 
ggsave('plots/clade_specific_order.pdf',width=16,height=4)

ggplot(curvedata_fam %>% filter(kmer == 29,rank == 'Class') %>% mutate(taxonomy = fct_reorder(taxonomy, cutoff * (covtype == "Total Coverage"), .desc = FALSE)), aes(y = cutoff, x = taxonomy, color = covtype, alpha = F1)) + geom_point(size=3) + theme_cowplot() + theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave('plots/clade_specific_class.pdf',width=8,height=5)

ggplot(curvedata_fam %>% filter(kmer == 29,rank == 'Phylum') %>% mutate(taxonomy = fct_reorder(taxonomy, cutoff * (covtype == "Total Coverage"), .desc = FALSE)), aes(y = cutoff, x = taxonomy, color = covtype, alpha = F1)) + geom_point(size=3) + theme_cowplot() + theme(axis.text.x = element_text(angle = 60,hjust=1))
ggsave('plots/clade_specific_phyla.pdf',width=6,height=4)


##### EXPECTED VS OBSERVED ABUNDANCE

## REP
foo = bind_rows(kraken2 %>% mutate(type = 'kraken2',kmer = 'default kraken2 kmer'),metaphlan %>% mutate(kmer = 'not kmer-based',observed_ra = observed_ra/100),parsed_tp %>% mutate(kmer = paste('k =',kmer),type = 'Xtree') )%>% select(sample,type,abset,genomeset,expected_ra,observed_ra,tp_fp,kmer) %>% filter(tp_fp == "TRUE POSITIVE") %>% group_by(sample,genomeset,type,kmer,abset) %>% mutate(expected_ra = expected_ra/sum(expected_ra, na.rm=TRUE), observed_ra = observed_ra/sum(observed_ra, na.rm=TRUE))

ggplot(foo, aes(x = expected_ra, y = observed_ra)) + geom_point(alpha=0.5, color = "black") + geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "black") + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") + facet_wrap(abset ~ type + kmer + genomeset, scales = 'free') + labs(x = "Expected Relative Abundance", y = "Observed Relative Abundance") + theme_cowplot(font_size = 14) + theme(strip.background = element_rect(fill = "gray90", color = "gray70"), strip.text = element_text(face = "bold", size = 16), panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2), axis.text = element_text(size = 14), axis.title = element_text(size = 16))
ggsave('./plots/SUPP_expected_observed_full.pdf',width=16,height=16)

## NONREP
foo = bind_rows(kraken2 %>% mutate(type = 'kraken2'),metaphlan %>% mutate(observed_ra = observed_ra/100),parsed_tp %>% mutate(type = 'Xtree')  %>% filter(kmer == 29)) %>% select(sample,type,abset,genomeset,expected_ra,observed_ra,tp_fp,kmer) %>% filter(tp_fp == "TRUE POSITIVE",genomeset == 'GTDB Non-Representative Species') %>% group_by(sample,genomeset,type,kmer,abset) %>% mutate(expected_ra = expected_ra/sum(expected_ra, na.rm=TRUE), observed_ra = observed_ra/sum(observed_ra, na.rm=TRUE))

corr_labels = foo %>% group_by(abset, type) %>% summarise(corr = cor(expected_ra, observed_ra, use = "pairwise.complete.obs"), pval = {mod = lm(observed_ra ~ expected_ra, data = cur_data()); summary(mod)$coefficients[2,4]}, .groups = 'drop') %>% mutate(label = paste0("Pearson = ", round(corr, 2), ", p<", formatC(pval, format = "e", digits = 1)))

axis_ranges = foo %>% group_by(abset, type) %>% summarise(x_min = min(expected_ra, na.rm = TRUE), x_max = max(expected_ra, na.rm = TRUE), y_min = min(observed_ra, na.rm = TRUE), y_max = max(observed_ra, na.rm = TRUE), .groups = 'drop') %>% mutate(range_min = pmin(x_min, y_min), range_max = pmax(x_max, y_max), range_min = range_min - (range_max - range_min) * 0.05, range_max = range_max + (range_max - range_min) * 0.05)

ggplot(foo, aes(x = expected_ra, y = observed_ra)) + geom_point(alpha=0.5, color = "black") + geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = "red", aes(group = 1)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "lightgreen", linewidth = 1) + facet_wrap(abset ~ type, scales = 'free', nrow = 1) + geom_text(data = corr_labels, aes(x = 0, y = Inf, label = label), hjust = 0, vjust = 1.5, inherit.aes = FALSE, fontface = "bold", color = "black") + geom_blank(data = axis_ranges %>% mutate(expected_ra = range_min, observed_ra = range_min), inherit.aes = FALSE, aes(x = expected_ra, y = observed_ra)) + geom_blank(data = axis_ranges %>% mutate(expected_ra = range_max, observed_ra = range_max), inherit.aes = FALSE, aes(x = expected_ra, y = observed_ra)) + labs(x = "Expected Relative Abundance", y = "Observed Relative Abundance") + theme_cowplot(font_size = 14) + theme(strip.background = element_rect(fill = "gray90", color = "gray70"), strip.text = element_text(face = "bold", size = 16), panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2), axis.text = element_text(size = 14), axis.title = element_text(size = 16))
ggsave('../../plots/expected_observed_nonrep.pdf',width=16,height=3)



