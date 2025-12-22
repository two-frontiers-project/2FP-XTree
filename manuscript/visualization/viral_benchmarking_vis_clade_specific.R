# standardized plots for true false positive rate PVC-X

library(tidyverse)
library(reshape2)
library(cowplot)
library(RColorBrewer)
library(ggplot2)
library(PRROC)


### FUNCTIONS
library(tidyverse)
library(PRROC)

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
      filter(Proportion_covered >= cutoff | tp_fp %in% c('TRUE NEGATIVE', 'FALSE NEGATIVE')) %>%
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
      filter(unique_coverage_union >= cutoff | tp_fp %in% c('TRUE NEGATIVE', 'FALSE NEGATIVE')) %>%
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

setwd('~/Dropbox (Mason Lab)/viraldb/')

parsed_tp = readRDS('viral_benchmarking_xtree_clade_optimized.rds')
parsed_tp = parsed_tp %>% mutate(unique_coverage_union = pmin(Adamantium_prop,Unique_proportion_covered))
parsed_tp$tp_fp <- factor(parsed_tp$tp_fp, levels = c("FALSE POSITIVE", "TRUE POSITIVE","FALSE NEGATIVE","TRUE NEGATIVE"))
#parsed_tp = parsed_tp %>% separate(lineage, into = c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep = ';', fill = 'right')
###### optimal coverage cutoff by clade by database

# get precision recall curves
curvedatalist1 = list()
curvedatalist2 = list()
for(num in c(24)){
  for(db in c('CHM')){
    for(p in unique(parsed_tp %>% filter(!is.na(sample)) %>% select(sample) %>% unlist %>% unname)){
      print(p)
      a = compute_pr_auc(parsed_tp %>% filter(kmer == num,dbqual == db,sample == p)) 
      a1 = a[[1]] %>% mutate(curvetype = 'P-R',dbqual = db,kmer = num,covtype = 'Total Coverage',Lineage = p)
      a2 = a[[2]] %>% mutate(curvetype = 'P-R',dbqual = db,kmer = num,covtype = 'Unique Coverage',Lineage = p)
      #   b = compute_tpr_fpr(parsed_tp %>% filter(kmer == num,dbqual == db,Phylum == p))
      #    b1 = b[[1]] %>% mutate(curvetype = 'TPR-FPR',dbqual = db,kmer = num,covtype = 'Total Coverage',phylum = p)
      #   b2 = b[[2]] %>% mutate(curvetype = 'TPR-FPR',dbqual = db,kmer = num,covtype = 'Unique Coverage',phylum = p)
      curvedatalist1[[paste(num,db,'a',p)]] = bind_rows(a1,a2)
      # curvedatalist2[[paste(num,db,'b',p)]] = bind_rows(b1,b2)
    }
  }
}

curvedata_phy = bind_rows(curvedatalist1)
#curvedata2_phy = bind_rows(curvedatalist2)


optimize_cutoff <- function(data) {
  data %>%
    mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>%
    arrange(desc(F1)) %>%
    slice(1)
}

optimized_df_phy <- curvedata_phy %>% group_by(curvetype, dbqual, kmer, covtype,Lineage) %>% do(optimize_cutoff(.))%>% select(curvetype,dbqual,kmer,covtype,cutoff,Lineage) %>% mutate(cutoffselect=TRUE)

optimized_df_old =optimized_df
colnames(optimized_df_old)[6] = paste0('Overall_',colnames(optimized_df_old)[6]) 
optimized_df_old = left_join(curvedata_phy2 %>% select(Class,covtype,kmer,curvetype,dbqual) %>% distinct,optimized_df_old)

curvedata_phy2 = left_join(curvedata_phy,optimized_df_phy) %>% mutate(cutoffselect = if_else(cutoffselect == TRUE,cutoffselect,FALSE))
curvedata_phy3 = left_join(curvedata_phy2,optimized_df_old) %>% mutate(Overall_cutoffselect = if_else(Overall_cutoffselect == TRUE,Overall_cutoffselect,FALSE))
























### ACTUAL/EXPECTED ABUNDANCE

