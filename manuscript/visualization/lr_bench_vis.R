# standardized plots for true false positive rate PVC-X

library(tidyverse)
library(reshape2)
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

dat = readRDS('data_packet/lrdat/lr_benchmark.rds') %>% dplyr::rename(expected_ra = expected_ab)
dat = dat %>% mutate(tp_fp = 'FALSE POSITIVE',tp_fp = if_else(expected_ra>0 & observed_ra>0,'TRUE POSITIVE',if_else(expected_ra==0 & observed_ra>0,'FALSE POSITIVE',if_else(expected_ra>0 &observed_ra == 0,'FALSE NEGATIVE',tp_fp)))) %>% mutate(unique_coverage_union = pmin(Adamantium_prop,Unique_proportion_covered))
dat = dat %>% mutate(technology = if_else(grepl('GridION',sample),'GridION',if_else(grepl('PromethION',sample),'PromethION','Hi-Fi')))

dat_tech = dat %>% select(sample,tp_fp,Proportion_covered,Unique_proportion_covered,Adamantium_prop) %>% melt()
dat_tech$technology = factor('Hi-Fi','PromethION','GridION')

ggplot(dat_tech,aes(x = value,y=sample,shape=tp_fp,color=variable)) + theme_cowplot() + geom_violin(alpha = 0.2, scale = "width", position = position_dodge(width = 0.8)) + ggbeeswarm::geom_quasirandom(width = 0.2, dodge.width = 0.8, alpha=0.5)

# compute FP/TP rates etc by points
dat2 = dat %>% mutate(tp_fp = if_else(unique_coverage_union<0.1 & tp_fp == 'TRUE POSITIVE','FALSE NEGATIVE',tp_fp))%>% mutate(tp_fp = if_else(unique_coverage_union<0.1 & tp_fp == 'FALSE POSITIVE','TRUE NEGATIVE',tp_fp))

ggplot(dat2 %>% filter(tp_fp %in% c('TRUE POSITIVE', 'FALSE NEGATIVE')), aes(y = observed_ra, x = expected_ra / 100, color = technology, shape = tp_fp)) + geom_point(size = 3, aes(alpha = unique_coverage_union)) + theme_cowplot() + facet_wrap(. ~ microbeset, nrow = 2) + theme(legend.position = 'bottom') + geom_abline(slope = 1, intercept = 0, color = "black", linetype = "dashed", size = 1, alpha = 0.3) + xlab('Expected Abundance') + ylab('Observed Abundance') + ylim(0,1) + xlim(0,1)+ scale_y_log10() + scale_x_log10() + scale_alpha_continuous(range = c(0.5, 1))  + scale_color_viridis_d()
ggsave('plots/lr_admt_exp_obs.pdf',width=6,height=3)

precision_recall_dat <- dat2 %>% group_by(sample,microbeset) %>% summarise(TP = sum(tp_fp == "TRUE POSITIVE"), FP = sum(tp_fp == "FALSE POSITIVE"),FN = sum(tp_fp == "FALSE NEGATIVE")) %>% mutate(Precision = TP / (TP + FP),Recall = TP / (TP + FN))

long_df_opt <- tidyr::pivot_longer(precision_recall_dat, cols = c("Precision", "Recall"), names_to = "metric", values_to = "value") %>% mutate(aligner = "xtree (10% coverage)")

combineddata= bind_rows(long_df_opt) 
combineddata$aligner = factor(combineddata$aligner, levels = c('xtree (10% coverage)'))

ggplot(data = combineddata, aes(x = metric, y = value, fill = microbeset, color = microbeset, group = paste(metric))) + stat_summary(fun = median, geom = "crossbar", width = 1, position = position_dodge(width = 0.8),aes(color = microbeset))+ theme_cowplot() + geom_violin(alpha = 0.2, scale = "width", position = position_dodge(width = 0.8)) + ggbeeswarm::geom_quasirandom(width = 0.2, dodge.width = 0.8, alpha=1,size=2)  + ylim(0, 1) + scale_color_viridis_d(option = 'H')  + theme(legend.position = 'none')
ggsave('plots/precision_recall_by_sample_nonrep_LR.pdf',width=2,height=3)

### FOR PUBLICATION
aligner_summary_table = bind_rows(precision_recall_dat %>% mutate(aligner = "xtree (10% cutoff)")) %>% group_by(microbeset,sample) %>% mutate(F1 = 2 * (Precision * Recall) / (Precision + Recall)) %>% filter(!is.na(TP))
write.csv(aligner_summary_table,'tables/aligner_summary_table_bac_LR.csv')

aligner_summary_table_avg = aligner_summary_table %>% group_by(aligner, microbeset) %>% summarise(mean_precision = mean(Precision, na.rm=TRUE), sd_precision = sd(Precision, na.rm=TRUE), mean_recall = mean(Recall, na.rm=TRUE), sd_recall = sd(Recall, na.rm=TRUE), mean_f1 = mean(F1, na.rm=TRUE), sd_f1 = sd(F1, na.rm=TRUE), .groups = 'drop')
write.csv(aligner_summary_table_avg,'tables/aligner_summary_table_bac_LR_AVG.csv')
print("Long Read F1 scores by microbeset:")
print(aligner_summary_table_avg)

aligner_summary_table_overall = aligner_summary_table %>% group_by(aligner) %>% summarise(mean_precision = mean(Precision, na.rm=TRUE), sd_precision = sd(Precision, na.rm=TRUE), mean_recall = mean(Recall, na.rm=TRUE), sd_recall = sd(Recall, na.rm=TRUE), mean_f1 = mean(F1, na.rm=TRUE), sd_f1 = sd(F1, na.rm=TRUE), .groups = 'drop')
write.csv(aligner_summary_table_overall,'tables/aligner_summary_table_bac_LR_OVERALL.csv')
print("Long Read F1 scores overall:")
print(aligner_summary_table_overall)
