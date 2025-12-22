setwd('/Users/bradentemp/Library/CloudStorage/GoogleDrive-braden@twofrontiers.org/Shared drives/The Two Frontiers Project/STUDIES/PUBLICATION_2024-11-12_xtree_pvc/data_packet/euk_bench')

library(tidyverse)
library(cowplot)
library(ggbeeswarm) 

SINGLE_COLOR <- "black"

load_summary_metrics <- function(base_dir, data_type) {
  read_tsv(file.path(base_dir, "metagenome_summary_metrics.tsv"), show_col_types = FALSE) %>%
    mutate(Data_Type = data_type)
}

load_per_genome_results <- function(base_dir, data_type) {
  read_tsv(file.path(base_dir, "per_genome_classification_results.tsv"), show_col_types = FALSE) %>%
    mutate(Data_Type = data_type)
}

plot_precision_recall_distribution <- function(summary_data, target_cutoff) {
  plot_data <- summary_data %>%
    filter(Cutoff == target_cutoff) %>%
    pivot_longer(
      cols = c(Sensitivity, Specificity),
      names_to = "Metric",
      values_to = "Value"
    ) %>%
    mutate(Metric = factor(Metric, 
                           levels = c("Specificity", "Sensitivity"),
                           labels = c("Precision", "Recall")))
  
  ggplot(plot_data, aes(x = Metric, y = Value)) +
    geom_violin(fill = SINGLE_COLOR, alpha = 0.2, draw_quantiles = FALSE) +
    geom_quasirandom(
      color = SINGLE_COLOR, 
      size = 2, 
      alpha = 0.8
    ) +
    facet_wrap(~ Data_Type, ncol = 2) +
    ylim(0, 1) +
    labs(
      y = "Value (0-1.0)",
      x = NULL
    ) +
    theme_cowplot(font_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold"),
      panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2)
    )
}


plot_abundance_correlation <- function(genome_data) {
  plot_data <- genome_data %>%
    filter(Classification == "True Positive (TP)") %>%
    mutate(
      Predicted_Relative_Abundance = pmax(Predicted_Relative_Abundance, 0),
      Actual_Coverage = pmax(Actual_Coverage, 0)
    ) %>%
    group_by(Metagenome_Key, Data_Type) %>%
    mutate(
      Predicted_Relative_Abundance = Predicted_Relative_Abundance/sum(Predicted_Relative_Abundance, na.rm=TRUE),
      Actual_Coverage = Actual_Coverage/sum(Actual_Coverage, na.rm=TRUE)
    ) %>%
    ungroup()
  
  corr_labels <- plot_data %>%
    group_by(Data_Type) %>%
    summarise(
      corr = cor(Predicted_Relative_Abundance, Actual_Coverage, use = "pairwise.complete.obs"),
      label = paste0("Pearson = ", round(corr, 3)),
      .groups = 'drop'
    )
  
  ggplot(plot_data, aes(x = Predicted_Relative_Abundance, y = Actual_Coverage)) +
    geom_point(color = SINGLE_COLOR, alpha = 0.5) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") + 
    geom_smooth(method = "lm", se = FALSE, linewidth = 1, color = SINGLE_COLOR) +
    facet_wrap(~ Data_Type, nrow = 1, scales = "free") +
    geom_text(data = corr_labels, 
              aes(x = 0, y = Inf, label = label),
              hjust = 0, vjust = 1.5,
              inherit.aes = FALSE,
              fontface = "bold") +
    labs(
      x = "Expected Relative Abundance",
      y = "Observed Relative Abundance"
    ) +
    theme_cowplot(font_size = 14) +
    theme(
      strip.background = element_rect(fill = "gray90", color = "gray70"),
      strip.text = element_text(face = "bold", size = 16),
      legend.position = "none",
      panel.grid.major.y = element_line(color = "gray80", linewidth = 0.2),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 16)
    )
}


TARGET_CUTOFF <- "0.01"
fungi_summary <- load_summary_metrics("./fungi_syndata", "Fungi")
fungi_genome <- load_per_genome_results("./fungi_syndata", "Fungi")
protozoa_summary <- load_summary_metrics("./protozoa_syndata", "Protozoa")
protozoa_genome <- load_per_genome_results("./protozoa_syndata", "Protozoa")
combined_summary <- bind_rows(fungi_summary, protozoa_summary)
combined_genome <- bind_rows(fungi_genome, protozoa_genome)

aligner_summary_table_euk <- combined_summary %>% filter(Cutoff == TARGET_CUTOFF) %>% mutate(aligner = "xtree (1% coverage cutoff)", F1 = 2 * (Specificity * Sensitivity) / (Specificity + Sensitivity)) %>% dplyr::rename(Precision = Specificity, Recall = Sensitivity, sample = Metagenome_Key)
write.csv(aligner_summary_table_euk, '../../tables/aligner_summary_table_euk.csv', row.names = FALSE)

mean_precision_recall <- aligner_summary_table_euk %>% group_by(Data_Type) %>% summarise(mean_precision = mean(Precision, na.rm=TRUE), mean_recall = mean(Recall, na.rm=TRUE), mean_f1 = mean(F1, na.rm=TRUE), sd_precision = sd(Precision, na.rm=TRUE), sd_recall = sd(Recall, na.rm=TRUE), sd_f1 = sd(F1, na.rm=TRUE), .groups = 'drop')
print(mean_precision_recall)
write.csv(mean_precision_recall, '../../tables/aligner_summary_table_euk_AVG.csv', row.names = FALSE)

plot_precision_recall_distribution(combined_summary, TARGET_CUTOFF)
ggsave("../../plots/xtree_precision_recall_distribution_faceted_EUK.pdf", width = 8, height = 3)

plot_abundance_correlation(combined_genome)
ggsave("../../plots/xtree_abundance_correlation_faceted_TP_only_EUK.pdf", width = 8, height = 3)






