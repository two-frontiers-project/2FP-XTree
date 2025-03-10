#### THIS IS AN EXAMPLE RSCRIPT THAT LOADS IN XTREE .COV FILES, FILTERS BASED ON USER-SPECIFIED COVERAGE, COMPUTES RELATIVE ABUNDANCE (NOT LENGTH NORMALIZED), AND OUTPUTS A WIDE-FORM MATRIX
### RUN THIS SCRIPT IN THE SAME FOLDER AS WHERE YOUR .COV FILES ARE
### COVERAGE CUTOFF, SET ACCORDING TO YOUR USE CASE
COVERAGE_CUTOFF = 0.05
###

library(dplyr)
library(readr)
library(stringr)
library(reshape2)

cov_files <- list.files(".", pattern = "\\.cov$", full.names = TRUE)
cov_files = cov_files[grepl('sub',cov_files)]

process_cov_file <- function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE) %>% mutate(sample = file)%>% mutate(min_coverage = pmin(Adamantium_covered,Unique_proportion_covered,Proportion_covered)) %>% mutate(RA = Coverage_est/sum(Coverage_est))
  return(df)
}

file_dfs <- lapply(cov_files, process_cov_file)
file_dfs <- Filter(function(df) !is.null(df) && nrow(df) > 0, file_dfs)  # Remove empty dataframes

all_data <- bind_rows(file_dfs)
all_data <- all_data %>% filter(coverage >= 0.05)

### EXAMPLE, COMMENTED OUT, AS TO HOW YOU COULD MERGE IN GTDB TAXONOMY                                                                                                                                                               
taxonomy_df <- read.delim("/mnt/b/gtdb_r220/bac120_taxonomy.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(taxonomy_df) <- c("Reference", "taxonomy")
taxonomy_df$Reference = gsub('RS_','',taxonomy_df$Reference)
taxonomy_df$Reference = gsub('GB_','',taxonomy_df$Reference)
all_data <- inner_join(all_data, taxonomy_df, by = "Reference")

all_data = all_data %>% reshape2::dcast(Reference ~ sample,value.var = 'RA')
all_data[is.na(RA)] = 0                                                                                                                                                                                        
write.csv(final_data, "merged_xtree.csv")


