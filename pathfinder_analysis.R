# Library calls and data #### 
# detach("package:RUVseq", unload=TRUE)

library(pathfindR)
library(dplyr)

?read.csv
?readxl::read_excel
results_vaso_sep_shrink_df <- read.csv("differential_expression/DEG_all_vasopressors_sepsis_only.csv", header = TRUE, row.names = 1)
results_vaso_nosep_shrink_df <- read.csv("differential_expression/DEG_all_vasopressors_no_sepsis.csv", header = TRUE, row.names = 1)

# Results prepared for pathfindr #### 
sepsis_input_df <- as.data.frame(results_vaso_sep_shrink_df) %>% select(gene_symbol, log2FoldChange, padj)
sepsis_input_df  

# change NA p values to 1 
sepsis_input_df$padj <- ifelse(is.na(sepsis_input_df$padj), 1, sepsis_input_df$padj)

# Run - sepsis ####

sepsis_output_df <- run_pathfindR(sepsis_input_df, output_dir = "pathfinder_sepsis",
                                gene_sets="GO-BP", min_gset_size = 5, max_gset_size = 1000)
