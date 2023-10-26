# Library calls #### 
# library(gProfileR)
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(pheatmap)
library(Cairo)
# None of my custom functions used in this file 

# Data #### 

results_vaso_shrink_DEG 

results_vaso_sep_shrink_DEG

# results_vaso_shrink_DEG <- read.csv("DEG_sig_vasopressors.csv")
# results_vaso_sep_shrink_DEG <- read.csv("DEG_sig_vasopressors_sepsis_only.csv")


# Pathway analysis ####

# ** Vasopressors, all groups ####
DEG_gene_names <- rownames(results_vaso_shrink_DEG)
DEG_gene_names

pathways_vaso <- gost(DEG_gene_names,
                          organism = "hsapiens", ordered_query = TRUE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = TRUE, 
                          user_threshold = 0.05, correction_method = "g_SCS", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC","KEGG","REAC"), as_short_link = FALSE)

writexl::write_xlsx(list(sheet1 = pathways_vaso$result[,!(names(pathways_vaso$result) %in% c("parents", "evidence_codes", "intersection"))]), 
                    "pathways_vasopressors_allgroups.xlsx")

head(pathways_vaso$result)
nrow(pathways_vaso$result %>% filter(precision > 0.35)) # 35 
pathways_vaso$result %>% filter(precision > 0.35) %>% dplyr::select(source, term_name, p_value)
# 19 GO:BP 

# ** Vasopressors, sepsis only ####

DEG_gene_names_sep <- rownames(results_vaso_sep_shrink_DEG)
DEG_gene_names_sep

pathways_vaso_sep <- gost(DEG_gene_names_sep,
                      organism = "hsapiens", ordered_query = TRUE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC","KEGG","REAC"), as_short_link = FALSE)

nrow(pathways_vaso_sep$result) # 357 

writexl::write_xlsx(list(sheet1 = pathways_vaso_sep$result[,!(names(pathways_vaso_sep$result) %in% c("parents", "evidence_codes", "intersection"))]), 
                    "pathways_vasopressors_sepsis_only.xlsx")

pathways_vaso_sep2 <- apply(pathways_vaso_sep$result[,!(names(pathways_vaso_sep$result) %in% c("parents", "evidence_codes"))],2,as.character)
write.csv(pathways_vaso_sep2, "pathways_vasopressors_sepsis_only.csv")

head(pathways_vaso_sep$result)
nrow(pathways_vaso_sep$result %>% filter(precision > 0.35)) # 30 
pathways_vaso_sep$result %>% filter(precision > 0.35) %>% dplyr::select(source, term_name, p_value)
# 14 GO:BP 

# ** Vasopressors, no sepsis only, top 1000 genes ####

results_vaso_nosep_shrink
DEG_gene_names_nosep <- as.data.frame(results_vaso_nosep_shrink) %>% arrange(padj) %>% dplyr::select(ensembl_id)
DEG_gene_names_nosep <- DEG_gene_names_nosep$ensembl_id[1:1000]
DEG_gene_names_nosep

pathways_vaso_nosep <- gost(DEG_gene_names_nosep,
                          organism = "hsapiens", ordered_query = TRUE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = TRUE, 
                          user_threshold = 0.05, correction_method = "g_SCS", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC","KEGG","REAC"), as_short_link = FALSE)

nrow(pathways_vaso_nosep$result) # 297 

writexl::write_xlsx(list(nosepsis = pathways_vaso_nosep$result[,!(names(pathways_vaso_sep$result) %in% c("parents", "evidence_codes", "intersection"))]), 
                    "pathways_vasopressors_nosepsis_top1000.xlsx")

pathways_vaso_nosep2 <- apply(pathways_vaso_nosep$result[,!(names(pathways_vaso_nosep$result) %in% c("parents", "evidence_codes"))],2,as.character)
write.csv(pathways_vaso_nosep2, "pathways_vasopressors_nosepsis_top1000.csv")

head(pathways_vaso_nosep$result)
nrow(pathways_vaso_nosep$result %>% filter(precision > 0.35)) # 30 
pathways_vaso_nosep$result %>% filter(precision > 0.35) %>% dplyr::select(source, term_name, p_value)
# 14 GO:BP 

# Common pathways in 2 sets #### 

intersect(pathways_vaso_sep$result$term_id, pathways_vaso_nosep$result$term_id) # 134 pathways in common 

common_pathways <- inner_join(pathways_vaso_sep$result, pathways_vaso_nosep$result, by = "term_id", suffix = c("_vaso_sep", "_vaso_nosep")) %>% 
  dplyr::select(-contains(c("parents", "evidence_codes", "query", "significant", 
                            "effective_domain_size", "source_order")))

# Add rank 

common_pathways$rank_sepsis <- ifelse(common_pathways$source_vaso_sep %in% c("GO:BP", "REAC"), 
                                      rank(common_pathways$p_value_vaso_sep[common_pathways$source_vaso_sep %in% c("GO:BP", "REAC")]), 
                                      NA)

common_pathways$rank_nosepsis <- ifelse(common_pathways$source_vaso_nosep %in% c("GO:BP", "REAC"), 
                                        rank(common_pathways$p_value_vaso_nosep[common_pathways$source_vaso_nosep %in% c("GO:BP", "REAC")]), 
                                        NA)
common_pathways <- common_pathways %>% relocate(rank_sepsis, rank_nosepsis)

nrow(common_pathways) # 134 
write.csv(common_pathways, "common_pathways_sep_nosep.csv")

writexl::write_xlsx(list(sheet1 = common_pathways),
                    "common_pathways_sep_nosep.xlsx")

# typeof(common_pathways$precision_vaso_sep)
# common_pathways_highsig <- common_pathways %>% filter(precision_vaso_sep > 0.35 & precision_vaso_nosep > 0.35) %>% 
#   dplyr::select(contains(c("precision", "source", "term_name", "p_value"))) # 11 rows 
# 
# write.csv(common_pathways_highsig, "common_pathways_high.csv")

# Specific pathways in sepsis but not in no sepsis #### 

setdiff(pathways_vaso_sep$result$term_id, pathways_vaso_nosep$result$term_id) # 223 pathways 
357 - 134 # makes sense :D 

specific_pathways <- anti_join(pathways_vaso_sep$result, pathways_vaso_nosep$result, by = "term_id", suffix = c("_vaso_sep", "_vaso_nosep")) %>% 
  dplyr::select(-contains(c("parents", "evidence_codes", "query", "significant", 
                            "effective_domain_size", "source_order")))

nrow(specific_pathways) # 134 

write.csv(specific_pathways, "specific_pathways_sepsis_only.csv")

writexl::write_xlsx(list(sheet1 = specific_pathways), 
                    "specific_pathways_sepsis_only.xlsx")


# ** PC8 negative loadings, top 1000 genes ####

gene_names_pc8_neg_100 <- get_loadings_gene_list(pca_data = top_genes_laterpcs,
                                             pc = "PC8",
                                             direction = "negative",
                                             mart = mart, number = 100)
gene_names_pc8_neg_100

gene_names_pc8_neg_100$ensembl_gene_id

pathways_pc8_neg <- gost(gene_names_pc8_neg_100$ensembl_gene_id,
                            organism = "hsapiens", ordered_query = TRUE, 
                            multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                            measure_underrepresentation = FALSE, evcodes = TRUE, 
                            user_threshold = 0.05, correction_method = "g_SCS", 
                            domain_scope = "annotated", custom_bg = NULL, 
                            numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC","KEGG","REAC"), as_short_link = FALSE)

nrow(pathways_pc8_neg$result) # 166 

writexl::write_xlsx(list(pc8_neg = pathways_pc8_neg$result[,!(names(pathways_pc8_neg$result) %in% c("parents", "evidence_codes", "intersection"))]), 
                    "pathways_pc8loadings_top100.xlsx")

pathways_pc8_neg2 <- apply(pathways_pc8_neg$result[,!(names(pathways_vaso_nosep$result) %in% c("parents", "evidence_codes"))],2,as.character)
write.csv(pathways_pc8_neg2, "pathways_pc8loadings_top100.csv")

# ** PC15 negative loadings, top 1000 genes ####
# there are only 496 loadings in total 
# Aaaah they are the same genes just ordered differently! 
# Because they PCA is done on top 500 most variable genes 
# Do 100 

gene_names_pc15_neg_100 <- get_loadings_gene_list(pca_data = top_genes_laterpcs,
                                                  pc = "PC15",
                                                  direction = "negative",
                                                  mart = mart, number = 100)
gene_names_pc15_neg_100

gene_names_pc15_neg_100$ensembl_gene_id

pathways_pc15_neg <- gost(gene_names_pc15_neg_100$ensembl_gene_id,
                         organism = "hsapiens", ordered_query = TRUE, 
                         multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                         measure_underrepresentation = FALSE, evcodes = TRUE, 
                         user_threshold = 0.05, correction_method = "g_SCS", 
                         domain_scope = "annotated", custom_bg = NULL, 
                         numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC","KEGG","REAC"), as_short_link = FALSE)

nrow(pathways_pc15_neg$result) # 78

writexl::write_xlsx(list(pc15_neg = pathways_pc15_neg$result[,!(names(pathways_pc15_neg$result) %in% c("parents", "evidence_codes", "intersection"))]), 
                    "pathways_pc15loadings_top100.xlsx")

pathways_pc15_neg2 <- apply(pathways_pc15_neg$result[,!(names(pathways_vaso_nosep$result) %in% c("parents", "evidence_codes"))],2,as.character)
write.csv(pathways_pc15_neg2, "pathways_pc15loadings_top100.csv")

# PC8 versus PC15? 

setdiff(pathways_pc8_neg$result$term_id, pathways_pc15_neg$result$term_id) # 125 
# why are they the same pathways? 

# Save gene lists and common genes #### 

rownames(gene_names_pc8_neg_100)
gene_names_pc8_neg_100$rank_pc8 <- rownames(gene_names_pc8_neg_100)
gene_names_pc15_neg_100$rank_pc15 <- rownames(gene_names_pc15_neg_100)

inner_join(gene_names_pc8_neg_100, gene_names_pc15_neg_100)

writexl::write_xlsx(list(top100_genes_pc8 = gene_names_pc8_neg_100,
                         top100_genes_pc15 = gene_names_pc15_neg_100,
                         common_genes = inner_join(gene_names_pc8_neg_100, gene_names_pc15_neg_100)), 
                    "top100_loadings_genes_pc8_pc15_common.xlsx")

# Common pathways in 2 sets #### 

# intersect(pathways_vaso_sep$result$term_id, pathways_vaso_nosep$result$term_id) # 134 pathways in common 

common_pathways_pc <- inner_join(pathways_pc8_neg$result, 
                                 pathways_pc15_neg$result, 
                                 by = "term_id", suffix = c("_pc8", "_pc15")) %>% 
  dplyr::select(-contains(c("parents", "evidence_codes", "query", "significant", 
                            "effective_domain_size", "source_order")))

# Add rank 
View(common_pathways_pc)

common_pathways_pc$rank_pc8 <- ifelse(common_pathways_pc$source_pc8 %in% c("GO:BP", "REAC"), 
                                      rank(common_pathways_pc$p_value_pc8[common_pathways_pc$source_pc8 %in% c("GO:BP", "REAC")]), 
                                      NA)

common_pathways_pc$rank_pc15 <- ifelse(common_pathways_pc$source_pc15 %in% c("GO:BP", "REAC"), 
                                        rank(common_pathways_pc$p_value_pc15[common_pathways_pc$source_pc15 %in% c("GO:BP", "REAC")]), 
                                        NA)
common_pathways_pc <- common_pathways_pc %>% relocate(rank_pc8, rank_pc15)

nrow(common_pathways_pc) # 41
write.csv(common_pathways_pc, "common_pathways_pc8_pc15.csv")

writexl::write_xlsx(list(common = common_pathways_pc,
                         specific_pc8 = specific_pathways_pc8,
                         specific_pc15 = specific_pathways_pc15),
                    "common_and_specific_pathways_pc8_pc15.xlsx")

# typeof(common_pathways$precision_vaso_sep)
# common_pathways_highsig <- common_pathways %>% filter(precision_vaso_sep > 0.35 & precision_vaso_nosep > 0.35) %>% 
#   dplyr::select(contains(c("precision", "source", "term_name", "p_value"))) # 11 rows 
# 
# write.csv(common_pathways_highsig, "common_pathways_high.csv")

# Specific pathways in sepsis but not in no sepsis #### 

setdiff(pathways_pc8_neg$result$term_id, pathways_pc15_neg$result$term_id) # 125 pathways 

specific_pathways_pc8 <- anti_join(pathways_pc8_neg$result, pathways_pc15_neg$result, by = "term_id") %>% 
  dplyr::select(-contains(c("parents", "evidence_codes", "query", "significant", 
                            "effective_domain_size", "source_order")))

nrow(specific_pathways_pc8) # 125 

write.csv(specific_pathways_pc8, "specific_pathways_pc8.csv")

# writexl::write_xlsx(list(sheet1 = specific_pathways), 
#                     "specific_pathways_sepsis_only.xlsx")

specific_pathways_pc15 <- anti_join(pathways_pc15_neg$result, pathways_pc8_neg$result, by = "term_id") %>% 
  dplyr::select(-contains(c("parents", "evidence_codes", "query", "significant", 
                            "effective_domain_size", "source_order")))

nrow(specific_pathways_pc15) # 37 

write.csv(specific_pathways_pc15, "specific_pathways_pc15.csv")


# Plotting #### 

gostplot(pathways_vaso, capped = TRUE, interactive = TRUE)

gostplot(pathways_vaso, capped = TRUE, interactive = TRUE)

gostplot(specific_pathways, capped = TRUE, interactive = TRUE) # cant get that, mus be a whole gprofiler output 
# Also neutrophil degranulation as a top pathway 

# Vasopressors all groups  
vaso_plot <- gostplot(pathways_vaso, capped = TRUE, interactive = FALSE)
vaso_plot2 <- publish_gostplot(vaso_plot, 
                          highlight_terms = c("GO:1901564",
                                              "GO:0016192",
                                              "GO:0008219",
                                              "GO:0030097",
                                              "REAC:R-HSA-6798695", 
                                              "REAC:R-HSA-168249", 
                                              "REAC:R-HSA-199991",
                                              "REAC:R-HSA-204005",
                                              "REAC:R-HSA-189451",
                                              "REAC:R-HSA-5687128"), 
                          width = NA, height = NA, filename = "pathways_vasopressors_allgroups.png")
# cluster 0 and 2 cartilage dev only - 1 pathway 
ggplot2::ggsave("pathways_vasopressors_allgroups.png", plot = vaso_plot2, scale = 1.5)

# Vasopressors sepsis only  
vaso_plot_sep <- gostplot(pathways_vaso_sep, capped = TRUE, interactive = FALSE)
vaso_plot2_sep <- publish_gostplot(vaso_plot_sep, 
                               highlight_terms = c("GO:1901564",
                                                   "GO:0006914",
                                                   "GO:0006950",
                                                   "REAC:R-HSA-6798695", 
                                                   "REAC:R-HSA-1234174", 
                                                   "REAC:R-HSA-1234176",
                                                   "REAC:R-HSA-8941858",
                                                   "REAC:R-HSA-180534",
                                                   "REAC:R-HSA-5687128"), 
                               width = NA, height = NA, filename = "pathways_vasopressors_sepsis_only.png")
# cluster 0 and 2 cartilage dev only - 1 pathway 
ggplot2::ggsave("pathways_vasopressors_sepsis_only.png", plot = vaso_plot2, scale = 1.25)

# ** Table #### 
publish_gosttable(pathways_vaso, highlight_terms = c("GO:1901564",
                                                    "GO:0016192",
                                                             "GO:0008219",
                                                             "GO:0030097",
                                                             "REAC:R-HSA-6798695", 
                                                             "REAC:R-HSA-168249", 
                                                             "REAC:R-HSA-199991",
                                                             "REAC:R-HSA-204005",
                                                             "REAC:R-HSA-189451",
                                                             "REAC:R-HSA-5687128"), 
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

ggsave("vaso_gosttable.jpeg", scale = 1.2)
# Saving 10 x 4.55 in image


publish_gosttable(pathways_vaso_sep, highlight_terms = c("GO:1901564",
                                                     "GO:0006914",
                                                     "GO:0006950",
                                                     "REAC:R-HSA-6798695", 
                                                     "REAC:R-HSA-1234174", 
                                                     "REAC:R-HSA-1234176",
                                                     "REAC:R-HSA-8941858",
                                                     "REAC:R-HSA-180534",
                                                     "REAC:R-HSA-5687128"), 
                  use_colors = TRUE, 
                  show_columns = c("source", "term_name", "term_size", "intersection_size"),
                  filename = NULL)

ggsave("vaso_sep_gosttable.jpeg", width = 12, height = 5, units = "in")


# Barchart #### 

# need to pick pathways I am interested in first 



# Genes within pathways - potentially separate file 
