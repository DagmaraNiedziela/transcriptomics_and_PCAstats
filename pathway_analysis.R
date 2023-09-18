# Library calls #### 
# library(gProfileR)
library(dplyr)
library(gprofiler2)
library(ggplot2)
library(pheatmap)
library(Cairo)

# Data #### 

results_vaso_shrink_DEG 

results_vaso_sep_shrink_DEG

# results_vaso_shrink_DEG <- read.csv("DEG_sig_vasopressors.csv")
# results_vaso_sep_shrink_DEG <- read.csv("DEG_sig_vasopressors_sepsis_only.csv")


# Pathway analysis ####

# Vasopressors, all groups 
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

# Vasopressors, sepsis only 

DEG_gene_names_sep <- rownames(results_vaso_sep_shrink_DEG)
DEG_gene_names_sep

pathways_vaso_sep <- gost(DEG_gene_names_sep,
                      organism = "hsapiens", ordered_query = TRUE, 
                      multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                      measure_underrepresentation = FALSE, evcodes = TRUE, 
                      user_threshold = 0.05, correction_method = "g_SCS", 
                      domain_scope = "annotated", custom_bg = NULL, 
                      numeric_ns = "", sources = c("GO:BP", "GO:MF", "GO:CC","KEGG","REAC"), as_short_link = FALSE)

writexl::write_xlsx(list(sheet1 = pathways_vaso_sep$result[,!(names(pathways_vaso_sep$result) %in% c("parents", "evidence_codes", "intersection"))]), 
                    "pathways_vasopressors_sepsis_only.xlsx")

head(pathways_vaso_sep$result)
nrow(pathways_vaso_sep$result %>% filter(precision > 0.35)) # 30 
pathways_vaso_sep$result %>% filter(precision > 0.35) %>% dplyr::select(source, term_name, p_value)
# 14 GO:BP 


# Common pathways in 2 sets #### 

intersect(pathways_vaso$result$term_id, pathways_vaso_sep$result$term_id) # 320 pathways in common 

common_pathways <- inner_join(pathways_vaso$result, pathways_vaso_sep$result, by = "term_id", suffix = c("_vaso", "_vaso_sep")) %>% 
  dplyr::select(-contains(c("parents", "evidence_codes", "intersection_v", "query", "significant", 
                            "effective_domain_size", "source_order")))

nrow(common_pathways) # 319 
write.csv(common_pathways, "common_pathways.csv")

writexl::write_xlsx(list(sheet1 = common_pathways), 
                    "common_pathways.xlsx")

typeof(common_pathways$precision_vaso)
common_pathways_highsig <- common_pathways %>% filter(precision_vaso > 0.35 & precision_vaso_sep > 0.35) %>% 
  dplyr::select(contains(c("precision", "source", "term_name", "p_value"))) # 27 rows 

write.csv(common_pathways_highsig, "common_pathways_high.csv")

# Plotting #### 

gostplot(pathways_vaso, capped = TRUE, interactive = TRUE)

gostplot(pathways_vaso_sep, capped = TRUE, interactive = TRUE)
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
