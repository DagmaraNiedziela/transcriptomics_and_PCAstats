# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
source("functions_PCA-transcriptomics.R")

# Revisit data #### 

# prcomp was done on vsd top 500 most variable genes, with default options - file PCA_prcomp.R 

# rv <- rowVars(assay(vsd))
# select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
#                                                    length(rv)))]
# 
# my_pca <- prcomp(t(assay(vsd)[select,]), scale = FALSE,
#                  center = TRUE, retx = TRUE)

# Prepare data to plot - revisit 

# head(my_pca$x)
# 
# pc_scores <- my_pca$x 
# 
# pc_scores <- pc_scores %>% 
#   as_tibble(rownames = "id2")
# 
pc_scores
# max(pc_scores$PC1)

# PCAdata with all PCs and all metadata #### 
pc_scores_toplot <- pc_scores %>% # id2 is at the end
  dplyr::full_join(Kalantar_metadata, by = "id2") 

# Add apacheiii and wbc LOGS 
pc_scores_toplot <- pc_scores_toplot %>% mutate(apacheiii_log = log10(apacheiii + 1),
                              wbc_max_log = log10(wbc_max + 1))

# Check if it works 
pc_scores_toplot %>%
  ggplot(aes(x = PC1, y = PC3, color = group)) +
  geom_point()

# get percent variance of each PC 

percent_Var <- round(propve * 100, 2)
# propve calculated in PCA_prcomp.R script 

# Plot PCs 3 4 and 5 #### 

# variables to plot
colnames(pc_scores_toplot)

metadata_for_newpcs <- c("group","on_pressors", "x28d_death", "sirs_total")
"apacheiii_log" 

# pc1 vs pc3 
plot_list_pc1v3 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, 
                          data = pc_scores_toplot, x = "PC1", y = "PC3", var1 = 1, var2 = 3)
plot_list_pc1v3[[1]]

# Plus single plot for apacheiii, include 
plot_list_pc1v3[[5]] <- plot_colored_PCA_freePC_cont(data = pc_scores_toplot,
                                         color = "apacheiii_log",
                                         x = "PC1", y = "PC3", var1 = 1, var2 = 3)

ggpubr::ggarrange(plotlist = plot_list_pc1v3)
ggsave("PCA/PCA_PC1vs3.jpeg", width = 18, height = 11.2, units = "in")

# pc3 vs 4 
plot_list_pc3v4 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, 
                          data = pc_scores_toplot, x = "PC3", y = "PC4", var1 = 3, var2 = 4)
plot_list_pc3v4[[1]]

# Plus single plot for apacheiii, include 
plot_list_pc3v4[[5]] <- plot_colored_PCA_freePC_cont(data = pc_scores_toplot,
                                                     color = "apacheiii_log",
                                                     x = "PC3", y = "PC4", var1 = 3, var2 = 4)

ggpubr::ggarrange(plotlist = plot_list_pc3v4)
ggsave("PCA/PCA_PC3vs4.jpeg", width = 18, height = 11.2, units = "in")

# pc4 vs 5 
plot_list_pc4v5 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, 
                          data = pc_scores_toplot, x = "PC4", y = "PC5", var1 = 4, var2 = 5)
plot_list_pc4v5[[1]]

# Plus single plot for apacheiii, include 
plot_list_pc4v5[[5]] <- plot_colored_PCA_freePC_cont(data = pc_scores_toplot,
                                                     color = "apacheiii_log",
                                                     x = "PC4", y = "PC5", var1 = 4, var2 = 5)

ggpubr::ggarrange(plotlist = plot_list_pc4v5)
ggsave("PCA/PCA_PC4vs5.jpeg", width = 18, height = 11.2, units = "in")

# pc2 vs 3 
plot_list_pc2v3 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, 
                          data = pc_scores_toplot, x = "PC2", y = "PC3", var1 = 2, var2 = 3)
plot_list_pc2v3[[2]]

# Plus single plot for apacheiii, include 
plot_list_pc2v3[[5]] <- plot_colored_PCA_freePC_cont(data = pc_scores_toplot,
                                                     color = "apacheiii_log",
                                                     x = "PC2", y = "PC3", var1 = 2, var2 = 3)

ggpubr::ggarrange(plotlist = plot_list_pc2v3)
ggsave("PCA/PCA_PC2vs3.jpeg", width = 18, height = 11.2, units = "in")


# Vasopressors in all PC combinations 
ggpubr::ggarrange(plot_list_pc1v3[[2]], plot_list_pc2v3[[2]], plot_list_pc3v4[[2]], 
                  plot_list_pc4v5[[2]])
ggsave("PCA/pressors_PCs2-5.jpeg", width = 12, height = 10, units = "in")


# THIS AFTER STATS #### 
# Gene loadings for later PCs #### 

# pc_loadings <- my_pca$rotation
# 
# pc_loadings <- pc_loadings %>% 
#   as_tibble(rownames = "gene")

pc_loadings


# Biomart Explain top genes ####

# Function used for this - straight from loadings data frame 
gene_names_pc5 <- get_loadings_gene_list(pca_data = pc_loadings,
                                             pc = "PC5",
                                             direction = "positive",
                                             mart = mart, number = 20)

gene_names_pc5

write.csv(gene_names_pc5, "PCA/top_loadings_genes_pc5.csv")

# Function used for this - straight from loadings data frame 
gene_names_pc5_neg <- get_loadings_gene_list(pca_data = pc_loadings,
                                              pc = "PC5",
                                              direction = "negative",
                                              mart = mart, number = 20)

gene_names_pc5_neg

write.csv(gene_names_pc5_neg, "PCA/top_loadings_genes_pc5_neg.csv")

# Excel with both 
writexl::write_xlsx(list(positive_pc5 = gene_names_pc5,
                         negative_pc5 = gene_names_pc5_neg),
                    "PCA/top_loadings_PC5_plasma.xlsx")



# Make data that includes gene loadings for selected PCs, log transformed #### 

# pc5 - this can be done using function above 

PCAdata_loadings_pc5 <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames, 
                                              pca_data = pc_scores_toplot,
                                              loadings_pc_names = gene_names_pc5)

PCAdata_loadings_pc5

# the hgnc symbol has an empty string, lapply throws error 
PCAdata_loadings_pc5_neg <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames, 
                                                  pca_data = pc_scores_toplot,
                                                  loadings_pc_names = gene_names_pc5_neg)

PCAdata_loadings_pc5_neg

# Plot top loadings of later PCs #### 

# pc1 vs pc5 positive 
plot_list_loadings_pc1v5 <- lapply(gene_names_pc5$hgnc_symbol, 
                                   FUN = plot_colored_PCA_freePC_cont, 
                                   data = PCAdata_loadings_pc5, 
                                   x = "PC1", y = "PC5", var1 = 1, var2 = 8)
plot_list_loadings_pc1v5[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings_pc1v5)
ggsave("PCA/PCA_loadings_PC1vs5.jpeg", width = 18, height = 11.2, units = "in")

# pc1 vs pc5 neg 
plot_list_loadings_pc1v5_neg <- lapply(stringi::stri_remove_empty(gene_names_pc5_neg$hgnc_symbol), 
                                       FUN = plot_colored_PCA_freePC_cont, 
                                      data = PCAdata_loadings_pc5_neg, 
                                      x = "PC1", y = "PC5", var1 = 1, var2 = 5)
plot_list_loadings_pc1v5_neg[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings_pc1v5_neg)
ggsave("PCA/PCA_loadings_PC1vs5_neg.jpeg", width = 18, height = 11.2, units = "in")


# Plot vasopressors on PC5 #### 

plot_colored_PCA_freePC(data = PCAdata_loadings_pc5, 
                        color = "on_pressors", x = "PC1", y = "PC5", var1 = 1, var2 = 5)
ggsave("PCA/PCA_pressors_PC1vs5.jpeg")

# plot_colored_PCA_freePC(data = PCAdata_loadings_pc5_neg, 
#                         color = "on_pressors", x = "PC1", y = "PC5", var1 = 1, var2 = 5)
# ggsave("PCA/PCA_pressors_PC1vs5_neg.jpeg")
