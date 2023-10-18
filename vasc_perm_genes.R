# Library calls #### 
library(dplyr)
library(biomaRt)

# Vascular permeability genes #### 

genes_vasc_perm <- c("PECAM1", "ICAM1", "VCAM1", "DPP3", "ANGPT2", "IL6", "ALB", 
                     "SELE", "CRP", "SDC1", "ADM", "VEGF", "VWF", "CD142", "F2", 
                     "SERPINE1", "PROC", "THBD", "ENG") 
length(genes_vasc_perm)
# SELE also LYAM2 
# CD142 also F3 - tissue factor 

# Biomart for vascular permeability genes 

gene_names_vasc_perm <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","description"),
                              values = genes_vasc_perm, mart= mart)
gene_names_vasc_perm # missing CD142 and VEGF 
write.csv(gene_names_vasc_perm, "genes_vascular_permeability.csv")


# Vascular permeability PCA data (log10 transformed) #### 


PCAdata_vasc_perm <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames, 
                                                  pca_data = pc_scores_toplot,
                                                  loadings_pc_names = gene_names_vasc_perm)

PCAdata_vasc_perm

# Look for vascular permeability genes in differential expression analysis #### 

# ** Vasopressors yes no - all groups #### 

results_vaso_shrink_DEG

vasc_perm_vasoDEG <- inner_join(gene_names_vasc_perm, as.data.frame(results_vaso_shrink_DEG), by = c("hgnc_symbol" = "gene_symbol"))
vasc_perm_vasoDEG

# ** Vasopressors yes no - sepsis only #### 

results_vaso_sep_shrink_DEG

vasc_perm_vasoDEG_sep <- inner_join(gene_names_vasc_perm, as.data.frame(results_vaso_sep_shrink_DEG), by = c("hgnc_symbol" = "gene_symbol"))
vasc_perm_vasoDEG_sep

# ** Condition - group x vasopressors #### 

subset(results_vaso_sepsis_nonbsi_cond_shrink, padj < 0.05)

vasc_perm_vaso_group_nonbsi <- inner_join(gene_names_vasc_perm, 
                                          as.data.frame(subset(results_vaso_sepsis_nonbsi_cond_shrink, padj < 0.05)), 
                                          by = c("hgnc_symbol" = "gene_symbol"))
vasc_perm_vaso_group_nonbsi

subset(results_vaso_no_sepsis_cond_shrink, padj < 0.05)

vasc_perm_vaso_group_nosepsis <- inner_join(gene_names_vasc_perm, 
                                            as.data.frame(subset(results_vaso_no_sepsis_cond_shrink, padj < 0.05)), 
                                            by = c("hgnc_symbol" = "gene_symbol"))

vasc_perm_vaso_group_nosepsis
