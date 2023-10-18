# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(biomaRt)

# Revisit data #### 

# prcomp was done on vsd top 500 most variable genes, with default options 
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]

my_pca <- prcomp(t(assay(vsd)[select,]), scale = FALSE,
                 center = TRUE, retx = TRUE)


# Prepare data to plot #### 

head(my_pca$x)

pc_scores <- my_pca$x 

pc_scores <- pc_scores %>% 
  as_tibble(rownames = "id2")

pc_scores$id2 <- pc_scores$ID2

pc_scores
max(pc_scores$PC1)

colnames(pc_scores)
colnames(Kalantar_metadata)

pc_scores$id2
Kalantar_metadata$id2

pc_scores_toplot <- pc_scores[,c(1:12,223)] %>% # id2 is at the end
  dplyr::full_join(Kalantar_metadata, by = "id2") 

# Clean names 
pc_scores_toplot <- pc_scores_toplot %>% janitor::clean_names()

# Clean up metadata columns (already done in metadata on second pass)

pc_scores_toplot$x28d_death <- ifelse(pc_scores_toplot$x28d_death == 0, "no", "yes")
pc_scores_toplot$intubated <- ifelse(pc_scores_toplot$intubated == 0, "no", "yes")
pc_scores_toplot$on_pressors <- ifelse(pc_scores_toplot$on_pressors == 0, "no", "yes")
pc_scores_toplot$immunocompromised <- ifelse(pc_scores_toplot$immunocompromised == 0, "no", "yes")
pc_scores_toplot$sirs_total <- factor(pc_scores_toplot$sirs_total)
pc_scores_toplot$age <- as.numeric(as.character(pc_scores_toplot$age))

# pc_scores_toplot %>% 
#   ggplot(aes(x = PC1, y = PC3, color = Group)) +
#   geom_point()

# get percent variance of each PC 

percent_Var <- round(propve * 100, 2)
# propve calculated in transcriptomics_PCA script 


# Take PCs 3 4 and 5 and plot #### 

# Plot scatter plots with PC3 4 and 5 against PC1 see what comes up, 3 vs 4, 4 vs 5 
# check vasopressors, death and sirs on these plots , add intubation 
# and group

# variables to plot
colnames(pc_scores_toplot)

metadata_for_newpcs <- c("group","on_pressors", "x28d_death", "sirs_total", "intubated")

# Formatted plot function 
# var1 number of a PC on x axis 
# var2 number of a PC on y axis 
plot_colored_PCA_freePC <- function(data, color, x, y, var1, var2){
  ggplot(data = data, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) +
    geom_point(size =3) + 
    xlab(paste0("PC", var1, ": ", percent_Var[var1], "% variance")) +
    ylab(paste0("PC", var2, ": ", percent_Var[var2], "% variance")) + 
    #coord_fixed() 
    scale_color_brewer(palette="Set1", direction = -1) + 
    scale_fill_brewer(palette = "Set1", direction = -1) + 
    theme(text = element_text(size = 14, family = "Calibri")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}

plot_colored_PCA_freePC_cont <- function(data, color, x, y, var1, var2){
  ggplot(data = data, aes(x = .data[[x]], y = .data[[y]], color = .data[[color]])) +
    geom_point(size =3) + 
    xlab(paste0("PC", var1, ": ", percent_Var[var1], "% variance")) +
    ylab(paste0("PC", var2, ": ", percent_Var[var2], "% variance")) + 
    #coord_fixed() 
    scale_color_viridis_b() + 
  #  scale_fill_brewer(palette = "Set1", direction = -1) + 
    theme(text = element_text(size = 14, family = "Calibri")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}


# lapply plots #### 

# pc1 vs pc3 
plot_list_pc1v3 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, data = pc_scores_toplot, x = "pc1", y = "pc2", var1 = 1, var2 = 3)
plot_list_pc1v3[[1]]

ggpubr::ggarrange(plotlist = plot_list_pc1v3)
ggsave("PCA_PC1vs3.jpeg", width = 18, height = 11.2, units = "in")

# pc3 vs 4 
plot_list_pc3v4 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, data = pc_scores_toplot, x = "pc3", y = "pc4", var1 = 3, var2 = 4)
plot_list_pc3v4[[1]]

ggpubr::ggarrange(plotlist = plot_list_pc3v4)
ggsave("PCA_PC3vs4.jpeg", width = 18, height = 11.2, units = "in")
# Saving 18 x 11.2 in image

# pc4 vs 5 
plot_list_pc4v5 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, data = pc_scores_toplot, x = "pc4", y = "pc5", var1 = 4, var2 = 5)
plot_list_pc4v5[[1]]

ggpubr::ggarrange(plotlist = plot_list_pc4v5)
ggsave("PCA_PC4vs5.jpeg", width = 18, height = 11.2, units = "in")

# pc2 vs 3 
plot_list_pc2v3 <- lapply(metadata_for_newpcs, FUN = plot_colored_PCA_freePC, data = pc_scores_toplot, x = "pc2", y = "pc3", var1 = 2, var2 = 3)
plot_list_pc2v3[[2]]

ggpubr::ggarrange(plotlist = plot_list_pc2v3)
ggsave("PCA_PC2vs3.jpeg", width = 18, height = 11.2, units = "in")

ggpubr::ggarrange(plot_list_pc1v3[[2]], plot_list_pc2v3[[2]], plot_list_pc3v4[[2]], plot_list_pc4v5[[2]])
ggsave("pressors_PCs2-5.jpeg", width = 12, height = 10, units = "in")


# Gene loadings for later PCs #### 

# ** Loadings #### 

pc_loadings <- my_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

pc_loadings

# top genes with highest loadings on PC1 and PC2 
top_genes_laterpcs <- pc_loadings %>% 
  select(gene, PC1, PC2, PC8, PC15) # %>% arrange(desc(PC8))
#  tidyr::pivot_longer(matches("PC"), names_to = "PC", values_to = "loading")  

pc_loadings %>% 
  select(gene, PC1, PC2, PC8, PC15) %>% arrange(PC15)

top_genes_laterpcs

top_genes_laterpcs

# positive 
top_loadings_pc8 <- top_genes_laterpcs$gene[with(top_genes_laterpcs, order(-PC8))][1:20]

# negative
# top_loadings_pc8_neg <- top_genes_laterpcs$gene[with(top_genes_laterpcs, order(PC8))][1:20]
# 
# top_loadings_pc8_neg

# positive 
top_loadings_pc15 <- top_genes_laterpcs$gene[with(top_genes_laterpcs, order(-PC15))][1:20]
top_genes_laterpcs$gene[with(top_genes_laterpcs, order(-top_genes_laterpcs[["PC15"]]))][1:20]

# negative 
# top_loadings_pc15_neg <- top_genes_laterpcs$gene[with(top_genes_laterpcs, order(top_genes_laterpcs[["PC15"]]))][1:20]
# 
# top_loadings_pc15_neg

#' Function to get top 20 loadings for a given pc and save a csv with them 
#' @param pca_data a data frame from pc loadings containing selected pc loadings and gene ensembl id in a "gene" column,
#' PC column names with capital letters "PC1" etc. 
#' @param pc a string - PC and number 
#' @param a string, values "positive" or "negative" for positive and negative loadings respectively 
#' @param mart a dataset from biomart, preloaded (mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl")))
#' @param number an integer, number of loading genes to get, default is 20, but for other analyses I might need more
#' mart should be pre-loaded for this function
#' 
#' @return a csv file with + a data frame result 
get_loadings_gene_list <- function(pca_data, pc, direction, mart = mart, number = 20){
  if (direction == "positive"){
    top_loadings <- pca_data$gene[with(pca_data, order(-pca_data[[pc]]))][1:number]
  }
  else {
    top_loadings <- pca_data$gene[with(pca_data, order(pca_data[[pc]]))][1:number]
  }
  # fix gene names 
  top_loadings <- gsub("\\.\\d{2,2}$", "", top_loadings)
  top_loadings <- gsub("\\.\\d{1,1}$", "", top_loadings)
  # get gene names, save and return 
  gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                          values = top_loadings, mart= mart)
  # with vsd object the top loadings changed!
  write.csv(gene_names, paste0("gene_names_",pc,"_",direction,".csv"))
  return(gene_names)
}


# ** Biomart Explain top genes ####

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

top_loadings_pc8_2 <- gsub("\\.\\d{2,2}$", "", top_loadings_pc8)
top_loadings_pc8_2 <- gsub("\\.\\d{1,1}$", "", top_loadings_pc8_2)
top_loadings_pc8_2

top_loadings_pc15_2 <- gsub("\\.\\d{2,2}$", "", top_loadings_pc15)
top_loadings_pc15_2 <- gsub("\\.\\d{1,1}$", "", top_loadings_pc15_2)
top_loadings_pc15_2

# top_loadings_pc15_neg <- gsub("\\.\\d{2,2}$", "", top_loadings_pc15_neg)
# top_loadings_pc15_neg <- gsub("\\.\\d{1,1}$", "", top_loadings_pc15_neg)
# top_loadings_pc15_neg

# ** gene names #### 
gene_names_pc8 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                    values = top_loadings_pc8_2, mart= mart)
gene_names_pc8
# with vsd object the top loadings changed! 
write.csv(gene_names_pc8, "top_loadings_genes_pc8.csv")

# Function used for this - straight from loadings data frame 
gene_names_pc8_neg <- get_loadings_gene_list(pca_data = top_genes_laterpcs,
                                             pc = "PC8",
                                             direction = "negative",
                                             mart = mart)

gene_names_pc8_neg

gene_names_pc15 <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                        values = top_loadings_pc15_2, mart= mart)
gene_names_pc15
# with vsd object the top loadings changed! 
write.csv(gene_names_pc15, "top_loadings_genes_pc15.csv")

# Function used for this - straight from loadings data frame 
gene_names_pc15_neg <- get_loadings_gene_list(pca_data = top_genes_laterpcs,
                                              pc = "PC15",
                                              direction = "negative",
                                              mart = mart)

gene_names_pc15_neg

# gene_names_pc15_neg <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
#                          values = top_loadings_pc15_neg, mart= mart)
# gene_names_pc15_neg
# # with vsd object the top loadings changed! 
# write.csv(gene_names_pc15, "top_loadings_genes_pc15_neg.csv")


# make data that includes gene loadings for pc8 and pc15, log transformed #### 

#' Function to make a data frame with genes with top loadings and medatata
#' 
#' @param counts_data a counts data frame with samples as columns and genes as additional column geneId 
#' (ensembl id) and hgnc_symbol (gene symbol)
#' @param pca_data a data frame like pc_scores_toplot - contains pcs, metadata, and sample name column named as "id2"
#' @param loadings_pc_names the gene_names for pc loadings - returned data frame of the get_loadings_gene_list function
#' 
#' @return a data frame with genes with top loadings and medatata - to plot on a PCA plot colored by gene count 
make_loadings_pcadata <- function(counts_data, pca_data, loadings_pc_names){
  # pick gene counts for all samples with only our gene loadings 
  counts_pc <- counts_data %>% 
    dplyr::filter(hgnc_symbol %in% loadings_pc_names[["hgnc_symbol"]])
  
  # Transpose and clean - top loadings 
  counts_pc_to_join <- t(counts_pc)
  colnames(counts_pc_to_join) <- counts_pc_to_join["hgnc_symbol",]
  counts_pc_to_join <- as.data.frame(counts_pc_to_join)
  counts_pc_to_join$name <- rownames(counts_pc_to_join)
  
  counts_pc_to_join <- slice(counts_pc_to_join, 1:(n() - 2))  

  pca_data$name <- pca_data$id2
  pcadata_loadings <- full_join(pca_data, counts_pc_to_join, by = "name")
  
  # Columns to numeric
  pcadata_loadings[, stringi::stri_remove_empty(loadings_pc_names[["hgnc_symbol"]])] <- lapply(stringi::stri_remove_empty(loadings_pc_names[["hgnc_symbol"]]), function(x) as.numeric(as.character(pcadata_loadings[[x]])))

  # Log transform counts
  pcadata_loadings[, stringi::stri_remove_empty(loadings_pc_names[["hgnc_symbol"]])] <- lapply(stringi::stri_remove_empty(loadings_pc_names[["hgnc_symbol"]]), function(x) log10(pcadata_loadings[[x]] + 1))
  
  return(pcadata_loadings)
}

# stringi::stri_remove_empty(gene_names_pc8_neg$hgnc_symbol)


# pc8 - this can be done using fucntion above 
Kalantar_counts_loadings_pc8 <- Kalantar_counts_fornames %>% 
  filter(hgnc_symbol %in% gene_names_pc8$hgnc_symbol)

Kalantar_counts_loadings_pc8

# Transpose and clean - top loadings 
Kalantar_counts_loadings_pc8_tojoin <- t(Kalantar_counts_loadings_pc8)
colnames(Kalantar_counts_loadings_pc8_tojoin) <- Kalantar_counts_loadings_pc8_tojoin["hgnc_symbol",]
Kalantar_counts_loadings_pc8_tojoin <- as.data.frame(Kalantar_counts_loadings_pc8_tojoin)
Kalantar_counts_loadings_pc8_tojoin$name <- rownames(Kalantar_counts_loadings_pc8_tojoin)

Kalantar_counts_loadings_pc8_tojoin <- slice(Kalantar_counts_loadings_pc8_tojoin, 1:(n() - 2))  

head(Kalantar_counts_loadings_pc8_tojoin)
pc_scores_toplot

pc_scores_toplot$name <- pc_scores_toplot$id2
PCAdata_loadings_pc8 <- full_join(pc_scores_toplot, Kalantar_counts_loadings_pc8_tojoin, by = "name")

PCAdata_loadings_pc8
colnames(PCAdata_loadings_pc8)[38]
lapply(PCAdata_loadings_pc8, class)

# Columns to numeric 
PCAdata_loadings_pc8[, 38:ncol(PCAdata_loadings_pc8)] <- lapply(38:ncol(PCAdata_loadings_pc8), function(x) as.numeric(as.character(PCAdata_loadings_pc8[[x]])))

PCAdata_loadings_pc8[, gene_names_pc8$hgnc_symbol] <- lapply(gene_names_pc8$hgnc_symbol, function(x) as.numeric(as.character(PCAdata_loadings_pc8[[x]])))
# Log transform counts 
PCAdata_loadings_pc8[, gene_names_pc8$hgnc_symbol] <- lapply(gene_names_pc8$hgnc_symbol, function(x) log10(PCAdata_loadings_pc8[[x]] + 1))

# pc15 - this all can be done using a function above now 
Kalantar_counts_loadings_pc15 <- Kalantar_counts_fornames %>% 
  filter(hgnc_symbol %in% gene_names_pc15$hgnc_symbol)

# Transpose and clean - top loadings 
Kalantar_counts_loadings_pc15_tojoin <- t(Kalantar_counts_loadings_pc15)
colnames(Kalantar_counts_loadings_pc15_tojoin) <- Kalantar_counts_loadings_pc15_tojoin["hgnc_symbol",]
Kalantar_counts_loadings_pc15_tojoin <- as.data.frame(Kalantar_counts_loadings_pc15_tojoin)
Kalantar_counts_loadings_pc15_tojoin$name <- rownames(Kalantar_counts_loadings_pc15_tojoin)

Kalantar_counts_loadings_pc15_tojoin <- slice(Kalantar_counts_loadings_pc15_tojoin, 1:(n() - 2))  

pc_scores_toplot2
pc_scores_toplot2$name <- pc_scores_toplot2$id2
PCAdata_loadings_pc15 <- full_join(pc_scores_toplot2, Kalantar_counts_loadings_pc15_tojoin, by = "name")

PCAdata_loadings_pc15
colnames(PCAdata_loadings_pc15)

# Columns to numeric 
PCAdata_loadings_pc15[, 81:ncol(PCAdata_loadings_pc15)] <- lapply(81:ncol(PCAdata_loadings_pc15), function(x) as.numeric(as.character(PCAdata_loadings_pc15[[x]])))

# Log transform counts 
PCAdata_loadings_pc15[, 81:ncol(PCAdata_loadings_pc15)] <- lapply(81:ncol(PCAdata_loadings_pc15), function(x) log10(PCAdata_loadings_pc15[[x]] + 1))

# the hgnc symbol has an empty string, lapply throws error 
PCAdata_loadings_pc8_neg <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames, 
                                                  pca_data = pc_scores_toplot,
                                                  loadings_pc_names = gene_names_pc8_neg)

PCAdata_loadings_pc8_neg

# the hgnc symbol has an empty string, lapply throws error 
PCAdata_loadings_pc15_neg <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames, 
                                                  pca_data = pc_scores_toplot2,
                                                  loadings_pc_names = gene_names_pc15_neg)

PCAdata_loadings_pc15_neg

# top loadings on a PCA #### 

# pc1 and PC8 
# plot_list_loadings_pc8_log_nof <- lapply(gene_names_pc8$hgnc_symbol, FUN = plot_colored_PCA_cont_nofacet, data = PCAdata_loadings_pc8) 
# plot_list_loadings_pc8_log_nof[[1]]
# 
# ggpubr::ggarrange(plotlist = plot_list_loadings_pc8_log_nof)
# ggsave("PCA_multi_color_topload_log_nof.jpeg", scale = 1.2)


# pc1 vs pc8 
plot_list_loadings_pc1v8 <- lapply(gene_names_pc8$hgnc_symbol, FUN = plot_colored_PCA_freePC_cont, data = PCAdata_loadings_pc8, x = "pc1", y = "pc8", var1 = 1, var2 = 8)
plot_list_loadings_pc1v8[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings_pc1v8)
ggsave("PCA_loadings_PC1vs8.jpeg", width = 18, height = 11.2, units = "in")

# pc1 and PC15 
plot_list_loadings_pc1v15 <- lapply(gene_names_pc15$hgnc_symbol, FUN = plot_colored_PCA_freePC_cont, data = PCAdata_loadings_pc15, x = "pc1", y = "pc15", var1 = 1, var2 = 15)
plot_list_loadings_pc1v15[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings_pc1v15)
ggsave("PCA_loadings_PC1vs15.jpeg", width = 18, height = 11.2, units = "in")

# pc1 vs pc8 neg 
plot_list_loadings_pc1v8_neg <- lapply(stringi::stri_remove_empty(gene_names_pc8_neg$hgnc_symbol), FUN = plot_colored_PCA_freePC_cont, 
                                   data = PCAdata_loadings_pc8_neg, x = "pc1", y = "pc8", var1 = 1, var2 = 8)
plot_list_loadings_pc1v8_neg[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings_pc1v8_neg)
ggsave("PCA_loadings_PC1vs8_neg.jpeg", width = 18, height = 11.2, units = "in")

# pc1 vs pc15 neg 
plot_list_loadings_pc1v15_neg <- lapply(stringi::stri_remove_empty(gene_names_pc15_neg$hgnc_symbol), FUN = plot_colored_PCA_freePC_cont, 
                                       data = PCAdata_loadings_pc15_neg, x = "pc1", y = "pc8", var1 = 1, var2 = 15)
plot_list_loadings_pc1v15_neg[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings_pc1v15_neg)
ggsave("PCA_loadings_PC1vs15_neg.jpeg", width = 18, height = 11.2, units = "in")


# Plot vasopressors on PC8 and PC15 #### 

plot_colored_PCA_freePC(data = PCAdata_loadings_pc8, color = "on_pressors", x = "pc1", y = "pc8", var1 = 1, var2 = 8)
ggsave("PCA_pressors_PC1vs8.jpeg")

plot_colored_PCA_freePC(data = PCAdata_loadings_pc15, color = "on_pressors", x = "pc1", y = "pc15", var1 = 1, var2 = 15)
ggsave("PCA_pressors_PC1vs15.jpeg")
