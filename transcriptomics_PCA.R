# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readxl)
library("biomaRt")
library(ggbiplot) # loads plyr!!!! need to specify dplyr functions possibly 

# Load data #### 

Kalantar_counts <- read.csv("GSE189400_earli_paxgene_counts_50k.csv", header = TRUE)
head(Kalantar_counts)

Kalantar_metadata <- read_excel("Kalantar_supplementary.xlsx", sheet = "Data_16")
head(Kalantar_metadata)

ncol(Kalantar_counts) # 379
nrow(Kalantar_metadata) # 221 
ncol(Kalantar_counts) == nrow(Kalantar_metadata) # FALSE

Kalantar_metadata$ID
unique(Kalantar_metadata$ID) # same, all unique

colnames(Kalantar_counts)
# duplicate columns for some samples - why??? 

Kalantar_counts_reduced <- Kalantar_counts %>% select(-ends_with(".1"))
colnames(Kalantar_counts_reduced)
length(unique(colnames(Kalantar_counts_reduced))) # 348 

# Need to split EARLI from colnames, or add EARLI to ID 
Kalantar_metadata$ID2 <- paste("EARLI", Kalantar_metadata$ID, sep = "_")

intersect(Kalantar_metadata$ID2, colnames(Kalantar_counts_reduced))

# Final counts only contained in metadata 
Kalantar_counts_intersected <- Kalantar_counts_reduced %>% select(X, Kalantar_metadata$ID2) 
# This also makes sure the counts are ordered the same way as metadata !!! 

# DESeq object #### 

rownames(Kalantar_counts_intersected) <- Kalantar_counts_intersected$X
Kalantar_counts_intersected <- Kalantar_counts_intersected[,-1]

dds <- DESeqDataSetFromMatrix(countData = Kalantar_counts_intersected,
                              colData = Kalantar_metadata, 
                              design = ~ Group) 

nrow(dds) #27097 

# DESeq PCA #### 

DESeq2::vst

# Normalise data with a vst function 
vsd <- vst(dds, blind=FALSE)

# assay(vsd)
DESeq2::plotPCA
BiocGenerics::plotPCA
showMethods(plotPCA)
selectMethod(plotPCA, signature = "DESeqTransform")

#PCA of samples 
plotPCA(vsd, intgroup=c("Group", "OnPressors")) # The simplest PCA 

colnames(Kalantar_metadata)

interesting_metadata <- colnames(Kalantar_metadata)[-c(1,3,8:14,23,24)]
interesting_metadata

PCAdata <- plotPCA(vsd, intgroup= interesting_metadata, returnData = TRUE) 

PCAdata <- PCAdata %>% select(-group)
colnames(PCAdata) <- colnames(PCAdata) %>% janitor::make_clean_names()
interesting_metadata <- interesting_metadata %>% janitor::make_clean_names()

# ** formatted PCA plot ####
percentVar <- round(100 * attr(PCAdata, "percentVar")) 

ggplot(PCAdata, aes(x = pc1, y = pc2, color = group)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
#coord_fixed() 
  scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA_group.jpeg") 

ggplot(PCAdata, aes(x = pc1, y = pc2, color = on_pressors)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() 
  scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA_vasopressors.jpeg") 

ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = pc1, y = pc2, color = group)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ factor(on_pressors), nrow = 1, scales = "free") +
  theme(text = element_text(size = 14, family = "Calibri")) + 
  scale_color_brewer(palette="Set1", direction = -1) + 
  ggtitle("Are patients on vasopressors?")
# more not on pressors in the no sepsis group 

ggsave("PCA_vasopressors_group.jpeg") 

ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = pc1, y = pc2, color = factor(sirs_total))) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ factor(OnPressors), nrow = 1, scales = "free") +
  theme(text = element_text(size = 14, family = "Calibri")) + 
  scale_color_brewer(palette="Set1", direction = -1) + 
  ggtitle("Are patients on vasopressors?")
# more not on pressors in the no sepsis group 

ggsave("PCA_vasopressors_SIRS.jpeg") 

# PCA colored by different metadata #### 

colnames(Kalantar_metadata)

# Make important columns factors, change -/1 to yes/no 

plot_colored_PCA_discrete <- function(data, color){
  p <- ggplot(data = data) + 
    geom_point(mapping = aes(x = pc1, y = pc2, color = .data[[color]])) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    facet_wrap(~ on_pressors, nrow = 1, scales = "free") +
    theme(text = element_text(size = 14, family = "Calibri")) + 
    scale_color_brewer(palette="Set1", direction = -1) + 
    ggtitle("Are patients on vasopressors?")
  return(p)
}

plot_colored_PCA_cont <- function(data, color){
  p <- ggplot(data = data) + 
    geom_point(mapping = aes(x = pc1, y = pc2, color = .data[[color]])) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    facet_wrap(~ on_pressors, nrow = 1, scales = "free") +
    theme(text = element_text(size = 14, family = "Calibri")) + 
    scale_color_viridis_b() + 
    ggtitle("Are patients on vasopressors?")
  return(p)
}

# Tidy up column types in PCAdata for better plotting 
colnames(PCAdata)
interesting_metadata

# age, temp_max, apacheIII, wbc_max are continuous 
# rest should be factors 
# for continuous data make separate
interesting_metadata_discrete <- interesting_metadata[-c(2,6,8,9)]
interesting_metadata_discrete
PCAdata$x28d_death <- factor(PCAdata$x28d_death)
PCAdata$x28d_death
PCAdata$x28d_death <- ifelse(PCAdata$x28d_death == 0, "no", "yes")

PCAdata$intubated <- factor(PCAdata$intubated)
PCAdata$intubated <- ifelse(PCAdata$intubated == 0, "no", "yes")

PCAdata$on_pressors <- factor(PCAdata$on_pressors)
PCAdata$on_pressors <- ifelse(PCAdata$on_pressors == 0, "no", "yes")

PCAdata$immunocompromised <- factor(PCAdata$immunocompromised)
PCAdata$immunocompromised <- ifelse(PCAdata$immunocompromised == 0, "no", "yes")

PCAdata$sirs_total <- factor(PCAdata$sirs_total)

interesting_metadata_cont <- interesting_metadata[c(2,6,8,9)]
interesting_metadata_cont
PCAdata$age <- as.numeric(as.character(PCAdata$age))

## prints out the plots as in the loop
lapply(interesting_metadata_discrete, FUN = plot_colored_PCA_discrete, data = PCAdata) 
lapply(interesting_metadata_cont, FUN = plot_colored_PCA_cont, data = PCAdata) 

## save them to a list instead
plot_list_cont <- lapply(interesting_metadata_cont, FUN = plot_colored_PCA_cont, data = PCAdata) 
plot_list_discrete <- lapply(interesting_metadata_discrete, FUN = plot_colored_PCA_discrete, data = PCAdata) 

## view the first plot
plot_list_cont[[1]]
plot_list_discrete[[1]]

# most likely use function + lapply --> ggarrange 
# Or save the plots as a pdf - same as I did in my first RNASeq 

ggpubr::ggarrange(plotlist = plot_list_discrete)
ggsave("PCA_multi_color_discrete.jpeg", scale = 2)
# Saving 25.4 x 15.4 in image

ggpubr::ggarrange(plotlist = plot_list_cont)
ggsave("PCA_multi_color_cont.jpeg", scale = 2)

# Log transform some metadata #### 

interesting_metadata_cont
PCAdata

PCAdata <- PCAdata %>% mutate(apacheiii_log = log10(apacheiii + 1),
                              wbc_max_log = log10(wbc_max + 1))

plot_list_cont_log <- lapply(c("apacheiii_log", "wbc_max_log"), FUN = plot_colored_PCA_cont, data = PCAdata) 
ggpubr::ggarrange(plotlist = plot_list_cont_log)
ggsave("PCA_multi_color_cont_log.jpeg", scale = 2)
# Saving 25.1 x 9.66 in image


# Color by expression of specific genes #### 

genes_vasc_perm <- c("PECAM1", "ICAM1", "VCAM1", "DPP3", "ANGPT2", "IL6", "ALB", 
                     "SELE", "CRP", "SDC1", "ADM", "VEGF", "VWF", "CD142", "F2", 
                     "SERPINE1", "PROC", "THBD") 
# SELE also LYAM2 
# CD142 also F3 - tissue factor 

top_genes2
gene_names
genes_vasc_perm

# Biomart for vascular permeability genes 

gene_names_vasc_perm <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","description"),
                    values = genes_vasc_perm, mart= mart)
gene_names_vasc_perm # missing CD142 and VEGF 
write.csv(gene_names_vasc_perm, "genes_vascular_permeability.csv")

# Create a gene counts matrix with gene symbols 

Kalantar_counts_intersected

Kalantar_counts_fornames <- Kalantar_counts_intersected
Kalantar_counts_fornames$gene_id <- rownames(Kalantar_counts_fornames)

head(Kalantar_counts_fornames)

Kalantar_counts_fornames$gene_id <- gsub("\\.\\d{2,2}$", "", Kalantar_counts_fornames$gene_id)
Kalantar_counts_fornames$gene_id <- gsub("\\.\\d{1,1}$", "", Kalantar_counts_fornames$gene_id)


Kalantar_counts_fornames_gene_symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                    values = Kalantar_counts_fornames$gene_id, mart= mart)
head(Kalantar_counts_fornames_gene_symbols)
nrow(Kalantar_counts_fornames_gene_symbols)
nrow(Kalantar_counts_fornames)

Kalantar_counts_fornames <- full_join(Kalantar_counts_fornames, Kalantar_counts_fornames_gene_symbols, by = c("gene_id" = "ensembl_gene_id"))

# Subset for genes in vascular permeability and top PCA genes 
Kalantar_counts_vascperm <- Kalantar_counts_fornames %>% 
  filter(hgnc_symbol %in% genes_vasc_perm)

Kalantar_counts_loadings <- Kalantar_counts_fornames %>% 
  filter(hgnc_symbol %in% gene_names$hgnc_symbol)

# Merge with PCA data 

PCAdata$name # sample names 

# Transpose and clean - vasc perm 
Kalantar_counts_vascperm_tojoin <- t(Kalantar_counts_vascperm)
colnames(Kalantar_counts_vascperm_tojoin) <- Kalantar_counts_vascperm_tojoin["hgnc_symbol",]
Kalantar_counts_vascperm_tojoin <- as.data.frame(Kalantar_counts_vascperm_tojoin)
Kalantar_counts_vascperm_tojoin$name <- rownames(Kalantar_counts_vascperm_tojoin)

Kalantar_counts_vascperm_tojoin <- slice(Kalantar_counts_vascperm_tojoin, 1:(n() - 2))  
PCAdata_vascperm <- full_join(PCAdata, Kalantar_counts_vascperm_tojoin, by = "name")

# Transpose and clean - top loadings 
Kalantar_counts_loadings_tojoin <- t(Kalantar_counts_loadings)
colnames(Kalantar_counts_loadings_tojoin) <- Kalantar_counts_loadings_tojoin["hgnc_symbol",]
Kalantar_counts_loadings_tojoin <- as.data.frame(Kalantar_counts_loadings_tojoin)
Kalantar_counts_loadings_tojoin$name <- rownames(Kalantar_counts_loadings_tojoin)

Kalantar_counts_loadings_tojoin <- slice(Kalantar_counts_loadings_tojoin, 1:(n() - 2))  
PCAdata_loadings <- full_join(PCAdata, Kalantar_counts_loadings_tojoin, by = "name")

# Plot vascular permeability genes #### 

plot_list_vascperm <- lapply(genes_vasc_perm[-c(12,14)], FUN = plot_colored_PCA_cont, data = PCAdata_vascperm) 
plot_list_vascperm[[1]]
# ! Binned scales only support continuous data 
# Need to check the types of data the gene counts are! 
lapply(PCAdata_vascperm, typeof)

genes_vasc_perm[-c(12,14)] %in% colnames(PCAdata_vascperm)[19:length(colnames(PCAdata_vascperm))]
PCAdata_vascperm[, 19:ncol(PCAdata_vascperm)] <- lapply(19:ncol(PCAdata_vascperm), function(x) as.numeric(PCAdata_vascperm[[x]]))

ggpubr::ggarrange(plotlist = plot_list_vascperm)
ggsave("PCA_multi_color_vscperm.jpeg", scale = 2)
# Saving 27.3 x 15.4 in image

# Seems with raw counts there is some outlier stuff - normalise counts? merge with vst normalised counts? 

PCAdata_vascperm_log <- PCAdata_vascperm
PCAdata_vascperm_log[, 19:ncol(PCAdata_vascperm_log)] <- lapply(19:ncol(PCAdata_vascperm_log), function(x) log10(PCAdata_vascperm_log[[x]] + 1))

plot_list_vascperm_log <- lapply(genes_vasc_perm[-c(12,14)], FUN = plot_colored_PCA_cont, data = PCAdata_vascperm_log) 

ggpubr::ggarrange(plotlist = plot_list_vascperm_log)
ggsave("PCA_multi_color_vscperm_log.jpeg", width = 27.3, height = 15.4, units = "in")

# Plot top genes #### 

plot_list_loadings <- lapply(gene_names$hgnc_symbol, FUN = plot_colored_PCA_cont, data = PCAdata_loadings) 
plot_list_loadings[[1]]
# ! Binned scales only support continuous data 
# Need to check the types of data the gene counts are! 
lapply(PCAdata_loadings, typeof)

colnames(PCAdata_loadings)
gene_names$hgnc_symbol %in% colnames(PCAdata_loadings)[19:length(colnames(PCAdata_loadings))] # all TRUE
PCAdata_loadings[, 19:ncol(PCAdata_loadings)] <- lapply(19:ncol(PCAdata_loadings), function(x) as.numeric(PCAdata_loadings[[x]]))

ggpubr::ggarrange(plotlist = plot_list_loadings)
ggsave("PCA_multi_color_topload.jpeg", width = 27.3, height = 15.4, units = "in")
# Saving 27.3 x 15.4 in image

# Log transform counts 

PCAdata_loadings_log <- PCAdata_loadings
PCAdata_loadings_log[, 19:ncol(PCAdata_loadings_log)] <- lapply(19:ncol(PCAdata_loadings_log), function(x) log10(PCAdata_loadings_log[[x]] + 1))

plot_list_loadings_log <- lapply(gene_names$hgnc_symbol, FUN = plot_colored_PCA_cont, data = PCAdata_loadings_log) 

ggpubr::ggarrange(plotlist = plot_list_loadings_log)
ggsave("PCA_multi_color_topload_log.jpeg", width = 27.3, height = 15.4, units = "in")

# Plot unfaceted lists #### 

# New functions without facet 

plot_colored_PCA_discrete_nofacet <- function(data, color){
  p <- ggplot(data = data) + 
    geom_point(mapping = aes(x = pc1, y = pc2, color = .data[[color]])) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
 #   facet_wrap(~ on_pressors, nrow = 1, scales = "free") +
    theme(text = element_text(size = 14, family = "Calibri")) + 
    scale_color_brewer(palette="Set1", direction = -1) #+ 
  #  ggtitle("Are patients on vasopressors?")
  return(p)
}

plot_colored_PCA_cont_nofacet <- function(data, color){
  p <- ggplot(data = data) + 
    geom_point(mapping = aes(x = pc1, y = pc2, color = .data[[color]])) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
 #   facet_wrap(~ on_pressors, nrow = 1, scales = "free") +
    theme(text = element_text(size = 14, family = "Calibri")) + 
    scale_color_viridis_b() #+ 
  #  ggtitle("Are patients on vasopressors?")
  return(p)
}

# Metadata 
## save them to a list instead
plot_list_cont_nof <- lapply(interesting_metadata_cont, FUN = plot_colored_PCA_cont_nofacet, data = PCAdata) 
plot_list_discrete_nof <- lapply(interesting_metadata_discrete, FUN = plot_colored_PCA_discrete_nofacet, data = PCAdata) 

## view the first plot
plot_list_cont_nof[[1]]
lapply(PCAdata, typeof)
plot_list_discrete_nof[[1]]

ggpubr::ggarrange(plotlist = plot_list_cont_nof)
ggsave("PCA_multi_color_metadata_cont_nof.jpeg", scale = 1.2)
# Saving 8.01 x 5.81 in image

ggpubr::ggarrange(plotlist = plot_list_discrete_nof)
ggsave("PCA_multi_color_metadata_discrete_nof.jpeg", scale = 1.2)
# Saving 13.9 x 9.23 in image 

# Log transform some continuous variables 

plot_list_cont_log_nof <- lapply(c("apacheiii_log", "wbc_max_log"), FUN = plot_colored_PCA_cont_nofacet, data = PCAdata) 
ggpubr::ggarrange(plotlist = plot_list_cont_log_nof)
ggsave("PCA_metadata_cont_log_nof.jpeg")
# Saving 10.7 x 7.7 in image

# Vasc perm 

plot_list_vascperm_log_nof <- lapply(genes_vasc_perm[-c(12,14)], FUN = plot_colored_PCA_cont_nofacet, data = PCAdata_vascperm_log) 
plot_list_vascperm_log_nof[[1]]

ggpubr::ggarrange(plotlist = plot_list_vascperm_log_nof)
ggsave("PCA_multi_color_vascperm_log_nofacet.jpeg", scale = 1.1)
# Saving 12 x 8.46 in image

# top loadings 

plot_list_loadings_log_nof <- lapply(gene_names$hgnc_symbol, FUN = plot_colored_PCA_cont_nofacet, data = PCAdata_loadings_log) 

ggpubr::ggarrange(plotlist = plot_list_loadings_log_nof)
ggsave("PCA_multi_color_topload_log_nof.jpeg", scale = 1.2)

# NOTES #### 
# additional thing I could do it plot this above with a vst normalised gene counts matrix 
# Also could plot these genes without dividing on vasopressors 

# Regular prcomp with arrows #### 
# Could do this without scaling, using vst normalised data? Then should be the same 

# Apply PCA using prcomp function
# Need to scale / Normalize as PCA depends on distance measure

# Using vst 
# To get it the same as DeSeq2 I need to select 500 most variable genes from the object! 
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]

my_pca <- prcomp(t(assay(vsd)[select,]), scale = FALSE,
                 center = TRUE, retx = TRUE)

names(my_pca)

# Summary
summary(my_pca)
my_pca

# View the principal component loading
# my_pca$rotation[1:5, 1:4]
my_pca$rotation

# See the principal components
dim(my_pca$x)
my_pca$x

# Plotting the resultant principal components
# The parameter scale = 0 ensures that arrows
# are scaled to represent the loadings
biplot(my_pca, main = "Biplot", scale = 0, var.axes = FALSE)
# this plot makes no sense, you can't see nothing 

# Compute standard deviation
my_pca$sdev

# For Scree plot #### 
# Compute variance 
my_pca.var <- my_pca$sdev ^ 2
my_pca.var

# Proportion of variance for a scree plot
propve <- my_pca.var / sum(my_pca.var)
propve

# Plot variance explained for each principal component
jpeg(file="scree_plot_prcomp_vsd.jpeg")
plot(propve[1:30], xlab = "principal component",
     ylab = "Proportion of Variance Explained",
     ylim = c(0, 1), type = "b",
     main = "Scree Plot - first 30 PCs")
dev.off()

# Plot the cumulative proportion of variance explained
plot(cumsum(propve),
     xlab = "Principal Component",
     ylab = "Cumulative Proportion of Variance Explained",
     ylim = c(0, 1), type = "b")

# Find Top n principal component
# which will atleast cover 90 % variance of dimension
which(cumsum(propve) >= 0.9)[1] # 56, with vsd 110, with vsd and top 500 now 36

# Plot with circles #### 
head(my_pca)

??ggbiplot
ggbiplot(my_pca, 
         obs.scale = 1,
         var.scale = 1, 
         groups = Kalantar_metadata$Group, # would need to be ordered the same? same for the other pca? 
         ellipse = TRUE, 
         circle = FALSE, 
         ellipse.prob = 0.68, 
         var.axes = FALSE) + 
  scale_color_discrete(name = '') + 
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggsave("prcomp_ellipses_vsd.jpeg")

ggbiplot(my_pca, 
         obs.scale = 1,
         var.scale = 1, 
         groups = factor(Kalantar_metadata$on_pressors), # would need to be ordered the same? same for the other pca? 
         ellipse = TRUE, 
         circle = FALSE, 
         ellipse.prob = 0.68, 
         var.axes = FALSE) + 
  scale_color_discrete(name = '') + 
  theme(legend.direction = 'horizontal', legend.position = 'top')

ggsave("prcomp_ellipses_vasopressors_vsd.jpeg")

# Other plots for prcomp #### 

# ** Variance (scree +) ####
pc_eigenvalues <- my_pca$sdev^2 

# get ready for plotting 
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                         variance = pc_eigenvalues) %>% 
  mutate(pct = variance/sum(variance)*100) %>% 
  mutate(pct_cum = cumsum(pct))

pc_eigenvalues

pc_eigenvalues[1:30,] %>% 
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) +
  labs(x = "Principal component", y = "Fraction variance explained")

ggsave("scree_ggplot_vsd.jpeg")

# ** Samples in PC space with ggplot #### 

pc_scores <- my_pca$x 

pc_scores <- pc_scores %>% 
  as_tibble(rownames = "ID2")

pc_scores
max(pc_scores$PC1)

pc_plot <- pc_scores %>% 
  full_join(Kalantar_metadata, by = "ID2") %>% 
  ggplot(aes(x = PC1, y = PC2, color = Group)) +
  geom_point()


# ** Loadings #### 

pc_loadings <- my_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

pc_loadings

# top genes with highest loadings on PC1 and PC2 
top_genes <- pc_loadings %>% 
  select(gene, PC1, PC2) %>%
  tidyr::pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
  group_by(PC) %>% 
  arrange(desc(abs(loading))) %>% 
  slice(1:10) %>% 
  pull(gene) %>% 
  unique()

top_genes

top_loadings <- pc_loadings %>% 
  filter(gene %in% top_genes)

top_loadings

# *** Biomart Explain top genes ####

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

top_genes2 <- gsub("\\.\\d{2,2}$", "", top_genes)
top_genes2 <- gsub("\\.\\d{1,1}$", "", top_genes2)
top_genes2

# ** gene names #### 
gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                    values = top_genes2, mart= mart)
gene_names
# with vsd object the top loadings changed! 
write.csv(gene_names, "top_loadings_genes_vsd.csv")

# join gene names 
top_loadings <- cbind(top_loadings, gene_names)
top_loadings

# plot loadings 
loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "brown") +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = hgnc_symbol),
            nudge_x = -0.005, nudge_y = 0, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))


# ** Arrange PCA and loadings together ####

ggpubr::ggarrange(pc_plot, loadings_plot, ncol = 2, widths = c(1.5,1))

ggsave("pca_loadings_vsd.jpeg")
# Saving 11 x 4.84 in image
# Saving 12.7 x 6.58 in image 