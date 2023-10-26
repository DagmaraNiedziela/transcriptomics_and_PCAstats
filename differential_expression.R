# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readxl)
library("biomaRt")
library(apeglm)
library(org.Hs.eg.db)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
source("functions_PCA-transcriptomics.R")

# Data 
# Load data #### 

# Counts and metadata loaded in PCA_DESeq.R 
head(Kalantar_counts)
head(Kalantar_metadata)

# Remove Sepsis.suspect and sepsis.indeterm groups - I don't need them 
# Remove samples with specific group in counts 
no_weird_samples <- Kalantar_metadata$id2[Kalantar_metadata$group %in% c("No-Sepsis", "SepsisBSI", "SepsisNon-BSI")]

Kalantar_counts_noweird <- Kalantar_counts %>% 
  dplyr::select(all_of(no_weird_samples))
ncol(Kalantar_counts_noweird) # 110

# Remove samples in coldata 

Kalantar_metadata_noweird <- Kalantar_metadata %>% filter(group %in% c("No-Sepsis", "SepsisBSI", "SepsisNon-BSI")) 
nrow(Kalantar_metadata_noweird) # 110

# check order of samples 
colnames(Kalantar_counts_noweird)
Kalantar_metadata_noweird$id2
colnames(Kalantar_counts_noweird) == Kalantar_metadata_noweird$id2 # yes 

# Load DESeq object 
dds_vaso <- DESeqDataSetFromMatrix(countData = Kalantar_counts_noweird,
                              colData = Kalantar_metadata_noweird, 
                              design = ~ on_pressors + group) 

nrow(dds_vaso) #19939 

# Outliers remove 

keep <- rowSums(counts(dds_vaso)) >= 10
dds_vaso <- dds_vaso[keep,]
nrow(dds_vaso) 
#[1] 15054  

# Quick PCA without weird groups ####

# Normalise data with a vst function 
vsd_vaso <- vst(dds_vaso, blind=FALSE)

#PCA of samples 
plotPCA(vsd_vaso, intgroup=c("group", "on_pressors")) # The simplest PCA 
plotPCA(vsd_vaso, intgroup="group")
plotPCA(vsd_vaso, intgroup="on_pressors")

?plyr::rename
PCAdata_vaso <- plotPCA(vsd_vaso, intgroup= colnames(Kalantar_metadata_noweird), returnData = TRUE) 
PCAdata_vaso <- PCAdata_vaso %>% dplyr::select(-group) 
PCAdata_vaso <- PCAdata_vaso %>% plyr::rename(c("group.1" = "group"))

# PCAdata <- PCAdata %>% select(-group)

# ** formatted PCA plot ####
percentVar <- round(100 * attr(PCAdata_vaso, "percentVar")) 

ggplot(PCAdata_vaso, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size =3) + 
  stat_ellipse() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() 
  scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA/PCA_group_lessgroups.jpeg") 
# Saving 4.98 x 3.55 in image

ggplot(PCAdata_vaso, aes(x = PC1, y = PC2, color = on_pressors)) +
  geom_point(size =3) + 
  stat_ellipse() +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() 
  scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA/PCA_vasopressors_lessgroups.jpeg") 

# Differential expression - vasopressors, sepsis only #### 

Kalantar_sepsisonly <- subset_counts_and_metadata(counts = Kalantar_counts,
                                                metadata = Kalantar_metadata,
                                                groups = c("SepsisBSI", "SepsisNon-BSI"))

colnames(Kalantar_sepsisonly$newcounts) == Kalantar_sepsisonly$newmetadata$id2 # yes

# Object 
dds_vaso_sep <- DESeqDataSetFromMatrix(countData = Kalantar_sepsisonly$newcounts,
                                   colData = Kalantar_sepsisonly$newmetadata, 
                                   design = ~ on_pressors) 

nrow(dds_vaso_sep) #19939 

# Outliers remove 
keep <- rowSums(counts(dds_vaso_sep)) >= 10
dds_vaso_sep <- dds_vaso_sep[keep,]
nrow(dds_vaso_sep) 
#[1] 14216  

# DESeq 

dds_vaso_sep <- estimateSizeFactors(dds_vaso_sep) 
dds_vaso_sep <- DESeq(dds_vaso_sep)

resultsNames(dds_vaso_sep) 
results_vaso_sep <- results(dds_vaso_sep, name="on_pressors_yes_vs_no")

# MA plot 
jpeg(file="differential_expression/MAplot_sepsisonly.jpeg", 
     width = 15, height = 10, units = "cm", res = 300)
plotMA(results_vaso_sep, ylim=c(-2,2))
dev.off()

plotMA(results_vaso_sep_shrink, ylim=c(-2,2))

# Shrinkage 

results_vaso_sep_shrink <- lfcShrink(dds_vaso_sep, coef="on_pressors_yes_vs_no", type="apeglm")
# will use these 

results_vaso_sep_shrink

# Annotate and export results 

results_vaso_nosep_shrink <- annotate_genes(results_vaso_sep_shrink)

# Check DE gene numbers 

sum(results_vaso_sep_shrink$padj < 0.4, na.rm=TRUE) # 13

write.csv(as.data.frame(results_vaso_sep_shrink), 
          "differential_expression/DEG_all_vasopressors_sepsisonly.csv")

writexl::write_xlsx(list(sepsis_all = as.data.frame(results_vaso_sep_shrink), 
                         sepsis_sig = as.data.frame(subset(results_vaso_sep_shrink, padj < 0.05))), 
                    "differential_expression/Kalantar_plasma_differential_expression_sepsis_only.xlsx")

# Top genes 

res_top10_sep <- make_top_genes(results_vaso_sep_shrink, rep_name = "nosepsis")
res_top10_sep

# Plot top genes 

ggplot(data = res_top10_sep) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - sepsis only")

ggsave("differential_expression/top_10_genes_sepsisonly.jpeg") 
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 


# Differential expression - vasopressors, no sepsis group #### 

Kalantar_nosepsis <- subset_counts_and_metadata(counts = Kalantar_counts,
                                                  metadata = Kalantar_metadata,
                                                  groups = "No-Sepsis")

colnames(Kalantar_nosepsis$newcounts) == Kalantar_nosepsis$newmetadata$id2 # yes 

# Object 
dds_vaso_nosep <- DESeqDataSetFromMatrix(countData = Kalantar_nosepsis$newcounts,
                                       colData = Kalantar_nosepsis$newmetadata, 
                                       design = ~ on_pressors) 

nrow(dds_vaso_nosep) #27097 

# Outliers remove 
keep <- rowSums(counts(dds_vaso_nosep)) >= 10
dds_vaso_nosep <- dds_vaso_nosep[keep,]
nrow(dds_vaso_nosep) 
#[1] 12833  

# DESeq 

dds_vaso_nosep <- estimateSizeFactors(dds_vaso_nosep) 
dds_vaso_nosep <- DESeq(dds_vaso_nosep)

resultsNames(dds_vaso_nosep) 
results_vaso_nosep <- results(dds_vaso_nosep, name="on_pressors_yes_vs_no")

# MA plot 
plotMA(results_vaso_nosep, ylim=c(-2,2))

plotMA(results_vaso_nosep_shrink, ylim=c(-2,2))

# Shrinkage 

results_vaso_nosep_shrink <- lfcShrink(dds_vaso_nosep, coef="on_pressors_yes_vs_no", type="apeglm")
# will use these 

results_vaso_nosep_shrink

# Annotate and export results 

results_vaso_nosep_shrink <- annotate_genes(results_vaso_nosep_shrink)

# Check DE gene numbers 

sum(results_vaso_nosep_shrink$padj < 0.4, na.rm=TRUE) # 13

write.csv(as.data.frame(results_vaso_nosep_shrink), 
          "differential_expression/DEG_all_vasopressors_no_sepsis.csv")

writexl::write_xlsx(list(nosepsis_all = as.data.frame(results_vaso_nosep_shrink), 
                         nosepsis_sig = as.data.frame(subset(results_vaso_nosep_shrink, padj < 0.05))), 
                    "Kalantar_differential_expression_nosepsis.xlsx")

# Top genes 

res_top10_nosep <- make_top_genes(results_vaso_nosep_shrink, rep_name = "nosepsis")
res_top10_nosep

# Plot top genes 

ggplot(data = res_top10_nosep) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - no sepsis")

ggsave("top_10_genes_nosep.jpeg") 
# Saving 6.65 x 5.15 in image
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 


# NOT DONE #### 
# Differential expression, Group  #### 

# Clean up metadata 

dds_group <- DESeqDataSetFromMatrix(countData = Kalantar_counts_intersected,
                                    colData = Kalantar_metadata, 
                                    design = ~ group) 

nrow(dds_group) #27097 

# Outliers remove 
keep <- rowSums(counts(dds_group)) >= 10
dds_group <- dds_group[keep,]
nrow(dds_group) 
#[1] 27065 

# DESeq 

dds_group <- estimateSizeFactors(dds_group) 
dds_group <- DESeq(dds_group)

resultsNames(dds_group) 
results_group_bsivsnonbsi <- results(dds_group, contrast=c("group","SepsisBSI","SepsisNon-BSI"))
results_group_bsivsnosep <- results(dds_group, name="group_SepsisBSI_vs_No.Sepsis")
results_group_nonbsivsnosep <- results(dds_group, name="group_SepsisNon.BSI_vs_No.Sepsis")

# MA plot 
plotMA(results_group_bsivsnonbsi, ylim=c(-2,2))
plotMA(results_group_bsivsnosep, ylim=c(-2,2))
plotMA(results_group_nonbsivsnosep, ylim=c(-2,2))

plotMA(results_group_bsivsnonbsi_shrink, ylim=c(-2,2))
plotMA(results_group_bsivsnosep_shrink, ylim=c(-2,2))
plotMA(results_group_nonbsivsnosep_shrink, ylim=c(-2,2))

# Shrinkage 

results_group_bsivsnonbsi_shrink <- lfcShrink(dds_group, contrast = c("group", "SepsisBSI", "SepsisNon-BSI"), type="ashr")
results_group_bsivsnosep_shrink <- lfcShrink(dds_group, coef="group_SepsisBSI_vs_No.Sepsis", type="ashr")
results_group_nonbsivsnosep_shrink <- lfcShrink(dds_group, coef="group_SepsisNon.BSI_vs_No.Sepsis", type="ashr")

# Annotate and export results 

results_group_bsivsnonbsi_shrink <- annotate_genes(results_group_bsivsnonbsi_shrink)
results_group_bsivsnosep_shrink <- annotate_genes(results_group_bsivsnosep_shrink)
results_group_nonbsivsnosep_shrink <- annotate_genes(results_group_nonbsivsnosep_shrink)
# results_group_bsivsnonbsi_shrink

# Check DE gene numbers 

sum(results_group_bsivsnonbsi_shrink$padj < 0.05, na.rm=TRUE) # 3757
sum(results_group_bsivsnosep_shrink$padj < 0.05, na.rm=TRUE) # 6162 
sum(results_group_nonbsivsnosep_shrink$padj < 0.05, na.rm=TRUE) # 1333 

writexl::write_xlsx(list(bsi_nonbsi_all = as.data.frame(results_group_bsivsnonbsi_shrink), 
                         bsi_nonbsi_sig = as.data.frame(subset(results_group_bsivsnonbsi_shrink, padj < 0.05)),
                         bsi_nosep_all = as.data.frame(results_group_bsivsnosep_shrink), 
                         bsi_nosep_sig = as.data.frame(subset(results_group_bsivsnosep_shrink, padj < 0.05)),
                         nonbsi_nosep_all = as.data.frame(results_group_nonbsivsnosep_shrink), 
                         nonbsi_nosep_sig = as.data.frame(subset(results_group_nonbsivsnosep_shrink, padj < 0.05))), 
                    "Kalantar_differential_expression_group.xlsx")

# Top genes 

res_top10_group1 <- make_top_genes(results_group_bsivsnonbsi_shrink, rep_name = "bsi_nonbsi")
res_top10_group2 <- make_top_genes(results_group_bsivsnosep_shrink, rep_name = "bsi_nosep")
res_top10_group3 <- make_top_genes(results_group_nonbsivsnosep_shrink, rep_name = "nonbsi_nosep")
res_top10_group <- rbind(res_top10_group1, res_top10_group2, res_top10_group3)

# Plot top genes 

ggplot(data = res_top10_group) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - all groups")

ggsave("top_10_genes_group.jpeg") 
# Saving 6.65 x 5.15 in image
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 