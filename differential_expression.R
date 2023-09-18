# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readxl)
library("biomaRt")
# BiocManager::install("apeglm")
library(apeglm)
# BiocManager::install("org.Hs.eg.db", update = FALSE)
library(org.Hs.eg.db)

# Data 
# Load data #### 

Kalantar_counts <- read.csv("GSE189400_earli_paxgene_counts_50k.csv", header = TRUE)
head(Kalantar_counts)

Kalantar_metadata <- read_excel("Kalantar_supplementary.xlsx", sheet = "Data_16")
head(Kalantar_metadata)

ncol(Kalantar_counts) # 379
nrow(Kalantar_metadata) # 221 

Kalantar_counts_reduced <- Kalantar_counts %>% select(-ends_with(".1"))
length(unique(colnames(Kalantar_counts_reduced))) # 348 

# Need to split EARLI from colnames, or add EARLI to ID 
Kalantar_metadata$ID2 <- paste("EARLI", Kalantar_metadata$ID, sep = "_")
# intersect(Kalantar_metadata$ID2, colnames(Kalantar_counts_reduced))

# Final counts only contained in metadata 
Kalantar_counts_intersected <- Kalantar_counts_reduced %>% select(X, Kalantar_metadata$ID2) 
# This also makes sure the counts are ordered the same way as metadata !!! 


# Clean up metadata #### 

Kalantar_metadata <- Kalantar_metadata %>% janitor::clean_names()

Kalantar_metadata$x28d_death <- ifelse(Kalantar_metadata$x28d_death == 0, "no", "yes")
Kalantar_metadata$intubated <- ifelse(Kalantar_metadata$intubated == 0, "no", "yes")
Kalantar_metadata$on_pressors <- ifelse(Kalantar_metadata$on_pressors == 0, "no", "yes")
Kalantar_metadata$immunocompromised <- ifelse(Kalantar_metadata$immunocompromised == 0, "no", "yes")
Kalantar_metadata$sirs_total <- factor(Kalantar_metadata$sirs_total)
Kalantar_metadata$age <- as.numeric(as.character(Kalantar_metadata$age))

# DESeq object - Group #### 

rownames(Kalantar_counts_intersected) <- Kalantar_counts_intersected$X
Kalantar_counts_intersected <- Kalantar_counts_intersected[,-1]

dds <- DESeqDataSetFromMatrix(countData = Kalantar_counts_intersected,
                              colData = Kalantar_metadata, 
                              design = ~ Group) 

nrow(dds) #27097 

# Outliers remove 

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
nrow(dds) 
#[1] 27065 

# Protein coding genes - maybe no need, do without for now ####

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

listAttributes(mart)
# rownames(assay(dds)) - cannot easily modify, will make new dds 
rownames(Kalantar_counts_intersected)

rownames(Kalantar_counts_intersected) <- gsub("\\.\\d{2,2}$", "", rownames(Kalantar_counts_intersected))
rownames(Kalantar_counts_intersected) <- gsub("\\.\\d{1,1}$", "", rownames(Kalantar_counts_intersected))
rownames(Kalantar_counts_intersected)

gene_biotypes <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","gene_biotype"),
                    values = top_genes2, mart= mart)


# Differential expression - vasopressors only #### 

# Clean up metadata 

dds_vaso <- DESeqDataSetFromMatrix(countData = Kalantar_counts_intersected,
                              colData = Kalantar_metadata, 
                              design = ~ on_pressors) 

nrow(dds_vaso) #27097 

# Outliers remove 
keep <- rowSums(counts(dds_vaso)) >= 10
dds_vaso <- dds_vaso[keep,]
nrow(dds_vaso) 
#[1] 27065 

# DESeq 

dds_vaso <- estimateSizeFactors(dds_vaso) 
dds_vaso <- DESeq(dds_vaso)

resultsNames(dds_vaso) 
results_vaso <- results(dds_vaso, name="on_pressors_yes_vs_no")

# MA plot 
plotMA(results_vaso, ylim=c(-2,2))

plotMA(results_vaso_shrink, ylim=c(-2,2))

# Shrinkage 

results_vaso_shrink <- lfcShrink(dds_vaso, coef="on_pressors_yes_vs_no", type="apeglm")
# will use these 

results_vaso_shrink

# Annotate and export results 

rownames(results_vaso_shrink) <- gsub("\\.\\d{2,2}$", "", rownames(results_vaso_shrink))
rownames(results_vaso_shrink) <- gsub("\\.\\d{1,1}$", "", rownames(results_vaso_shrink))
rownames(results_vaso_shrink)

results_vaso_shrink$gene_symbol <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = rownames(results_vaso_shrink),
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

results_vaso_shrink$gene_name <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = rownames(results_vaso_shrink),
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "GENENAME", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

write.csv(results_vaso_shrink, "DEG_all_vasopressors.csv")

# Check DE gene numbers 

sum(results_vaso$padj < 0.05, na.rm=TRUE) # 2410 
results_vaso_shrink_DEG <- subset(results_vaso_shrink, padj < 0.05)
results_vaso_shrink_DEG

write.csv(results_vaso_shrink_DEG, "DEG_sig_vasopressors.csv")

# Top genes 

head(as.data.frame(results_vaso_shrink_DEG))

res_vaso_top10_up <- as.data.frame(results_vaso_shrink_DEG) %>% filter(log2FoldChange > 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
  tibble::add_column(condition = (rep("On_vasopressors_yes-no",10)), Direction = rep("Up", 10)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_vaso_top10_down <- as.data.frame(results_vaso_shrink_DEG) %>% filter(log2FoldChange < 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
  tibble::add_column(condition = (rep("On_vasopressors_yes-no",10)), Direction = rep("Down", 10))

res_top10_vaso <- rbind(res_vaso_top10_up, res_vaso_top10_down)
res_top10_vaso <- res_top10_vaso %>%
  mutate(gene_symbol2 = coalesce(gene_symbol, rownames(.)))
res_top10_vaso 
res_top10_vaso$log2FoldChange <- as.numeric(as.character(res_top10_vaso$log2FoldChange))

# Plot top genes 

ggplot(data = res_top10_vaso) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - all groups")

ggsave("top_10_genes_vaso.jpeg") 
# Saving 6.65 x 5.15 in image
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 

# Differential expression - vasopressors, remove no sepsis group #### 

# Remove samples with specific group in counts 
no_sepsis_samples <- Kalantar_metadata$id2[Kalantar_metadata$group == "No-Sepsis"]

Kalantar_counts_sepsis <- Kalantar_counts_intersected %>% 
  dplyr::select(-all_of(no_sepsis_samples))
ncol(Kalantar_counts_sepsis) # 129 

# Remove samples in coldata 

Kalantar_metadata_sepsis <- Kalantar_metadata %>% filter(group != "No-Sepsis") 
nrow(Kalantar_metadata_sepsis) # 129 

# check order of samples 
colnames(Kalantar_counts_sepsis)
Kalantar_metadata_sepsis$id2
colnames(Kalantar_counts_sepsis) == Kalantar_metadata_sepsis$id2

# Object 
dds_vaso_sep <- DESeqDataSetFromMatrix(countData = Kalantar_counts_sepsis,
                                   colData = Kalantar_metadata_sepsis, 
                                   design = ~ on_pressors) 

nrow(dds_vaso_sep) #27097 

# Outliers remove 
keep <- rowSums(counts(dds_vaso_sep)) >= 10
dds_vaso_sep <- dds_vaso_sep[keep,]
nrow(dds_vaso_sep) 
#[1] 27009  

# DESeq 

dds_vaso_sep <- estimateSizeFactors(dds_vaso_sep) 
dds_vaso_sep <- DESeq(dds_vaso_sep)

resultsNames(dds_vaso_sep) 
results_vaso_sep <- results(dds_vaso_sep, name="on_pressors_yes_vs_no")

# MA plot 
plotMA(results_vaso_sep, ylim=c(-2,2))

plotMA(results_vaso_sep_shrink, ylim=c(-2,2))

# Shrinkage 

results_vaso_sep_shrink <- lfcShrink(dds_vaso_sep, coef="on_pressors_yes_vs_no", type="apeglm")
# will use these 

results_vaso_sep_shrink

# Annotate and export results 

rownames(results_vaso_sep_shrink) <- gsub("\\.\\d{2,2}$", "", rownames(results_vaso_sep_shrink))
rownames(results_vaso_sep_shrink) <- gsub("\\.\\d{1,1}$", "", rownames(results_vaso_sep_shrink))
rownames(results_vaso_sep_shrink)

results_vaso_sep_shrink$gene_symbol <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = rownames(results_vaso_sep_shrink),
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

results_vaso_sep_shrink$gene_name <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = rownames(results_vaso_sep_shrink),
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "GENENAME", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

write.csv(results_vaso_sep_shrink, "DEG_all_vasopressors_sepsis_only.csv")

# Check DE gene numbers 

sum(results_vaso_sep$padj < 0.05, na.rm=TRUE) # 1158 
results_vaso_sep_shrink_DEG <- subset(results_vaso_sep_shrink, padj < 0.05)
results_vaso_sep_shrink_DEG

write.csv(results_vaso_sep_shrink_DEG, "DEG_sig_vasopressors_sepsis_only.csv")

# Top genes 

head(as.data.frame(results_vaso_sep_shrink_DEG))

res_vaso_sep_top10_up <- as.data.frame(results_vaso_sep_shrink_DEG) %>% filter(log2FoldChange > 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
  tibble::add_column(condition = (rep("On_vasopressors_yes-no",11)), Direction = rep("Up", 11)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_vaso_sep_top10_down <- as.data.frame(results_vaso_sep_shrink_DEG) %>% filter(log2FoldChange < 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
  tibble::add_column(condition = (rep("On_vasopressors_yes-no",10)), Direction = rep("Down", 10))

res_top10_vaso_sep <- rbind(res_vaso_sep_top10_up, res_vaso_sep_top10_down)
res_top10_vaso_sep <- res_top10_vaso_sep %>%
  mutate(gene_symbol2 = coalesce(gene_symbol, rownames(.)))
res_top10_vaso_sep 
res_top10_vaso_sep$log2FoldChange <- as.numeric(as.character(res_top10_vaso_sep$log2FoldChange))

# Plot top genes 

ggplot(data = res_top10_vaso_sep) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - sepsis only")

ggsave("top_10_genes_vaso_sep.jpeg") 
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 


# ** Differential expression, Group x Vasopressors #### 



