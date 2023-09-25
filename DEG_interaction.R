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

# Differential expression - vasopressors only #### 

# Clean up metadata 

dds_vaso_group <- DESeqDataSetFromMatrix(countData = Kalantar_counts_intersected,
                                   colData = Kalantar_metadata, 
                                   design = ~ group * on_pressors) # issue with group having - sign 

nrow(dds_vaso_group) #27097 

# Outliers remove 
keep <- rowSums(counts(dds_vaso_group)) >= 10
dds_vaso_group <- dds_vaso_group[keep,]
nrow(dds_vaso_group) 
#[1] 27065 

# DESeq 

dds_vaso_group <- estimateSizeFactors(dds_vaso_group) 
dds_vaso_group <- DESeq(dds_vaso_group)

resultsNames(dds_vaso_group) 

# for Group * on_pressors 
# [1] "Intercept"                         "group_SepsisBSI_vs_No.Sepsis"      "group_SepsisNon.BSI_vs_No.Sepsis" 
# [4] "on_pressors_yes_vs_no"             "groupSepsisBSI.on_pressorsyes"     "groupSepsisNon.BSI.on_pressorsyes"

# The condition effect for on pressors for No-Sepsis 
results_vaso_no_sepsis <- results(dds_vaso_group, name="on_pressors_yes_vs_no") # for no sepsis? 
sum(results_vaso_no_sepsis$padj < 0.05, na.rm=TRUE) # 16 

# The condition (on_pressors) effect for genotype (group 2)
# This is listed as the "extra" effect in group 3 compared to group 1 - noooo 
results_vaso_sepsis_bsi <- results(dds_vaso_group, contrast=list( c("on_pressors_yes_vs_no","groupSepsisBSI.on_pressorsyes") ))
sum(results_vaso_sepsis_bsi$padj < 0.05, na.rm=TRUE) # 0 

# The condition (on_pressors) effect for genotype (group 3)
results_vaso_sepsis_nonbsi <- results(dds_vaso_group, contrast=list( c("on_pressors_yes_vs_no","groupSepsisNon.BSI.on_pressorsyes") ))
sum(results_vaso_sepsis_nonbsi$padj < 0.05, na.rm=TRUE) # 121 
head(results_vaso_sepsis_nonbsi)


# Use the analysis below, results are the same, but model looks easier 

# Top genes 

head(as.data.frame(results_vaso_group_shrink_DEG))

####### fix add_column here, wrong condition
res_vaso_group_top10_up <- as.data.frame(results_vaso_group_shrink_DEG) %>% filter(log2FoldChange > 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
  tibble::add_column(condition = (rep("On_vaso_grouppressors_yes-no",10)), Direction = rep("Up", 10)) 
# AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
res_vaso_group_top10_down <- as.data.frame(results_vaso_group_shrink_DEG) %>% filter(log2FoldChange < 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
  tibble::add_column(condition = (rep("On_vaso_grouppressors_yes-no",10)), Direction = rep("Down", 10))

res_top10_vaso_group <- rbind(res_vaso_group_top10_up, res_vaso_group_top10_down)
res_top10_vaso_group <- res_top10_vaso_group %>%
  mutate(gene_symbol2 = coalesce(gene_symbol, rownames(.)))
res_top10_vaso_group 
res_top10_vaso_group$log2FoldChange <- as.numeric(as.character(res_top10_vaso_group$log2FoldChange))

# Plot top genes 

ggplot(data = res_top10_vaso_group) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - all groups")

ggsave("top_10_genes_vaso_group.jpeg") 
# Saving 6.65 x 5.15 in image
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 

# DEG Condition - new column for group & vasopressors #### 

# Make a condition column 
Kalantar_metadata
Kalantar_metadata$condition <- paste0(Kalantar_metadata$group,"_",Kalantar_metadata$on_pressors)
unique(Kalantar_metadata$condition)

# Deseq object 
dds_condition <- DESeqDataSetFromMatrix(countData = Kalantar_counts_intersected,
                                         colData = Kalantar_metadata, 
                                         design = ~ condition) # issue with group having - sign 

nrow(dds_condition) #27097 

# Outliers remove 
keep <- rowSums(counts(dds_condition)) >= 10
dds_condition <- dds_condition[keep,]
nrow(dds_condition) 
#[1] 27065 

# DESeq 

dds_condition <- estimateSizeFactors(dds_condition) 
dds_condition <- DESeq(dds_condition)

resultsNames(dds_condition) 

# [1] "Intercept"                                   "condition_No.Sepsis_yes_vs_No.Sepsis_no"    
# [3] "condition_SepsisBSI_no_vs_No.Sepsis_no"      "condition_SepsisBSI_yes_vs_No.Sepsis_no"    
# [5] "condition_SepsisNon.BSI_no_vs_No.Sepsis_no"  "condition_SepsisNon.BSI_yes_vs_No.Sepsis_no"


# Get results per group 
results_vaso_sepsis_bsi_cond <- results(dds_condition, contrast=c("condition","SepsisBSI_yes","SepsisBSI_no"))
sum(results_vaso_sepsis_bsi_cond$padj < 0.05, na.rm=TRUE) # 0 

results_vaso_sepsis_nonbsi_cond <- results(dds_condition, contrast=c("condition","SepsisNon-BSI_yes","SepsisNon-BSI_no"))
sum(results_vaso_sepsis_nonbsi_cond$padj < 0.05, na.rm=TRUE) #121 
head(results_vaso_sepsis_nonbsi_cond)

results_vaso_no_sepsis_cond <- results(dds_condition, name="condition_No.Sepsis_yes_vs_No.Sepsis_no")
sum(results_vaso_no_sepsis_cond$padj < 0.05, na.rm=TRUE) # 16 
# The same as the interaction results!!!!! 

# MA plots 
plotMA(results_vaso_sepsis_bsi_cond, ylim=c(-2,2))
plotMA(results_vaso_sepsis_nonbsi_cond, ylim=c(-2,2))
plotMA(results_vaso_no_sepsis_cond, ylim=c(-2,2))

# plotMA(results_vaso_sepsis_bsi_cond_shrink, ylim=c(-2,2)) - no need, no genes 
plotMA(results_vaso_sepsis_nonbsi_cond_shrink, ylim=c(-2,2))
plotMA(results_vaso_no_sepsis_cond_shrink, ylim=c(-2,2))

# Shrink results 
# results_vaso_sepsis_bsi_cond_shrink <- lfcShrink(dds_condition, contrast = c("condition", "SepsisBSI_yes", "SepsisBSI_no"), type="ashr") # no genes 
# install.packages("ashr")
results_vaso_sepsis_nonbsi_cond_shrink <- lfcShrink(dds_condition, contrast = c("condition", "SepsisNon-BSI_yes", "SepsisNon-BSI_no"), type="ashr")
results_vaso_no_sepsis_cond_shrink <- lfcShrink(dds_condition, coef="condition_No.Sepsis_yes_vs_No.Sepsis_no", type="ashr")


# Annotate genes 

annotate_genes <- function(results_shrunk){
  rownames(results_shrunk) <- gsub("\\.\\d{2,2}$", "", rownames(results_shrunk))
  rownames(results_shrunk) <- gsub("\\.\\d{1,1}$", "", rownames(results_shrunk))

  results_shrunk$gene_symbol <- mapIds(
    org.Hs.eg.db, # Replace with annotation package for your organism
    keys = rownames(results_shrunk),
    keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
    column = "SYMBOL", # The type of gene identifiers you would like to map to
    multiVals = "first"
  )

  results_shrunk$gene_name <- mapIds(
    org.Hs.eg.db, # Replace with annotation package for your organism
    keys = rownames(results_shrunk),
    keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
    column = "GENENAME", # The type of gene identifiers you would like to map to
    multiVals = "first"
  )
  
  results_shrunk$ensembl_id <- rownames(results_shrunk)
  return(results_shrunk)
  
}

results_vaso_sepsis_nonbsi_cond_shrink <- annotate_genes(results_vaso_sepsis_nonbsi_cond_shrink)
results_vaso_no_sepsis_cond_shrink <- annotate_genes(results_vaso_no_sepsis_cond_shrink)

results_vaso_sepsis_nonbsi_cond_shrink

# Save files 

# write.csv(results_vaso_sepsis_nonbsi_cond_shrink, "DEG_conditions_sepsis_nonbsi.csv")
# 
# write.csv(subset(results_vaso_sepsis_nonbsi_cond_shrink, padj < 0.05), "DEG_conditions_sepsis_nonbsi_sig.csv")

writexl::write_xlsx(list(nonbsi_all = as.data.frame(results_vaso_sepsis_nonbsi_cond_shrink), 
                         nonbsi_sig = as.data.frame(subset(results_vaso_sepsis_nonbsi_cond_shrink, padj < 0.05)),
                         nosepsis_all = as.data.frame(results_vaso_no_sepsis_cond_shrink),
                         nosepsis_sig = as.data.frame(subset(results_vaso_no_sepsis_cond_shrink, padj < 0.05))), 
                    "Kalantar_differential_expression_condition.xlsx")


# Plot top genes #### 

head(as.data.frame(results_vaso_shrink_DEG))

#' make_top_genes 
#' function to generate a data frame of top genes for a results table - one per result, can be rbound after
#' @param results_shrunk a deseq results object generated using results function and lfcShrink 
#' @param rep_name a string, name I would like 
#' @param ngenes an integer, to change the number of repetitions of the label - sometimes top_n gives more or less than 10 
#' results because some genes are the same padj 
#' 
#' @return res_top10 a data frame with top10 up and down regulated genes by padj, and their log2 fold changes 
#' for a top10 genes plot 
make_top_genes <- function(results_shrunk, rep_name, ngenes = 10){
  res_top10_up <- as.data.frame(results_shrunk) %>% filter(log2FoldChange > 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
    tibble::add_column(condition = (rep(rep_name,ngenes)), Direction = rep("Up", ngenes)) 
  # AT 0 WPI there is only 1 gene up, so change the new vectors to 1  
  res_top10_down <- as.data.frame(results_shrunk) %>% filter(log2FoldChange < 0) %>% top_n(-10, padj) %>% dplyr::select(-baseMean, -lfcSE, -pvalue) %>% 
    tibble::add_column(condition = (rep(rep_name,ngenes)), Direction = rep("Down", ngenes))
  
  res_top10 <- rbind(res_top10_up, res_top10_down)
  res_top10 <- res_top10 %>%
    mutate(gene_symbol2 = coalesce(gene_symbol, rownames(.))) # this could fail - as.data.frame will remove rownames? 
  res_top10 
  res_top10$log2FoldChange <- as.numeric(as.character(res_top10$log2FoldChange))
  
  return(res_top10)
}

res_top10_sepsis_nonbsi_cond <- make_top_genes(results_vaso_sepsis_nonbsi_cond_shrink, rep_name = "sepsis_non_bsi_yes_vs_no")
res_top10_no_sepsis_cond <- make_top_genes(results_vaso_no_sepsis_cond_shrink, rep_name = "no_sepsis_yes_vs_no")

res_top10_condition <- rbind(res_top10_sepsis_nonbsi_cond, res_top10_no_sepsis_cond)

# Plot top genes 

ggplot(data = res_top10_condition) + 
  geom_point(mapping = aes(x = condition, y = log2FoldChange, color = Direction, group = Direction)) + 
  scale_color_brewer(palette = "Set1", name="Direction", breaks = c("Up", "Down"), labels=c("Upregulated", "Downregulated"), direction = -1) + 
  xlab("Condition") + ylab("Log2 fold change") + 
  geom_hline(mapping = aes(yintercept = 0), linetype = "dotted", size = 1) + 
  ggrepel::geom_text_repel(mapping = aes(x = condition, y = log2FoldChange, group = Direction, label = gene_symbol2), fontface = "italic", max.overlaps = Inf) +
  ggtitle("Top 10 most significant genes - condition")

ggsave("top_10_genes_vaso_group_cond.jpeg") 
# Saving 6.65 x 5.15 in image
# font italics doesn't work, Ubuntu does not have the fonts working, would need to install 


# RESULTS NOTE #### 
# ~~~ Using interaction terms ~~~

dds <- makeExampleDESeqDataSet(n=100,m=18)
dds$genotype <- factor(rep(rep(c("I","II","III"),each=3),2))
design(dds) <- ~ genotype + condition + genotype:condition
dds <- DESeq(dds)
resultsNames(dds)

# the condition effect for genotype I (the main effect)
results(dds, contrast=c("condition","B","A"))

# the condition effect for genotype III.
# this is the main effect *plus* the interaction term
# (the extra condition effect in genotype III compared to genotype I).
results(dds, contrast=list( c("condition_B_vs_A","genotypeIII.conditionB") ))

# the interaction term for condition effect in genotype III vs genotype I.
# this tests if the condition effect is different in III compared to I
results(dds, name="genotypeIII.conditionB")

# the interaction term for condition effect in genotype III vs genotype II.
# this tests if the condition effect is different in III compared to II
results(dds, contrast=list("genotypeIII.conditionB", "genotypeII.conditionB"))