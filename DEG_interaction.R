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

# for on_pressors * Group
# [1] "Intercept"                         "on_pressors_yes_vs_no"             "group_SepsisBSI_vs_No.Sepsis"     
# [4] "group_SepsisNon.BSI_vs_No.Sepsis"  "on_pressorsyes.groupSepsisBSI"     "on_pressorsyes.groupSepsisNon.BSI"

# for Group * on_pressors 
# [1] "Intercept"                         "group_SepsisBSI_vs_No.Sepsis"      "group_SepsisNon.BSI_vs_No.Sepsis" 
# [4] "on_pressors_yes_vs_no"             "groupSepsisBSI.on_pressorsyes"     "groupSepsisNon.BSI.on_pressorsyes"

?results

# The key point to remember about designs with interaction terms is that, unlike for a design ~genotype + condition, 
# where the condition effect represents the overall effect controlling for differences due to genotype, 
# by adding genotype:condition, the main condition effect only represents the effect of condition for the reference 
# level of genotype (I, or whichever level was defined by the user as the reference level). 
# The interaction terms genotypeII.conditionB and genotypeIII.conditionB give the difference between the condition effect 
# for a given genotype and the condition effect for the reference genotype.


results_vaso_group1 <- results(dds_vaso_group, name="on_pressorsyes.groupSepsisBSI")
sum(results_vaso_group1$padj < 0.05, na.rm=TRUE)

# This should be the correct option 
results_vaso_group2 <- results(dds_vaso_group, name="groupSepsisBSI.on_pressorsyes")
sum(results_vaso_group2$padj < 0.05, na.rm=TRUE) # 0 
results_vaso_group3 <- results(dds_vaso_group, name="groupSepsisNon.BSI.on_pressorsyes")
sum(results_vaso_group3$padj < 0.05, na.rm=TRUE) # 0 


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


# MA plot 
plotMA(results_vaso_group4, ylim=c(-2,2))

plotMA(results_vaso_group_shrink, ylim=c(-2,2))

# Shrinkage 

results_vaso_group_shrink <- lfcShrink(dds_vaso_group, coef="on_pressors_yes_vs_no", type="apeglm")
# will use these 

results_vaso_group_shrink

# Annotate and export results 

rownames(results_vaso_group_shrink) <- gsub("\\.\\d{2,2}$", "", rownames(results_vaso_group_shrink))
rownames(results_vaso_group_shrink) <- gsub("\\.\\d{1,1}$", "", rownames(results_vaso_group_shrink))
rownames(results_vaso_group_shrink)

results_vaso_group_shrink$gene_symbol <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = rownames(results_vaso_group_shrink),
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

results_vaso_group_shrink$gene_name <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = rownames(results_vaso_group_shrink),
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "GENENAME", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

write.csv(results_vaso_group_shrink, "DEG_group_vasopressors.csv")

# Check DE gene numbers 

sum(results_vaso_group$padj < 0.05, na.rm=TRUE) # 2410 
results_vaso_group_shrink_DEG <- subset(results_vaso_group_shrink, padj < 0.05)
results_vaso_group_shrink_DEG

write.csv(results_vaso_group_shrink_DEG, "DEG_sig_group_vasopressors.csv")

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

# TO DO #### 
# Split groups and run vasopressors on separate objects #### 



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