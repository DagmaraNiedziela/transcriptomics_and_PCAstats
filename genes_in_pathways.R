# Library calls #### 
library(dplyr)
library(pheatmap)
library("biomaRt")

# Prepare data #### 

results_vaso_shrink_fh <- as.data.frame(results_vaso_shrink) %>% tibble::rownames_to_column("gene_id")

results_vaso_sep_shrink_fh <- as.data.frame(results_vaso_sep_shrink) %>% tibble::rownames_to_column("gene_id")

# A heatmap of genes within pathways #### 

marker_genes_onegene <- full_join(results_vaso_shrink_fh, results_vaso_sep_shrink_fh, by = "gene_id", suffix = c("_vaso", "_vaso_sep"))
marker_genes_onegene
colnames(marker_genes_onegene)

# Make not significant genes NA? 
marker_genes_onegene[c("padj_vaso", "padj_vaso_sep")][(marker_genes_onegene[c("padj_vaso", "padj_vaso_sep")]) > 0.05] <- NA
marker_genes_onegene$log2FoldChange_vaso <- ifelse(is.na(marker_genes_onegene$padj_vaso), NA, marker_genes_onegene$log2FoldChange_vaso)
marker_genes_onegene$log2FoldChange_vaso_sep <- ifelse(is.na(marker_genes_onegene$padj_vaso_sep), 
                                                       NA, marker_genes_onegene$log2FoldChange_vaso_sep)

# ** Neutrophil degranulation #### 

# Get gene symbols - intersection - ENSEMBL ids 
neutrophil_genes <- unique(c(unlist(strsplit(pathways_vaso$result$intersection[pathways_vaso$result$term_id == "REAC:R-HSA-6798695"], ",")),
                            unlist(strsplit(pathways_vaso_sep$result$intersection[pathways_vaso_sep$result$term_id == "REAC:R-HSA-6798695"], ","))))
                  

# Get gene symbols 
neutrophil_genes
netrophil_gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                    values = neutrophil_genes, mart= mart)

netrophil_gene_names
neutrophil_genes_name <- sort(netrophil_gene_names$hgnc_symbol)
length(neutrophil_genes_name)

neutrophil_LFC <- marker_genes_onegene %>% filter(gene_id %in% neutrophil_genes) %>% 
  arrange(gene_symbol_vaso) %>% select(gene_symbol_vaso, log2FoldChange_vaso, log2FoldChange_vaso_sep)

# color palette 
paletteLength <- 50
myColor <- colorRampPalette(c("blue", "white", "red2"))(paletteLength)
# length(breaks) == length(paletteLength) + 1
# use floor and ceiling to deal with even/odd length pallettelengths
myBreaks <- c(seq(min(neutrophil_LFC$log2FoldChange_vaso, na.rm = TRUE), 0, length.out=ceiling(paletteLength/2 + 1)), 
              seq((max(neutrophil_LFC$log2FoldChange_vaso, na.rm = TRUE)/paletteLength), max(neutrophil_LFC$log2FoldChange_vaso, na.rm = TRUE), length.out=floor(paletteLength/2)))
# make breaks the same for all heatmaps to compare expression

# get rownames 
neutrophil_LFC2 <- neutrophil_LFC
head(neutrophil_LFC2)
rownames(neutrophil_LFC2) <- neutrophil_LFC2$gene_symbol_vaso
neutrophil_LFC2 <- neutrophil_LFC2[,-1]

# heatmap 
neutrophil_heatmap <- pheatmap(neutrophil_LFC2, na_col = "grey", 
                              cluster_cols = FALSE, cluster_rows = FALSE, 
                              cellwidth = 15, fontsize = 7, angle_col = 315,
                              color=myColor, breaks = myBreaks) 


tiff("neutrophil_heatmap.tiff", height = 38, width = 20, units = 'cm', 
     compression = "lzw", res = 600) 
neutrophil_heatmap
dev.off() 

Cairo(file="neutrophil_heatmap.png", 
      type="png",
      units="cm", 
      width=24, 
      height=20, 
      pointsize=60, 
      dpi=600)
neutrophil_heatmap
dev.off() 