# Library calls #### 
library(dplyr)
library(pheatmap)
library("biomaRt")
library(ggplot2)

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

# Geom_tile - get genes within GO terms #### 

head(specific_pathways)

specific_pathways %>% filter(source == "GO:BP") 

strsplit(specific_pathways$intersection,",")  # that makes it a list of character vectors 

specific_pathways$intersection[1]
unlist(strsplit(specific_pathways$intersection[1],","))

# Need a long format data frame with GO term and genes for each term 

GO_BP_genes <- specific_pathways %>% filter(source == "GO:BP") %>% 
  select(term_id, intersection_size, term_name, intersection)

# GO_BP_genes$new_intersection <- strsplit(GO_BP_genes$intersection,",")
# GO_BP_genes$new_intersection # it is a list of vectors 

colnames(GO_BP_genes)

# GO_BO_genes_long <- GO_BP_genes %>% tidyr::separate_rows(new_intersection, convert = TRUE) 
# GO_BO_genes_long$new_intersection[1:50]

GO_BO_genes_long <- GO_BP_genes %>% mutate(intersection = stringr::str_split(intersection, ",")) %>% 
  unnest(cols = c(intersection)) # stringr::str_split

# Get a gene_symbol column 
library(org.Hs.eg.db)

GO_BO_genes_long$gene_symbol <- mapIds(
  org.Hs.eg.db, # Replace with annotation package for your organism
  keys = GO_BO_genes_long$intersection,
  keytype = "ENSEMBL", # Replace with the type of gene identifiers in your data
  column = "SYMBOL", # The type of gene identifiers you would like to map to
  multiVals = "first"
)

length(unique(GO_BO_genes_long$gene_symbol)) # 903 


# Plot with geom_tile 

plot <- ggplot(GO_BO_genes_long, aes(x = gene_symbol, y = term_name, fill= intersection_size)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  theme(axis.text.x = element_text(angle = 90))+
  scale_y_discrete(limits = rev(levels(factor(GO_BO_genes_long$term_name)))) +
#  scale_fill_discrete(na.value = 'white') +
  scale_color_discrete(na.value = 'black')

plot

ggsave("GO_BP_genes_overlapping.jpeg", width = 13, height = 10, units = "in")

# ** Interactive plot #### 
plotly_geomtile <- plotly::ggplotly(plot)
file_name <- "GO_BP_genes_heatmap"
save_directory <- "plotly_plots"

htmlwidgets::saveWidget(widget = plotly_geomtile, 
                        file = glue::glue("{save_directory}/{file_name}.html"), 
                        selfcontained = FALSE, 
                      #  libdir = glue::glue("{save_directory}/{library_name}"), 
                        title = file_name)
?htmlwidgets::saveWidget

GO_BO_genes_long %>% filter(intersection_size < 100) %>% 
  ggplot(aes(x = gene_symbol, y = term_name, fill= intersection_size)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  scale_y_discrete(limits = rev(levels(factor(GO_BO_genes_long$term_name)))) +
  theme(axis.text.x = element_text(angle = 90))+
  #  scale_fill_discrete(na.value = 'white') +
  scale_color_discrete(na.value = 'black') 

ggsave("GO_BP_smaller_pathways.jpeg", width = 28, height = 10, units = "in")

# smaller pathways but also remove non-common genes - TO DO 
GO_BO_genes_long %>% filter(intersection_size < 50) %>% 
  ggplot(aes(x = gene_symbol, y = term_name, fill= intersection_size)) + 
  geom_tile(colour = "black") + xlab("") + ylab("") +
  scale_y_discrete(limits = rev(levels(factor(GO_BO_genes_long$term_name)))) +
  theme(axis.text.x = element_text(angle = 90))+
  #  scale_fill_discrete(na.value = 'white') +
  scale_color_discrete(na.value = 'black') 

ggsave("GO_BP_smallest_pathways.jpeg", width = 28, height = 10, units = "in")
