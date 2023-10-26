
#' Function to add annotated gene symbol and name columns to DESeq results 
#' requires org.Hs.eg.db package
#' 
#' @param results_shrunk a DESeq2 results object, can be after lfcShrink or before 
#' gene names are in rownames and must be ensembl IDs 
#' 
#' @return DESeq results object with additional gene_symbol and gene_name columns 
#' and rownames pulled into a column ensembl_id 
#' Note that this only works for human genes at the moment 
annotate_genes <- function(results_shrunk){
  # rownames(results_shrunk) <- gsub("\\.\\d{2,2}$", "", rownames(results_shrunk))
  # rownames(results_shrunk) <- gsub("\\.\\d{1,1}$", "", rownames(results_shrunk))
  # the above is optional for whole blood and should probably be removed (and fixed for raw counts)
  
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


#' Function to generate a plot of differentially expressed genes 
#' 
#' @param results_object a list of results objects that we should count DEGs in - shrunk or not 
#' @param grouping a string for grouping on the plot x axis (time, or other group that is used)
#' @param padj adjusted p value cut-off 
#' @param lfc log2foldchange cut-off for up- and downregulated genes 
#' 
#' @return gene_counts - a data frame 
count_deg_numbers <- function(results_object = NULL, grouping = NULL,
                              padj = 0.05, lfc = 0){
  
  a <- nrow(subset(results_object, padj < 0.05 & log2FoldChange > lfc))
  b <- nrow(subset(results_object, padj < 0.05 & log2FoldChange < lfc))
  count <- c(a,b)
  count <- as.numeric(as.character(count))
  direction <- c("Up", "Down")
  condition <- rep(grouping, 2)
  
  genecounts <- data.frame(condition, direction, count)
  return(genecounts)
  
}


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
  # fix gene names - not needed in the plasma dataset 
  # top_loadings <- gsub("\\.\\d{2,2}$", "", top_loadings)
  # top_loadings <- gsub("\\.\\d{1,1}$", "", top_loadings)
  # get gene names, save and return 
  gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                      values = top_loadings, mart= mart)
  # with vsd object the top loadings changed!
  write.csv(gene_names, paste0("gene_names_",pc,"_",direction,".csv"))
  return(gene_names)
}


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


#' A function to create a paired t test for different PCs and models 
#' can be used with lapply to iterate over PCc/variables 
#' 
#' @param pc_no an integer - number of a PC of interest 
#' @param data PCA data which includes PCs of interest and metadata for the model (formula2)
#' @param formula2 string of a model formula, can be one factor or an interaction
#' formula string is then placed into a full formula of PCn ~ formula2 
#' 
#' @return rstatix tibble with t test results 
#' You can use PC number vector as input to lapply and then combine multiple ttest results with 
#' ttest_list <- do.call(rbind, ttest_list)
pc_general_vaso_ttest <- function(pc_no, data, formula2 = "on_pressors"){
  model_formula <- formula(paste(paste0("PC",pc_no), '~', formula2))
  t_test(data, formula = model_formula) %>%
    add_significance()
}


#' Plot PCA generated by DESeq2 - colored by chosen discrete variables, no faceting 
#' percent_Var must be predefined 
#' percentVar <- round(100 * attr(PCAdata, "percentVar")) 
#' 
#' @param data PCAdata which contains PC1 and PC2, column names must be in capital letters 
#' and metadata you need to plot by (discrete groupings)
#' @param color a string with the grouping to color by (must be a discrete variable)
#' 
#' @return ggplot object - PCA plot with PCs 1 and 2 colored by a discrete grouping
#' 
plot_colored_PCA_discrete_nofacet <- function(data, color){
  p <- ggplot(data = data) + 
    geom_point(mapping = aes(x = PC1, y = PC2, color = .data[[color]])) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    #   facet_wrap(~ on_pressors, nrow = 1, scales = "free") +
    theme(text = element_text(size = 14, family = "Calibri")) + 
    scale_color_brewer(palette="Set1", direction = -1) #+ 
  #  ggtitle("Are patients on vasopressors?")
  return(p)
}

#' Plot PCA generated by DESeq2 - colored by chosen continuous variables, no faceting 
#' percent_Var must be predefined 
#' percentVar <- round(100 * attr(PCAdata, "percentVar")) 
#' 
#' @param data PCAdata which contains PC1 and PC2, column names must be in capital letters 
#' and metadata you need to plot by (continous groupings, top loadings, gene expression)
#' PCAdata with genes can be made using the make_loadings_pcadata function 
#' @param color a string with the grouping to color by (must be a continuous variable)
#' 
#' @return ggplot object - PCA plot with PCs 1 and 2 colored by a continuous grouping
#' 
plot_colored_PCA_cont_nofacet <- function(data, color){
  p <- ggplot(data = data) + 
    geom_point(mapping = aes(x = PC1, y = PC2, color = .data[[color]])) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    #   facet_wrap(~ on_pressors, nrow = 1, scales = "free") +
    theme(text = element_text(size = 14, family = "Calibri")) + 
    scale_color_viridis_b() #+ 
  #  ggtitle("Are patients on vasopressors?")
  return(p)
}


#' Plot PCA with the ability to select PCs - discrete variables, no faceting 
#' #' percent_Var must be predefined 
#' percent_Var <- round(propve * 100, 2)
#' propve calculated in PCA_prcomp script (my_pca.var <- my_pca$sdev ^ 2 ; 
#' propve <- my_pca.var / sum(my_pca.var)  )
#' 
#' @param data PCAdata which contains the PCs you want (not the DESeq PCA which only had 2 PCs) 
#' and metadata you need to plot by (discrete groupings)
#' @param color a string with the grouping to color by (must be a discrete variable)
#' @param x a string, "PC2" or another PC to be plotted on the x axis 
#' @param y a string, PC to be plotted on the y axis 
#' @param var1 numeric, the number of PC on the x axis 
#' @param var2 numeric, the number of PC on the y axis 
#' 
#' @return ggplot object - PCA plot with PCs colored by a discrete grouping
#' 
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


#' Plot PCA with the ability to select PCs - continuous variable, no facets 
#' percent_Var must be predefined 
#' percent_Var <- round(propve * 100, 2)
#' propve calculated in PCA_prcomp script (my_pca.var <- my_pca$sdev ^ 2 ; 
#' propve <- my_pca.var / sum(my_pca.var)  )
#' 
#' @param data PCAdata which contains the PCs you want (not the DESeq PCA which only had 2 PCs) 
#' and metadata you need to plot by (continous groupings, top loadings, gene expression)
#' PCAdata with genes can be made using the make_loadings_pcadata function 
#' @param color a string with the grouping to color by (must be a continuous variable, gene, wbc etc.)
#' @param x a string, "PC2" or another PC to be plotted on the x axis 
#' @param y a string, PC to be plotted on the y axis 
#' @param var1 numeric, the number of PC on the x axis 
#' @param var2 numeric, the number of PC on the y axis 
#' 
#' @return ggplot object - PCA plot with PCs colored by a discrete grouping
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


#' Subset counts and metadata 
#' 
#' @param counts raw gene counts for DESeq object - rownames are genes, columns are samples 
#' @param metadada coldata for DESeq obect - rows are samples, 
#' sample names (eqivalent to column names of counts) are in column "id2"
#' groups to subset are in column "group" 
#' @param groups a string or character vector of group names to keep 
#'
#' @return a list with newcounts and newmetadata - subsetted values 
#' OR a warning if ncol(newcounts) != nrow(newmetadata)

# Or make a complete subsetted dds? cause otherwise I have to return a list of objects 
# (however this is kind of ok I guess)
# Will need to make checks for 

subset_counts_and_metadata <- function(counts, metadata, groups){
  if(length(groups) == 1){
    # Remove samples with specific group in counts 
    samples <- Kalantar_metadata$id2[Kalantar_metadata$group == groups]
    
    newcounts <- counts %>% 
      dplyr::select(all_of(samples))
    
    # Remove samples in coldata 
    newmetadata <- metadata %>% filter(group == groups) 
  }
  else
  {
    # Remove Sepsis.suspect and sepsis.indeterm groups - I don't need them 
    # Remove samples with specific group in counts 
    samples <- metadata$id2[metadata$group %in% groups]
    
    newcounts <- counts %>% 
      dplyr::select(all_of(samples))
    
    # Remove samples in coldata 
    
    newmetadata <- metadata %>% filter(group %in% groups) 

  }
  if(ncol(newcounts) == nrow(newmetadata)){
    return(list(newcounts = newcounts,
                newmetadata = newmetadata))
  }
    else{
      print("Ncol counts don't equal nrow metadata")
    }
}
