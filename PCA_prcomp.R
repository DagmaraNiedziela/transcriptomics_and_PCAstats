# Library calls #### 
library(dplyr)
library(ggplot2)
library("biomaRt")
library(ggbiplot) # loads plyr!!!! need to specify dplyr functions possibly 
# Can add ellipses onto ggplot with stat_ellipse()
source("functions_PCA-transcriptomics.R")

# Regular prcomp with arrows #### 
# Could do this without scaling, using vst normalised data? Then should be the same 

# Apply PCA using prcomp function
# Need to scale / Normalize as PCA depends on distance measure

# Using vst 
# To get it the same as DeSeq2 I need to select 500 most variable genes from the object 
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, 
                                                   length(rv)))]

my_pca <- prcomp(t(assay(vsd)[select,]), scale = FALSE,
                 center = TRUE, retx = TRUE)

names(my_pca)

# Summary
summary(my_pca) # 138 PCs 
my_pca

# View the principal component loading
# my_pca$rotation[1:5, 1:4]
my_pca$rotation

# See the principal components
dim(my_pca$x)
my_pca$x

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
jpeg(file="PCA/scree_plot_prcomp_vsd.jpeg")
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
# which will at least cover 90 % variance of dimension
which(cumsum(propve) >= 0.9)[1] # 80 

# Ellipses plot #### 
ellipses_group <- ggbiplot(my_pca, 
                           obs.scale = 1,
                           var.scale = 1, 
                           groups = Kalantar_metadata$group, # would need to be ordered the same? same for the other pca? 
                           ellipse = TRUE, 
                           circle = FALSE, 
                           ellipse.prob = 0.68, 
                           var.axes = FALSE) + 
  scale_color_discrete(name = '') + 
  theme(legend.direction = 'horizontal', legend.position = 'top')

ellipses_group
ggsave("PCA/prcomp_ellipses_vsd_group.jpeg")

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

ggsave("PCA/prcomp_ellipses_vasopressors_vsd.jpeg")

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

ggsave("PCA/scree_ggplot_vsd.jpeg")

# ** PC scores #### 

pc_scores <- my_pca$x 

pc_scores <- pc_scores %>% 
  as_tibble(rownames = "id2")

# pc_scores
# max(pc_scores$PC1)
# 
# pc_plot <- pc_scores %>% 
#   full_join(Kalantar_metadata, by = "id2") %>% 
#   ggplot(aes(x = PC1, y = PC2, color = group)) + # NA
#   geom_point()
# 
# pc_plot


# ** Loadings #### 

pc_loadings <- my_pca$rotation

pc_loadings <- pc_loadings %>% 
  as_tibble(rownames = "gene")

pc_loadings

# # This would work if I only wanted 1 PC
# top_genes <- get_loadings_gene_list(pca_data = pc_loadings,
#                                     pc = "PC1", direction = "positive", number = 10)

# # top genes with highest loadings on PC1 and PC2 
top_genes <- pc_loadings %>%
  dplyr::select(gene, PC1, PC2) %>%
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

# ** Biomart Explain top genes ####
top_gene_names <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol","description"),
                        values = top_genes, mart= mart)
top_gene_names

write.csv(top_gene_names, "top_loadings_genes_vsd.csv")

# # join gene names 
top_loadings <- cbind(top_loadings, top_gene_names)
top_loadings

# plot loadings 
loadings_plot <- ggplot(data = top_loadings) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
               arrow = arrow(length = unit(0.1, "in")),
               color = "brown") +
  ggrepel::geom_text_repel(aes(x = PC1, y = PC2, label = hgnc_symbol),
                           nudge_x = -0.005, nudge_y = 0, size = 3) +
  scale_x_continuous(expand = c(0.02, 0.02))

loadings_plot

# ** Arrange PCA and loadings together ####

ggpubr::ggarrange(ellipses_group, loadings_plot, ncol = 2, widths = c(1.5,1))

ggsave("PCA/pca_loadings_vsd.jpeg")
# Saving 10.6 x 4.84 in image

# Plot top loading genes #### 
PCAdata_loadings <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames,
                                          pca_data = PCAdata,
                                          loadings_pc_names = top_gene_names)

# Log transformed 
plot_list_loadings <- lapply(top_gene_names$hgnc_symbol, 
                             FUN = plot_colored_PCA_cont_nofacet, data = PCAdata_loadings) 
plot_list_loadings[[1]]

ggpubr::ggarrange(plotlist = plot_list_loadings)
ggsave("PCA/PCA_multi_color_toploading.jpeg", width = 24, height = 12, units = "in")
# Saving 27.3 x 15.4 in image