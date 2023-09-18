# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)

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
