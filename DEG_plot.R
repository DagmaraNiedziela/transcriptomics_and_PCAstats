library(dplyr)

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

results_vaso_sepsis_nonbsi_cond_shrink 
results_vaso_no_sepsis_cond_shrink

genecounts1 <- count_deg_numbers(results_vaso_sepsis_nonbsi_cond_shrink, 
                                 grouping = "sepsis_non_bsi_yes_vs_no")


genecounts2 <- count_deg_numbers(results_vaso_no_sepsis_cond_shrink, 
                                 grouping = "no_sepsis_yes_vs_no")

genecounts <- rbind(genecounts1, genecounts2)  
genecounts_up <- dplyr::filter(genecounts, direction == "Up") 
genecounts_down <- subset(genecounts, direction == "Down") 

# plot ####
# Can turn this into a function as well 

ggplot(data = genecounts_up) + 
  geom_col(data = genecounts_up, 
           mapping = aes(x = condition, y = count, fill = direction, group = 1), position = "dodge") +
  geom_text(data = genecounts_up, 
            mapping = aes(x = condition, y = count, label = count, group = 1), vjust = -0.6) +
  geom_col(data = genecounts_down, 
           mapping = aes(x = condition, y = -count, fill = direction), position = "dodge") + 
  geom_text(data = genecounts_down, 
            mapping = aes(x = condition, y = -count, label = count), vjust = 1.2) + 
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name="Direction",
                    breaks=c("Up", "Down"),
                    labels=c("Upregulated", "Downregulated")) + 
  scale_y_continuous(expand = c(.1, .1)) +
  # xlab("Weeks post infection (wpi)") + 
  ylab("Number of genes") + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("DEG_numbers_condition.jpeg")

# ggsave(plot, paste0(plot_title, ".jpeg"), dpi = 300)
# return(plot)


# For group #### 

results_group_bsivsnonbsi_shrink 
results_group_bsivsnosep_shrink

genecounts1 <- count_deg_numbers(results_group_bsivsnonbsi_shrink, 
                                 grouping = "bsi_nonbsi")


genecounts2 <- count_deg_numbers(results_group_bsivsnosep_shrink, 
                                 grouping = "bsi_nosep")

genecounts3 <- count_deg_numbers(results_group_nonbsivsnosep_shrink, 
                                 grouping = "nonbsi_nosep")

genecounts <- rbind(genecounts1, genecounts2, genecounts3)  
genecounts_up <- dplyr::filter(genecounts, direction == "Up") 
genecounts_down <- subset(genecounts, direction == "Down") 

# plot ####
# Can turn this into a function as well 

ggplot(data = genecounts_up) + 
  geom_col(data = genecounts_up, 
           mapping = aes(x = condition, y = count, fill = direction, group = 1), position = "dodge") +
  geom_text(data = genecounts_up, 
            mapping = aes(x = condition, y = count, label = count, group = 1), vjust = -0.6) +
  geom_col(data = genecounts_down, 
           mapping = aes(x = condition, y = -count, fill = direction), position = "dodge") + 
  geom_text(data = genecounts_down, 
            mapping = aes(x = condition, y = -count, label = count), vjust = 1.2) + 
  scale_fill_brewer(palette = "Set1", direction = -1,
                    name="Direction",
                    breaks=c("Up", "Down"),
                    labels=c("Upregulated", "Downregulated")) + 
  scale_y_continuous(expand = c(.1, .1)) +
  # xlab("Weeks post infection (wpi)") + 
  ylab("Number of genes") + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("DEG_numbers_group.jpeg")

# ggsave(plot, paste0(plot_title, ".jpeg"), dpi = 300)
# return(plot)

