# Later PCs and vasopressors #### 

library(tidyverse)
library(ggpubr)
library(rstatix)

head(pc_scores_toplot) # this is a PCA table with more PCs 

# From PCA_laterPCS.R

pc_scores_toplot2 <- pc_scores[,c(1:55,223)] %>% # id2 is at the end
  dplyr::full_join(Kalantar_metadata, by = "id2") 
pc_scores_toplot2 <- pc_scores_toplot2 %>% janitor::clean_names()

# Lapply a t test and make a table (do.call, rbind)? 

pc_general_vaso_ttest <- function(pc_no, data, formula2){
  model_formula <- formula(paste(paste0("pc",pc_no), '~', "on_pressors"))
  t_test(data, formula = model_formula) %>%
    add_significance()
}

ttest_list <- lapply(1:50, pc_general_vaso_ttest, data = pc_scores_toplot2, formula2 = "on_pressors")

ttest_list <- do.call(rbind, ttest_list)
ttest_list

write.csv(ttest_list, "ttests_vasopressors_pcs1-50.csv")

# PC15 is significant 


# ** pc15 #### 
pc15_vaso_ttest <- pc_scores_toplot2 %>% 
  t_test(pc15 ~ on_pressors) %>%
  add_significance()
pc15_vaso_ttest

# Add p-value and significance levels
pc15_vaso_ttest <- pc15_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot2, x = "on_pressors", y = "pc15", color = "on_pressors", 
  ylab = "pc15", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc15_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc15_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc15.jpeg")

# Look into genes with strongest negative loadings for PC15 #### 

pc_loadings %>% select(gene,PC15) %>% arrange(desc(PC15))

pc_loadings %>% select(gene,PC15) %>% arrange(PC15)

# TO DO #### 
# Generate some tables, sub the genes and explain using biomart / organismdbi 

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")

# ** PC2 #### 
pc2_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc2 ~ on_pressors) %>%
  add_significance()
pc2_vaso_ttest

# Add p-value and significance levels
pc2_vaso_ttest <- pc2_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc2", color = "on_pressors", 
  ylab = "PC2", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc2_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc2_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc2.jpeg")