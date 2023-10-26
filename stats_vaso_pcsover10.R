# Later PCs and vasopressors #### 

library(tidyverse)
library(ggpubr)
library(rstatix)
source("functions_PCA-transcriptomics.R")

head(pc_scores_toplot) # this is a PCA table with all PCs 

# From PCA_laterPCS.R

# pc_scores_toplot <- pc_scores %>% # id2 is at the end
#   dplyr::full_join(Kalantar_metadata, by = "id2") 


# Lapply a t test and make a table (do.call, rbind)? 

ttest_list <- lapply(1:50, pc_general_vaso_ttest, 
                     data = pc_scores_toplot, formula2 = "on_pressors")

ttest_list <- do.call(rbind, ttest_list)
ttest_list

write.csv(ttest_list, "stats/ttests_vasopressors_pcs1-50.csv")
# PC 5 and 43 are significant on t test for pressors 

# PC15 is significant 


# ** PC5 #### 
pc5_vaso_ttest <- pc_scores_toplot %>% 
  t_test(PC5 ~ on_pressors) %>%
  add_significance()
pc5_vaso_ttest

# Add p-value and significance levels
pc5_vaso_ttest <- pc5_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "PC5", color = "on_pressors", 
  ylab = "PC5", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc5_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc5_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("stats/boxplot_pressors_pc5.jpeg")

# ** PC43 #### 
pc43_vaso_ttest <- pc_scores_toplot %>% 
  t_test(PC43 ~ on_pressors) %>%
  add_significance()
pc43_vaso_ttest

# Add p-value and significance levels
pc43_vaso_ttest <- pc43_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "PC43", color = "on_pressors", 
  ylab = "PC43", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc43_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc43_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("stats/boxplot_pressors_pc43.jpeg")
