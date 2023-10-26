# Library calls #### 
library(DESeq2)
library(readxl)
library("biomaRt")

library(tidyverse)
library(ggpubr)
library(rstatix)
# None of my custom functions used in this file 

# Data #### 

head(PCAdata)
# I don't need the genes, only metadata so this will be enough 
PCAdata %>% group_by(group) %>% dplyr::count()

# Pairwise t test - PC1, vasopressors #### 
pc1_vaso_ttest <- PCAdata %>% 
  t_test(PC1 ~ on_pressors) %>%
  add_significance()
pc1_vaso_ttest

# Add p-value and significance levels
pc1_vaso_ttest <- pc1_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  PCAdata, x = "on_pressors", y = "PC1", color = "on_pressors", 
  ylab = "PC1", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc1_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc1_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("stats/boxplot_pressors_ttest.jpeg")


# Anova with t tests #### 
# ** PC1 on Group ####  

# Anova 
anova_pc1_group <- PCAdata %>% anova_test(PC1 ~ group)
anova_pc1_group # NS 

# Pairwise comparisons
pca1_group_ttest <- PCAdata %>%
  pairwise_t_test(PC1 ~ group, p.adjust.method = "bonferroni")
pca1_group_ttest
write.csv(pca1_group_ttest %>% dplyr::select(-groups), "stats/pairwise_ttests_group.csv")

# Show adjusted p-values
pca1_group_ttest <- pca1_group_ttest %>% add_xy_position(x = "group")
ggboxplot(PCAdata, x = "group", y = "PC1", color = "group", 
          ylab = "PC1", xlab = "Group", add = "jitter") +
  stat_pvalue_manual(pca1_group_ttest, label = "p.adj", tip.length = 0, step.increase = 0.1) + # hide.ns = TRUE, 
  labs(
    subtitle = get_test_label(anova_pc1_group, detailed = TRUE),
    caption = get_pwc_label(pca1_group_ttest)
  ) + scale_color_brewer(palette = "Set1", direction = -1)
# all NS! 

ggsave("stats/boxplot_group_pariwise_ttests.jpeg")

# ** groups where suspect and indeterm are removed ####
anova_pc1_group2 <- PCAdata %>% filter(group %in% c("No-Sepsis", "SepsisBSI", "SepsisNon-BSI")) %>% 
  anova_test(PC1 ~ group)
anova_pc1_group2 # yes, signif 

# Pairwise comparisons
pca1_group_ttest2 <- PCAdata %>% filter(group %in% c("No-Sepsis", "SepsisBSI", "SepsisNon-BSI")) %>% 
  pairwise_t_test(PC1 ~ group, p.adjust.method = "bonferroni") %>% add_xy_position(x = "group")
pca1_group_ttest2
write.csv(pca1_group_ttest2 %>% dplyr::select(-groups), "stats/pairwise_ttests_group_lessgroups.csv")

# adding xy position does not get filtered groups, so manually change xmin and xmax 
pca1_group_ttest2$xmin <- c(1, 1, 2)
pca1_group_ttest2$xmax <- c(2, 3, 3)
pca1_group_ttest2

# Show adjusted p-values
ggboxplot(PCAdata %>%  filter(group %in% c("No-Sepsis", "SepsisBSI", "SepsisNon-BSI")), 
          x = "group", y = "PC1", color = "group", 
          ylab = "PC1", xlab = "Group", add = "jitter") +
  stat_pvalue_manual(pca1_group_ttest2, label = "p.adj", tip.length = 0, step.increase = 0.1) + # hide.ns = TRUE, 
  labs(
    subtitle = get_test_label(anova_pc1_group2, detailed = TRUE),
    caption = get_pwc_label(pca1_group_ttest2)
  ) + scale_color_brewer(palette = "Set1", direction = -1)
# all NS! 

ggsave("stats/boxplot_group_pariwise_ttests_lessgroups.jpeg")

# ** PC1 on Group and vasopressors #### 
# Is there a difference between the use of vasopressors in sepsis vs non-sepsis 

# Anova 
anova_pc1_group_vaso <- PCAdata %>% anova_test(PC1 ~ group + on_pressors + group:on_pressors)
anova_pc1_group_vaso
# Interaction non significant , and neither are single terms 

# Sliced Tukey tests #### 
# ** Pressors and group test - for a grouped plot ####

# Box plot: comparison against reference
stat.test <- PCAdata %>%
  group_by(on_pressors) %>%
  tukey_hsd(PC1 ~ group) # ref.group = "0.5" , previously t_test()
stat.test

# Box plots
stat.test <- stat.test %>% 
  add_xy_position(x = "on_pressors", dodge = 0.8)

ggboxplot(PCAdata, x = "on_pressors", y = "PC1", color = "group", add = "jitter") + 
  stat_pvalue_manual(stat.test,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE)
ggsave("stats/boxplot_vasopressors_group_sliced_tukey.jpeg") # all NS 

# ** reverse - group and pressors  ####
stat.test <- PCAdata %>%
  group_by(group) %>%
  tukey_hsd(PC1 ~ on_pressors) # ref.group = "0.5"
stat.test

# Box plots
stat.test <- stat.test %>% 
  add_xy_position(x = "group", dodge = 0.8)

ggboxplot(PCAdata, x = "group", y = "PC1", color = "on_pressors", add = "jitter") + 
  stat_pvalue_manual(stat.test,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE) + 
  scale_color_brewer(palette = "Set1", direction = -1)
# no p adj in a comparison of 2 levels 
ggsave("stats/boxplot_vasopressors_group_rev_sliced_tukey.jpeg")

# PC1 on SIRS and vasopressors ####

# ** Anova #### 
anova_pc1_sirs_vaso <- PCAdata %>% anova_test(PC1 ~ sirs_total + on_pressors + sirs_total:on_pressors)
anova_pc1_sirs_vaso # SIRS is significant 
# Interaction 0.252 - non significant 

# ** Sliced - pressors / sirs ####

# Box plot: comparison against reference
stat.test2 <- PCAdata %>%
  group_by(on_pressors) %>%
  tukey_hsd(PC1 ~ sirs_total) # ref.group = "0.5", previously t_test()
stat.test2

# Box plots
stat.test2 <- stat.test2 %>% 
  add_xy_position(x = "on_pressors", dodge = 0.8)

ggboxplot(PCAdata, x = "on_pressors", y = "PC1", color = "sirs_total", add = "jitter") + 
  stat_pvalue_manual(stat.test2,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE)
ggsave("stats/boxplot_pressors_sirs_sliced_tukey.jpeg")


# ** Sliced - sirs / pressors #### 
stat.test2 <- PCAdata %>% filter(sirs_total %in% c(2,3,4)) %>%
  group_by(sirs_total) %>%
  tukey_hsd(PC1 ~ on_pressors) # ref.group = "0.5"
stat.test2
# cannot test, because there is 1 level without yes-no

# Box plots
stat.test2 <- stat.test2 %>% 
  add_xy_position(x = "sirs_total", dodge = 0.8)

ggboxplot(PCAdata, x = "sirs_total", y = "PC1", color = "on_pressors", add = "jitter") + 
  stat_pvalue_manual(stat.test2,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE) + 
  scale_color_brewer(palette = "Set1", direction = -1)
# no p adj in a comparison of 2 levels 
ggsave("stats/boxplot_sirs_pressors_sliced_tukey.jpeg")


# NOT IN USE, Later PCs and vasopressors - placeholder #### 

head(pc_scores_toplot) # this is a PCA table with more PCs 

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




