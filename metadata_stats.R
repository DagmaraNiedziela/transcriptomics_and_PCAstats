# Library calls #### 
library(DESeq2)
library(readxl)
library("biomaRt")
library(ggbiplot) # loads plyr!!!! need to specify dplyr functions possibly 

library(tidyverse)
library(ggpubr)
library(rstatix)

library(lsmeans)
library(multcomp)
library(multcompView)

# Data #### 

head(PCAdata)
# I don't need the genes, only metadata so this will be enough 


# Pairwise t test - PC1, vasopressors #### 
pc1_vaso_ttest <- PCAdata %>% 
  t_test(pc1 ~ on_pressors) %>%
  add_significance()
pc1_vaso_ttest

# Add p-value and significance levels
pc1_vaso_ttest <- pc1_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  PCAdata, x = "on_pressors", y = "pc1", color = "on_pressors", 
  ylab = "PC1", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc1_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc1_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors.jpeg")

# Anova with t tests #### 
# ** PC1 on Group ####  

# Anova 
anova_pc1_group <- PCAdata %>% anova_test(pc1 ~ group)
anova_pc1_group

# Pairwise comparisons
pca1_group_ttest <- PCAdata %>%
  pairwise_t_test(pc1 ~ group, p.adjust.method = "bonferroni")
pca1_group_ttest

# Show adjusted p-values
pca1_group_ttest <- pca1_group_ttest %>% add_xy_position(x = "group")
ggboxplot(PCAdata, x = "group", y = "pc1", color = "group", 
          ylab = "PC1", xlab = "Group", add = "jitter") +
  stat_pvalue_manual(pca1_group_ttest, hide.ns = TRUE, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(
    subtitle = get_test_label(anova_pc1_group, detailed = TRUE),
    caption = get_pwc_label(pca1_group_ttest)
  ) + scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_group.jpeg")

# ** PC1 on Group and vasopressors #### 
# Is there a difference between the use of vasopressors in sepsis vs non-sepsis 

# Anova 
anova_pc1_group_vaso <- PCAdata %>% anova_test(pc1 ~ group + on_pressors + group:on_pressors)
anova_pc1_group_vaso
# Interaction 0.4 - non significant 

# Create an interaction column 
PCAdata <- PCAdata %>% 
  mutate(group_vaso = paste(group, on_pressors, sep = "_"))

# Pairwise comparisons
pca1_group_vaso_ttest <- PCAdata %>%
  pairwise_t_test(pc1 ~ group_vaso, p.adjust.method = "bonferroni")
pca1_group_vaso_ttest
write.csv(pca1_group_vaso_ttest %>% select(-groups), "pca1_group_vaso_ttest.csv")

# Show adjusted p-values
pca1_group_vaso_ttest <- pca1_group_vaso_ttest %>% add_xy_position(x = "group_vaso")
ggboxplot(PCAdata, x = "group_vaso", y = "pc1", color = "group",
          ylab = "PC1", xlab = "Group and vasopressor use", add = "jitter") +
  stat_pvalue_manual(pca1_group_vaso_ttest, hide.ns = TRUE, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(
    subtitle = get_test_label(anova_pc1_group_vaso, detailed = TRUE),
    caption = get_pwc_label(pca1_group_vaso_ttest)
  ) + scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_group_vaso.jpeg")

# *** Tukey comparisons #### 

tukey_vaso <- PCAdata %>% tukey_hsd(pc1 ~ group + on_pressors + group:on_pressors, which = "group")
tukey_vaso

# Slicing 
fittukey_groupsvaso <- lsmeans(lm(pc1 ~ group + on_pressors + group:on_pressors, data = PCAdata), 
                               pairwise ~ group | on_pressors)[[2]]
fittukey_groupsvaso

fittukey_groupsvaso2 <- lsmeans(lm(pc1 ~ group + on_pressors + group:on_pressors, data = PCAdata), 
                               pairwise ~ on_pressors | group)[[2]]
fittukey_groupsvaso2

# tukey_swim_temp <- as.data.frame(cld(fit.tukey_swim))
# tukey_swim_temp
# write.csv(tukey_swim_temp, "Tukey_swim_sliced_by_temp.csv") 

# Both lsmeans and emmeans work, but the wrapper does not 
# fittukey_groupsvaso3 <- emmeans_test(formula = pc1 ~ group * on_pressors, 
#                                      data = PCAdata)
# pairwise_emmeans_test()

# Ls means will be deprecated, this is emmeans package equivalent 
EMM <- emmeans(lm(pc1 ~ group + on_pressors + group:on_pressors, data = PCAdata), 
               ~ group * on_pressors)
EMM    # display the cell means

### Simple pairwise comparisons...
pairs(EMM, simple = "group")    # compare treats for each dose -- "simple effects"
pairs(EMM, simple = "on_pressors")     # compare doses for each treat
# comparisons = pairwise ~ group | on_pressors

# *** Group a test - for a grouped plot ####

# Box plot: comparison against reference
stat.test <- PCAdata %>%
  group_by(on_pressors) %>%
  tukey_hsd(pc1 ~ group) # ref.group = "0.5" , previously t_test()
stat.test

# Box plots
stat.test <- stat.test %>% 
  add_xy_position(x = "on_pressors", dodge = 0.8)

ggboxplot(PCAdata, x = "on_pressors", y = "pc1", color = "group") + 
  stat_pvalue_manual(stat.test,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE)
ggsave("boxplot_vasopressors_group_sliced_tukey.jpeg")

# reverse 
stat.test <- PCAdata %>%
  group_by(group) %>%
  tukey_hsd(pc1 ~ on_pressors) # ref.group = "0.5"
stat.test

# Box plots
stat.test <- stat.test %>% 
  add_xy_position(x = "group", dodge = 0.8)

ggboxplot(PCAdata, x = "group", y = "pc1", color = "on_pressors") + 
  stat_pvalue_manual(stat.test,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE) + 
  scale_color_brewer(palette = "Set1", direction = -1)
# no p adj in a comparison of 2 levels 
ggsave("boxplot_vasopressors_group_rev_sliced_tukey.jpeg")

# ** PC1 on SIRS and vasopressors ####

# Anova 
anova_pc1_sirs_vaso <- PCAdata %>% anova_test(pc1 ~ sirs_total + on_pressors + sirs_total:on_pressors)
anova_pc1_sirs_vaso
# Interaction 0.5 - non significant 

# Create an interaction column 
PCAdata <- PCAdata %>% 
  mutate(sirs_vaso = paste(sirs_total, on_pressors, sep = "_"))

# Pairwise comparisons
pca1_sirs_vaso_ttest <- PCAdata %>%
  pairwise_t_test(pc1 ~ sirs_vaso, p.adjust.method = "bonferroni")
pca1_sirs_vaso_ttest
write.csv(pca1_sirs_vaso_ttest, "pca1_sirs_vaso_ttest.csv")

# Show adjusted p-values
pca1_sirs_vaso_ttest <- pca1_sirs_vaso_ttest %>% add_xy_position(x = "sirs_vaso")
ggboxplot(PCAdata, x = "sirs_vaso", y = "pc1", color = "sirs_total",
          ylab = "PC1", xlab = "SIRS and vasopressor use", add = "jitter") +
  stat_pvalue_manual(pca1_sirs_vaso_ttest, hide.ns = TRUE, label = "p.adj", tip.length = 0, step.increase = 0.1) +
  labs(
    subtitle = get_test_label(anova_pc1_sirs_vaso, detailed = TRUE),
    caption = get_pwc_label(pca1_sirs_vaso_ttest)
  ) + scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_sirs_vaso.jpeg")

# *** Tukey comparisons #### 

# Slicing 
fittukey_sirsvaso <- lsmeans(lm(pc1 ~ sirs_total + on_pressors + on_pressors:sirs_total, data = PCAdata), 
                               pairwise ~ sirs_total | on_pressors)[[2]]
fittukey_sirsvaso
??cld
cld(fittukey_sirsvaso) 

fittukey_sirsvaso2 <- lsmeans(lm(pc1 ~ sirs_total + on_pressors + sirs_total:on_pressors, data = PCAdata), 
                                pairwise ~ on_pressors | sirs_total)[[2]]
fittukey_sirsvaso2

# *** Group a test - for a grouped plot ####

# Box plot: comparison against reference
stat.test2 <- PCAdata %>%
  group_by(on_pressors) %>%
  tukey_hsd(pc1 ~ sirs_total) # ref.group = "0.5", previously t_test()
stat.test2

# Box plots
stat.test2 <- stat.test2 %>% 
  add_xy_position(x = "on_pressors", dodge = 0.8)

ggboxplot(PCAdata, x = "on_pressors", y = "pc1", color = "sirs_total") + 
  stat_pvalue_manual(stat.test2,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE)
ggsave("Vasopressors_sirs_sliced_tukey.jpeg")


# reverse 
stat.test2 <- PCAdata %>% filter(sirs_total %in% c(2,3,4)) %>%
  group_by(sirs_total) %>%
  tukey_hsd(pc1 ~ on_pressors) # ref.group = "0.5"
stat.test2
# cannot test, because there is 1 level without yes-no

# Box plots
stat.test2 <- stat.test2 %>% 
  add_xy_position(x = "sirs_total", dodge = 0.8)

ggboxplot(PCAdata, x = "sirs_total", y = "pc1", color = "on_pressors") + 
  stat_pvalue_manual(stat.test2,   label = "p.adj", tip.length = 0.01, hide.ns = TRUE) + 
  scale_color_brewer(palette = "Set1", direction = -1)
# no p adj in a comparison of 2 levels 
ggsave("boxplot_vasopressors_sirs_rev_sliced_tukey.jpeg")


# Later PCs and vasopressors #### 

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

# ** PC3 #### 
pc3_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc3 ~ on_pressors) %>%
  add_significance()
pc3_vaso_ttest

# Add p-value and significance levels
pc3_vaso_ttest <- pc3_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc3", color = "on_pressors", 
  ylab = "PC3", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc3_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc3_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc3.jpeg")

# ** PC4 #### 
pc4_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc4 ~ on_pressors) %>%
  add_significance()
pc4_vaso_ttest

# Add p-value and significance levels
pc4_vaso_ttest <- pc4_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc4", color = "on_pressors", 
  ylab = "PC4", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc4_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc4_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc4.jpeg")

# ** PC5 #### 
pc5_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc5 ~ on_pressors) %>%
  add_significance()
pc5_vaso_ttest

# Add p-value and significance levels
pc5_vaso_ttest <- pc5_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc5", color = "on_pressors", 
  ylab = "PC5", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc5_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc5_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc5.jpeg")

# ** pc6 #### 
pc6_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc6 ~ on_pressors) %>%
  add_significance()
pc6_vaso_ttest

# Add p-value and significance levels
pc6_vaso_ttest <- pc6_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc6", color = "on_pressors", 
  ylab = "pc6", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc6_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc6_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc6.jpeg")

# ** pc7 #### 
pc7_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc7 ~ on_pressors) %>%
  add_significance()
pc7_vaso_ttest

# Add p-value and significance levels
pc7_vaso_ttest <- pc7_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc7", color = "on_pressors", 
  ylab = "pc7", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc7_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc7_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc7.jpeg")

# ** pc8 #### 
pc8_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc8 ~ on_pressors) %>%
  add_significance()
pc8_vaso_ttest

# Add p-value and significance levels
pc8_vaso_ttest <- pc8_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc8", color = "on_pressors", 
  ylab = "pc8", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc8_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc8_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc8.jpeg")

# ** pc9 #### 
pc9_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc9 ~ on_pressors) %>%
  add_significance()
pc9_vaso_ttest

# Add p-value and significance levels
pc9_vaso_ttest <- pc9_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc9", color = "on_pressors", 
  ylab = "pc9", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc9_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc9_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc9.jpeg")

# ** pc10 #### 
pc10_vaso_ttest <- pc_scores_toplot %>% 
  t_test(pc10 ~ on_pressors) %>%
  add_significance()
pc10_vaso_ttest

# Add p-value and significance levels
pc10_vaso_ttest <- pc10_vaso_ttest %>% add_xy_position(x = "on_pressors")

# Create a box-plot
ggboxplot(
  pc_scores_toplot, x = "on_pressors", y = "pc10", color = "on_pressors", 
  ylab = "pc10", xlab = "Vasopressor use", add = "jitter"
) + 
  stat_pvalue_manual(pc10_vaso_ttest, tip.length = 0) +
  labs(subtitle = get_test_label(pc10_vaso_ttest, detailed = TRUE)) + 
  scale_color_brewer(palette = "Set1", direction = -1)

ggsave("boxplot_pressors_pc10.jpeg")


# Categorical variable tests - contignency tables #### 

# Number of people on vasopressors in each group 


# Number of people with specific SIRS scores in each Group 



