# Library calls #### 

library(RUVSeq)
library(RColorBrewer)
# detach(RUVSeq)
# My custom functions not used in this file 

# Data - LFC DEG for sepsis and non sepsis #### 

results_vaso_sep_shrink 

results_vaso_nosep_shrink

# Diagnostic plots #### 
# Boxplot - RUVSeq ####

# ** Sepsis #### 
set <- newSeqExpressionSet(as.matrix(Kalantar_sepsisonly$newcounts),
                           phenoData = data.frame(Kalantar_sepsisonly$newmetadata, row.names=colnames(Kalantar_sepsisonly$newcounts)))
set

x <- as.factor(Kalantar_sepsisonly$newmetadata$on_pressors)
names(x) <- colnames(Kalantar_sepsisonly$newcounts)

colors <- brewer.pal(3, "Set2")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], 
        main = "Relative log expr - sepsis")
# No counts!!!! 
plotPCA(set, col=colors[x], cex=1.2, 
        main = "PCA - sepsis")

uq <- betweenLaneNormalization(set, which="upper")
plotRLE(uq, outline=FALSE, ylim=c(-4, 4), col=colors[x],
        main = "Relative log expr - sepsis, normalised")
warnings()
plotPCA(uq, col=colors[x], cex=1.2,
        main = "PCA - sepsis, normalised")

dir.create("batch_effect")

jpeg(file="batch_effect/RUV_boxplots_PCA_sepsis.jpeg", 
     width = 15, height = 10, units = "cm", res = 300) # width = 30, height = 20, units = "cm", 
# par(mfrow=c(2,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], 
        main = "Relative log expr - sepsis")
# plotPCA(set, col=colors[x], cex=1.2, 
#         main = "PCA - sepsis")
# plotRLE(uq, outline=FALSE, ylim=c(-4, 4), col=colors[x],
#         main = "Relative log expr - sepsis, normalised")
# plotPCA(uq, col=colors[x], cex=1.2,
#         main = "PCA - sepsis, normalised")
dev.off()

# ** no sepsis #### 
set2 <- newSeqExpressionSet(as.matrix(Kalantar_nosepsis$newcounts),
                            phenoData = data.frame(Kalantar_nosepsis$newmetadata, 
                                                   row.names=colnames(Kalantar_nosepsis$newcounts)))
set2

x2 <- as.factor(Kalantar_nosepsis$newmetadata$on_pressors)
names(x2) <- colnames(Kalantar_nosepsis$newcounts)

# colors <- brewer.pal(3, "Set2")

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x2], 
        main = "Relative log expr - sepsis")
plotPCA(set2, col=colors[x2], cex=1.2, 
        main = "PCA - sepsis")

uq2 <- betweenLaneNormalization(set2, which="upper")
plotRLE(uq2, outline=FALSE, ylim=c(-4, 4), col=colors[x2],
        main = "Relative log expr - sepsis, normalised")
plotPCA(uq2, col=colors[x2], cex=1.2,
        main = "PCA - sepsis, normalised")


jpeg(file="batch_effect/RUV_boxplots_PCA_no_sepsis.jpeg", width = 30, height = 20, units = "cm", res = 300)
par(mfrow=c(2,2))
plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x2], 
        main = "Relative log expr - no sepsis")
plotPCA(set2, col=colors[x2], cex=1.2, 
        main = "PCA - no sepsis")
plotRLE(uq2, outline=FALSE, ylim=c(-4, 4), col=colors[x2],
        main = "Relative log expr - no sepsis, normalised")
plotPCA(uq2, col=colors[x2], cex=1.2,
        main = "PCA - no sepsis, normalised")
dev.off()

# Plot p values - unadjusted ####
# https://github.com/drisso/peixoto2015_tutorial/blob/master/Peixoto_Additional_File_1.Rmd 
# tutorial from the normalisation/batch effect publication 

jpeg(file="pvalues_histogram_sepsis.jpeg")
hist(results_vaso_sep_shrink$pvalue, main="Sepsis unadjusted p values", 
     xlab="p-value", breaks=100) #  ylim=c(0, 900)
dev.off()

jpeg(file="pvalues_histogram_nosepsis.jpeg")
hist(results_vaso_nosep_shrink$pvalue, main="No sepsis unadjusted p values", 
     xlab="p-value", breaks=100, ylim=c(0, 1500)) # 
dev.off()

jpeg(file="pvalues_histogram_2panels.jpeg", width = 20, height = 10, units = "cm", res = 300)
par(mfrow=c(1,2))
hist(results_vaso_sep_shrink$pvalue, main="Sepsis unadjusted p values", 
     xlab="p-value", breaks=100) #  ylim=c(0, 900)
hist(results_vaso_nosep_shrink$pvalue, main="No sepsis unadjusted p values", 
     xlab="p-value", breaks=100, ylim=c(0, 1500))
dev.off()

jpeg(file="pvalues_histogram_both.jpeg")
hist(results_vaso_shrink$pvalue, main="All groups unadjusted p values", 
     xlab="p-value", breaks=100) # ylim=c(0, 1600)
dev.off()


# Batch correction for no sepsis #### 

#  ** Use RUVr - residuals from differential expression - no sepsis ####

genes <- rownames(Kalantar_counts_nosepsis) # all genes

design <- model.matrix(~x2, data=pData(set2))
y <- DGEList(counts=counts(set2), group=x2)
y <- calcNormFactors(y, method="upperquartile") # Uq normalisation is already done here 

y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
res <- residuals(fit, type="deviance")

set_batch_nosepsis <- RUVr(set2, genes, k=1, res) # ran with k = 1, 2, 3, 4, 5 and 8
pData(set_batch_nosepsis)
str(pData(set_batch_nosepsis))
pData(set_batch_nosepsis)$W_1 # a "batch" column is a set of numbers, different for each sample, what do you do with that? 

jpeg(file="RUV_boxplots_PCA_no_sepsis_RUVr.jpeg", width = 30, height = 20, units = "cm", res = 300)
par(mfrow=c(2,2))

plotRLE(set2, outline=FALSE, ylim=c(-4, 4), col=colors[x2], 
        main = "Relative log expr - no sepsis")
plotPCA(set2, col=colors[x2], cex=1.2, 
        main = "PCA - no sepsis")
# These are new: 
plotRLE(set_batch_nosepsis, outline=FALSE, ylim=c(-4, 4), col=colors[x2], 
        main = "Relative log expr - no sepsis, RUVr")
plotPCA(set_batch_nosepsis, col=colors[x2], cex=1.2, 
        main = "PCA - no sepsis, RUVr")

dev.off()


# New differential expression with batch effect - not done 
# Would need to include batch in DESeq2, need to add model.matrix to DESeq2? 

design <- model.matrix(~x + W_1, data=pData(set1))
y <- DGEList(counts=counts(set), group=x)
y <- calcNormFactors(y, method="upperquartile")
y <- estimateGLMCommonDisp(y, design)
y <- estimateGLMTagwiseDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef=2)
topTags(lrt)