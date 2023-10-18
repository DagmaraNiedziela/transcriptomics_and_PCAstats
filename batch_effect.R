# Library calls #### 

# BiocManager::install("RUVSeq")
library(RUVSeq)
library(RColorBrewer)

# Data - LFC DEG for sepsis and non sepsis #### 

results_vaso_sep_shrink 

results_vaso_nosep_shrink

# Diagnostic plots #### 
# Boxplot - RUVSeq ####

# ** Sepsis #### 
set <- newSeqExpressionSet(as.matrix(Kalantar_counts_sepsis),
                           phenoData = data.frame(Kalantar_metadata_sepsis, row.names=colnames(Kalantar_counts_sepsis)))
set

library(RColorBrewer)

x <- as.factor(Kalantar_metadata_sepsis$on_pressors)
names(x) <- colnames(Kalantar_counts_sepsis)

colors <- brewer.pal(3, "Set2")

plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], 
        main = "Relative log expr - sepsis")
plotPCA(set, col=colors[x], cex=1.2, 
        main = "PCA - sepsis")

uq <- betweenLaneNormalization(set, which="upper")
plotRLE(uq, outline=FALSE, ylim=c(-4, 4), col=colors[x],
        main = "Relative log expr - sepsis, normalised")
plotPCA(uq, col=colors[x], cex=1.2,
        main = "PCA - sepsis, normalised")


jpeg(file="RUV_boxplots_PCA_sepsis.jpeg", width = 30, height = 20, units = "cm", res = 300)
par(mfrow=c(2,2))
plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x], 
        main = "Relative log expr - sepsis")
plotPCA(set, col=colors[x], cex=1.2, 
        main = "PCA - sepsis")
plotRLE(uq, outline=FALSE, ylim=c(-4, 4), col=colors[x],
        main = "Relative log expr - sepsis, normalised")
plotPCA(uq, col=colors[x], cex=1.2,
        main = "PCA - sepsis, normalised")
dev.off()


# ** no sepsis #### 
set2 <- newSeqExpressionSet(as.matrix(Kalantar_counts_nosepsis),
                            phenoData = data.frame(Kalantar_metadata_nosepsis, 
                                                   row.names=colnames(Kalantar_counts_nosepsis)))
set2

x2 <- as.factor(Kalantar_metadata_nosepsis$on_pressors)
names(x2) <- colnames(Kalantar_counts_nosepsis)

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


jpeg(file="RUV_boxplots_PCA_no_sepsis.jpeg", width = 30, height = 20, units = "cm", res = 300)
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

# Use RUVr - residuals from differential expression 

