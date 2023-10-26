# Library calls #### 
library(DESeq2)
library(dplyr)
library(ggplot2)
library(readxl)
library("biomaRt")
library(ggbiplot) # loads plyr!!!! need to specify dplyr functions possibly 
# Can add ellipses onto ggplot with stat_ellipse()
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
source("functions_PCA-transcriptomics.R")

# Load data #### 

Kalantar_counts <- read.csv("GSE189401_earli_plasma_counts_50k.csv", header = TRUE, row.names = 1)
head(Kalantar_counts)
# contains a gene symbol column, currently "X.1" - remove 
Kalantar_gene_symbols <- Kalantar_counts$X.1
Kalantar_gene_symbols # could also add rownames as names to this vector, but for now no 

# remove gene symbol column 
Kalantar_counts <- Kalantar_counts[,-1]

Kalantar_metadata <- read_excel("Kalantar_supplementary.xlsx", sheet = "Data_17")
head(Kalantar_metadata)

ncol(Kalantar_counts) # 138 after removing gene symbols 
nrow(Kalantar_metadata) # 138 
ncol(Kalantar_counts) == nrow(Kalantar_metadata) # TRUE

Kalantar_metadata$ID
# Need to split EARLI from colnames, or add EARLI to ID 
Kalantar_metadata$ID2 <- paste("EARLI", Kalantar_metadata$ID, sep = "_")

colnames(Kalantar_counts) == Kalantar_metadata$ID2
length(intersect(colnames(Kalantar_counts),Kalantar_metadata$ID2)) # 138 
# so, all are the same, just in different order 

# Order Kalantar metadata 
Kalantar_metadata <- Kalantar_metadata[match(colnames(Kalantar_counts), Kalantar_metadata$ID2),]
colnames(Kalantar_metadata) <- colnames(Kalantar_metadata) %>% janitor::make_clean_names()
colnames(Kalantar_metadata)

# Tidy up column of metadata for better plotting ####

# age, temp_max, apacheIII, wbc_max are continuous 
# rest should be factors 
# for continuous data make separate

Kalantar_metadata$x28d_death <- factor(Kalantar_metadata$x28d_death)
Kalantar_metadata$x28d_death
Kalantar_metadata$x28d_death <- ifelse(Kalantar_metadata$x28d_death == 0, "no", "yes")

Kalantar_metadata$intubated <- factor(Kalantar_metadata$intubated)
Kalantar_metadata$intubated <- ifelse(Kalantar_metadata$intubated == 0, "no", "yes")

Kalantar_metadata$on_pressors <- factor(Kalantar_metadata$on_pressors)
Kalantar_metadata$on_pressors <- ifelse(Kalantar_metadata$on_pressors == 0, "no", "yes")

Kalantar_metadata$immunocompromised <- factor(Kalantar_metadata$immunocompromised)
Kalantar_metadata$immunocompromised <- ifelse(Kalantar_metadata$immunocompromised == 0, "no", "yes")

Kalantar_metadata$sirs_total <- factor(Kalantar_metadata$sirs_total)
Kalantar_metadata$age <- as.numeric(as.character(Kalantar_metadata$age))
Kalantar_metadata$age

# transform all here 
# PCAdata <- PCAdata %>% mutate(apacheiii_log = log10(apacheiii + 1),
#                               wbc_max_log = log10(wbc_max + 1))


# DESeq object #### 

dds <- DESeqDataSetFromMatrix(countData = Kalantar_counts,
                              colData = Kalantar_metadata, 
                              design = ~ group + on_pressors) 

nrow(dds) #19939 

# DESeq PCA #### 

# Normalise data with a vst function 
vsd <- vst(dds, blind=FALSE)

#PCA of samples 
plotPCA(vsd, intgroup=c("group", "on_pressors")) # The simplest PCA 
plotPCA(vsd, intgroup="group")
plotPCA(vsd, intgroup="on_pressors")

?plyr::rename
PCAdata <- plotPCA(vsd, intgroup= colnames(Kalantar_metadata), returnData = TRUE) 
PCAdata <- PCAdata %>% dplyr::select(-group) %>% rename(c("group.1" = "group"))

# PCAdata <- PCAdata %>% select(-group)

# ** formatted PCA plot ####
percentVar <- round(100 * attr(PCAdata, "percentVar")) 

ggplot(PCAdata, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
#coord_fixed() 
  scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA/PCA_group.jpeg") 

ggplot(PCAdata, aes(x = PC1, y = PC2, color = on_pressors)) +
  geom_point(size =3) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  #coord_fixed() 
  scale_color_brewer(palette="Set1", direction = -1) + 
  scale_fill_brewer(palette = "Set1", direction = -1) + 
  theme(text = element_text(size = 14, family = "Calibri")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggsave("PCA/PCA_vasopressors.jpeg") 

ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = PC1, y = PC2, color = group)) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ factor(on_pressors), nrow = 1, scales = "free") +
  theme(text = element_text(size = 14, family = "Calibri")) + 
  scale_color_brewer(palette="Set1", direction = -1) + 
  ggtitle("Are patients on vasopressors?")

ggsave("PCA/PCA_vasopressors_group.jpeg") 

ggplot(data = PCAdata) + 
  geom_point(mapping = aes(x = PC1, y = PC2, color = factor(sirs_total))) + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  facet_wrap(~ factor(on_pressors), nrow = 1, scales = "free") +
  theme(text = element_text(size = 14, family = "Calibri")) + 
  scale_color_brewer(palette="Set1", direction = -1) + 
  ggtitle("Are patients on vasopressors?")

ggsave("PCA_vasopressors_SIRS.jpeg") 

# PCA colored by different metadata #### 

# ** Log transform some metadata #### 

PCAdata <- PCAdata %>% mutate(apacheiii_log = log10(apacheiii + 1),
                              wbc_max_log = log10(wbc_max + 1))

# Get data for plotting #### 
colnames(PCAdata)
interesting_metadata_cont <- c("age", "apacheiii", "temp_max","wbc_max",
                               "apacheiii_log", "wbc_max_log")
interesting_metadata_cont
interesting_metadata_discrete <- c("group", "gender", "race",
                                  "ethnicity", "x28d_death",
                                  "sirs_total", "on_pressors",
                                  "immunocompromised", "antibiotics")
interesting_metadata_discrete

# Plot discrete and continuous lists #### 

plot_list_cont <- lapply(interesting_metadata_cont, FUN = plot_colored_PCA_cont_nofacet, data = PCAdata) 
plot_list_discrete <- lapply(interesting_metadata_discrete, FUN = plot_colored_PCA_discrete_nofacet, data = PCAdata) 

## view the first plot
plot_list_cont[[1]]
plot_list_discrete[[1]]

ggpubr::ggarrange(plotlist = plot_list_discrete)
ggsave("PCA/PCA_multi_color_discrete.jpeg", width = 15, height = 12, units = "in")
# Saving 25.4 x 15.4 in image

ggpubr::ggarrange(plotlist = plot_list_cont)
ggsave("PCA/PCA_multi_color_cont.jpeg", width = 15, height = 8, units = "in")
# Severity on a diagonal of PC1 and PC2
# ApacheIII log on the severity diagonal  

# Vascular permeability #### 

genes_vasc_perm <- c("PECAM1", "ICAM1", "VCAM1", "DPP3", "ANGPT2", "IL6", "ALB", 
                     "SELE", "CRP", "SDC1", "ADM", "VEGF", "VWF", "CD142", "F2", 
                     "SERPINE1", "PROC", "THBD") 

# Biomart for vascular permeability genes 

gene_names_vasc_perm <- getBM(filters= "hgnc_symbol", attributes= c("hgnc_symbol","description"),
                    values = genes_vasc_perm, mart= mart)
gene_names_vasc_perm # missing CD142 and VEGF 
write.csv(gene_names_vasc_perm, "genes_vascular_permeability.csv")

# Create a gene counts matrix with gene symbols (counts_data for function make_loadings_pcadata) #### 

Kalantar_counts_fornames <- Kalantar_counts
Kalantar_counts_fornames$gene_id <- rownames(Kalantar_counts_fornames)

head(Kalantar_counts_fornames)

Kalantar_counts_fornames_gene_symbols <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),
                    values = Kalantar_counts_fornames$gene_id, mart= mart)
head(Kalantar_counts_fornames_gene_symbols)

Kalantar_counts_fornames <- full_join(Kalantar_counts_fornames, Kalantar_counts_fornames_gene_symbols, by = c("gene_id" = "ensembl_gene_id"))

# Subset for genes top PCA genes (vascular permeability using function, in a separate file)

PCAdata_vascperm <- make_loadings_pcadata(counts_data = Kalantar_counts_fornames,
                                          pca_data = PCAdata,
                                          loadings_pc_names = gene_names_vasc_perm)

# Plot vascular permeability genes #### 

# Already log transformed 
plot_list_vascperm <- lapply(gene_names_vasc_perm[["hgnc_symbol"]], 
                             FUN = plot_colored_PCA_cont_nofacet, data = PCAdata_vascperm) 
plot_list_vascperm[[1]]
# ! Binned scales only support continuous data 
# Need to check the types of data the gene counts are! 

ggpubr::ggarrange(plotlist = plot_list_vascperm)
ggsave("PCA/PCA_multi_color_vascperm.jpeg", width = 20, height = 12, units = "in")
# Saving 27.3 x 15.4 in image

# Expression nearly non-existent 
