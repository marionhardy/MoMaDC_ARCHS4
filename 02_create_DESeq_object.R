
## This assumes you have a count matrix and metadata dataframe

library(limma)
library(DESeq2)
library(tidyverse)
library(ggrepel)
library(RColorBrewer)

## Expression analysis of

counts = readRDS('./data_output/counts_AML_THP1_DC_M0_M1_M2.Rds')
coldata = read.csv("./data_output/metadata_AML_THP1_DC_M0_M1_M2.csv")

# Correct for batch effects using the GSE GEO ID

# mat = assay(vsd)
# mm = model.matrix(~condition, colData(vsd))
# mat =- limma::removeBatchEffect(mat, batch=vsd$batch, design=mm)
# assay(vsd) = mat
# plotPCA(vsd)

# Comparing all conditions
# Create the full model for comparison of samples
# This process will take multiple hours because of the number of samples!
# Should be run on a computer that can be occupied for that long
# Does not have to be powerful as R is not made to parallelize.


# treatment + ID
dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                              design = ~treatment+ID) 

# Generate a linear model

dds$condition = relevel(dds$treatment, "DMSO")
dds$condition = relevel(dds$ID, "Macrophage")
dds = estimateSizeFactors(dds)
nc = counts(dds, normalized=TRUE)
filter = rowSums(nc >= 10) >= 4 # keeping genes with at least 10 counts in at least 4 samples
dds = dds[filter,]
dds = DESeq(dds)
dds = DESeq(dds)

# Save it

saveRDS(dds, "./data_output/dds_treatment_ID.Rds")


# treatment + ID
dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                             design = ~treatment) 

# Generate a linear model

dds$condition = relevel(dds$treatment, "DMSO")
dds = estimateSizeFactors(dds)
nc = counts(dds, normalized=TRUE)
filter = rowSums(nc >= 10) >= 4 # keeping genes with at least 10 counts in at least 4 samples
dds = dds[filter,]
dds = DESeq(dds)
dds = DESeq(dds)

# Save it

saveRDS(dds, "./data_output/dds_treatment.Rds")

# ID
dds = DESeqDataSetFromMatrix(countData = counts, colData = coldata,
                             design = ~ID) 

# Generate a linear model

dds$condition = relevel(dds$ID, "Macrophage")

dds = estimateSizeFactors(dds)
nc = counts(dds, normalized=TRUE)
filter = rowSums(nc >= 10) >= 4 # keeping genes with at least 10 counts in at least 4 samples
dds = dds[filter,]
dds = DESeq(dds)

# Save it

saveRDS(dds, "./data_output/dds_ID.Rds")















# Filter out the last unwanted samples based on the metadata annotations



# Checking size factors and dispersion

sizeFactors(dds) # only takes into account the sequencing depth, looks ok
plotDispEsts(dds) # verifies normalization, graph looks a-ok

# Checking PCA

rld <- vst(dds)

p1 <- plotPCA(rld,intgroup="condition") + 
  geom_text_repel(max.overlaps = 15,
                  box.padding = 0.25,
                  segment.color = 'grey50',
                  fontface = "italic",
                  label = rld$condition)+
  labs(title = 'PCA per condition')

print(p1)  

# Checking sample similarity

sampleDists <- dist(t(assay(dds)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(dds$condition, sep="-")
colnames(sampleDistMatrix) <- paste(dds$condition, sep="-")
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

