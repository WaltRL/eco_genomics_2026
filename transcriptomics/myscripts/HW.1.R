### Ecological Genomics - Homework #1 - 10/13/2026

### WGNCA

setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")
options(stringsAsFactors = FALSE);
library(WGCNA);
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot) 
library(gridExtra)
library(Rmisc) 

###Setting Up and Running DESeq
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 

# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# All of our samples look good; the genes that don't look good shouldn't cluster well and 
# should end up in the "grey", catch all module for uncorrelated genes.

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order

# perform variance stabilization
dds_norm <- vst(dds)
#dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices 

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power_6 <- 6 # should try running with 7 based on R^2
soft_power_8 <- 8
soft_power_4 <- 4
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package

###WGNCA Visualizations

bwnet_4 <- readRDS("bwnet.4.rds")
bwnet_6 <- readRDS("bwnet.6.rds")
bwnet_8 <- readRDS("bwnet.8.rds")

###Bwent with threshold of 4 ---------------------------------------------------------
module_eigengenes_4 <- bwnet_4$MEs

head(module_eigengenes_4)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_4$dendrograms[[1]], cbind(bwnet_4$unmergedColors, bwnet_4$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

##All Together
nets <- list("Power 4" = bwnet_4,
             "Power 6" = bwnet_6,
             "Power 8" = bwnet_8)

# Function: dendrogram + unmerged/merged module colors
plot_dendro <- function(net, title, block = 1) {
  # genes belonging to this block (handles multi-block runs safely)
  idx <- net$blockGenes[[block]]
  
  plotDendroAndColors(
    net$dendrograms[[block]],
    cbind(net$unmergedColors[idx], net$colors[idx]),
    groupLabels = c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Gene Dendrogram and Module Colors -", title)
  )
}

###Bwent with threshold of 6 ---------------------------------------------------------
module_eigengenes_6 <- bwnet_6$MEs

head(module_eigengenes_6)

# get number of genes for each module
table(bwnet_6$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_6$dendrograms[[1]], cbind(bwnet_6$unmergedColors, bwnet_6$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

###Bwent with threshold of 8 ---------------------------------------------------------
module_eigengenes_8 <- bwnet_8$MEs

head(module_eigengenes_8)

# get number of genes for each module
table(bwnet_8$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_8$dendrograms[[1]], cbind(bwnet_8$unmergedColors, bwnet_8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#Stitching them together
for (nm in names(nets)) {
  png(paste0("dendrogram_", gsub(" ", "_", nm), ".png"),
      width = 2400, height = 1200, res = 200)
  plot_dendro(nets[[nm]], nm)
  dev.off()
}

install.packages("magick")
library(magick)
setwd("~/projects/eco_genomics_2026/transcriptomics/images")

files <- c("dendrogram_Power_4.png",
           "dendrogram_Power_6.png",
           "dendrogram_Power_8.png")

imgs <- image_read(files)

# Stack vertically (one above the other), which suits wide dendrograms
stacked <- image_append(imgs, stack = TRUE)
image_write(stacked, "dendrograms_stacked.png", format = "png")

# View it
print(stacked)

setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")

###Stats for different thresholds
#Number of modules
length(unique(bwnet_4$colors)) #10
length(unique(bwnet_6$colors)) #11
length(unique(bwnet_8$colors)) #15

#Percentage of genes changing module assignment
mean(bwnet_4$colors != bwnet_6$colors) #0.1945368
mean(bwnet_4$colors != bwnet_8$colors) #0.2805226
mean(bwnet_6$colors != bwnet_8$colors) #0.3032858

#Module eigengene similarity
cor(bwnet_4$MEs, bwnet_6$MEs)
cor(bwnet_4$MEs, bwnet_8$MEs)

#Graph comparing eigengenes
### Ecological Genomics - Homework #1 - 10/13/2026

### WGNCA

setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")
options(stringsAsFactors = FALSE);
library(WGCNA);
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(CorLevelPlot) 
library(gridExtra)
library(Rmisc) 

###Setting Up and Running DESeq
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

sample_metadata = read.table(file = "Ahud_trait_data.txt",header=T, row.names = 1)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=sample_metadata, 
                              design= ~ 1)

dim(dds)

# Filter out genes with too few reads - remove all genes with counts < 15 in more than 75% of samples, so ~28)
## suggested by WGCNA on RNAseq FAQ

dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 
# [1] 25260, that have at least 15 reads (a.k.a counts) in 75% of the samples

# Run the DESeq model to test for differential gene expression
dds <- DESeq(dds)


# 2. QC - outlier detection ------------------------------------------------
# detect outlier genes

gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK

table(gsg$goodGenes)
table(gsg$goodSamples)

# detect outlier samples - hierarchical clustering - method 1
htree <- hclust(dist(t(countsTable)), method = "average")
plot(htree) 

# pca - method 2

pca <- prcomp(t(countsTable))
pca.dat <- pca$x

pca.var <- pca$sdev^2
pca.var.percent <- round(pca.var/sum(pca.var)*100, digits = 2)

pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('PC1: ', pca.var.percent[1], ' %'),
       y = paste0('PC2: ', pca.var.percent[2], ' %'))

# All of our samples look good; the genes that don't look good shouldn't cluster well and 
# should end up in the "grey", catch all module for uncorrelated genes.

# 3. Normalization ----------------------------------------------------------------------

colData <- row.names(sample_metadata)

# making the rownames and column names identical
all(rownames(colData) %in% colnames(countsTableRound)) # to see if all samples are present in both
all(rownames(colData) == colnames(countsTableRound))  # to see if all samples are in the same order

# perform variance stabilization
dds_norm <- vst(dds)
#dds_norm <- vst(normalized_counts)

# get normalized counts
norm.counts <- assay(dds_norm) %>% 
  t()

# 4. Network Construction  ---------------------------------------------------
# Choose a set of soft-thresholding powers
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Call the network topology analysis function; this step takes a couple minutes
sft <- pickSoftThreshold(norm.counts,
                         powerVector = power,
                         networkType = "signed",
                         verbose = 5)

sft.data <- sft$fitIndices 

# visualization to pick power

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()


a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()


grid.arrange(a1, a2, nrow = 2)
# based on this plot, choose a soft power to maximize R^2 (above 0.8) and minimize connectivity
# for these ahud data: 6-8; Higher R2 should yield more modules.


# convert matrix to numeric
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power_6 <- 6 # should try running with 7 based on R^2
soft_power_8 <- 8
soft_power_4 <- 4
temp_cor <- cor
cor <- WGCNA::cor # use the 'cor' function from the WGCNA package

###WGNCA Visualizations

bwnet_4 <- readRDS("bwnet.4.rds")
bwnet_6 <- readRDS("bwnet.6.rds")
bwnet_8 <- readRDS("bwnet.8.rds")

###Bwent with threshold of 4 ---------------------------------------------------------
module_eigengenes_4 <- bwnet_4$MEs

head(module_eigengenes_4)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_4$dendrograms[[1]], cbind(bwnet_4$unmergedColors, bwnet_4$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

###Bwent with threshold of 6 ---------------------------------------------------------
module_eigengenes_6 <- bwnet_6$MEs

head(module_eigengenes_6)

# get number of genes for each module
table(bwnet_6$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_6$dendrograms[[1]], cbind(bwnet_6$unmergedColors, bwnet_6$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

###Bwent with threshold of 8 ---------------------------------------------------------
module_eigengenes_8 <- bwnet_8$MEs

head(module_eigengenes_8)

# get number of genes for each module
table(bwnet_8$colors)

# Plot the dendrogram and the module colors before and after merging underneath
plotDendroAndColors(bwnet_8$dendrograms[[1]], cbind(bwnet_8$unmergedColors, bwnet_8$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

##All Together
nets <- list("Power 4" = bwnet_4,
             "Power 6" = bwnet_6,
             "Power 8" = bwnet_8)

# Function: dendrogram + unmerged/merged module colors
plot_dendro <- function(net, title, block = 1) {
  # genes belonging to this block (handles multi-block runs safely)
  idx <- net$blockGenes[[block]]
  
  plotDendroAndColors(
    net$dendrograms[[block]],
    cbind(net$unmergedColors[idx], net$colors[idx]),
    groupLabels = c("Unmerged", "Merged"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste("Gene Dendrogram and Module Colors -", title)
  )
}

for (nm in names(nets)) {
  png(paste0("dendrogram_", gsub(" ", "_", nm), ".png"),
      width = 2400, height = 1200, res = 200)
  plot_dendro(nets[[nm]], nm)
  dev.off()
}

##Eigengene Visualization
ME4 <- bwnet_4$MEs
ME6 <- bwnet_6$MEs
ME8 <- bwnet_8$MEs

# Modules present at all three powers (drop grey, it's the unassigned bin)
shared <- Reduce(intersect, list(colnames(ME4), colnames(ME6), colnames(ME8)))
shared <- setdiff(shared, "MEgrey")
shared

# Eigengene signs are arbitrary, so flip any that correlate negatively with
# power 6 so the bars are comparable
align_sign <- function(ME, ref, mods) {
  for (m in mods) {
    if (WGCNA::cor(ME[[m]], ref[[m]]) < 0) ME[[m]] <- -ME[[m]]
  }
  ME
}
ME4a <- align_sign(ME4, ME6, shared)
ME8a <- align_sign(ME8, ME6, shared)

# Long format
to_long <- function(ME, power) {
  d <- as.data.frame(ME[, shared, drop = FALSE])
  d$Sample <- rownames(ME)
  d$Power  <- power
  tidyr::pivot_longer(d, cols = all_of(shared),
                      names_to = "Module", values_to = "Eigengene")
}

me_long <- rbind(to_long(ME4a, "4"),
                 to_long(ME6,  "6"),
                 to_long(ME8a, "8"))
me_long$Module <- sub("^ME", "", me_long$Module)

# 1. Per-sample eigengene values, one panel per module ------------------------
p_me <- ggplot(me_long, aes(x = Sample, y = Eigengene, fill = Power)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.85) +
  facet_wrap(~ Module, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Module eigengene values across WGCNA soft powers",
       x = "Sample", y = "Eigengene value", fill = "Soft power") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7))

p_me

ggsave("~/projects/eco_genomics_2026/transcriptomics/images/eigengene_values_by_power.png",
       p_me, width = 14, height = 3 * length(shared), dpi = 200, limitsize = FALSE)

# 2. How similar is each module's eigengene to power 6? -----------------------
cor_df <- rbind(
  data.frame(Module = shared,
             Comparison = "Power 4 vs 6",
             Correlation = sapply(shared, function(m) WGCNA::cor(ME4[[m]], ME6[[m]]))),
  data.frame(Module = shared,
             Comparison = "Power 8 vs 6",
             Correlation = sapply(shared, function(m) WGCNA::cor(ME8[[m]], ME6[[m]])))
)
cor_df$Module <- sub("^ME", "", cor_df$Module)

p_cor <- ggplot(cor_df, aes(x = Module, y = abs(Correlation), fill = Comparison)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  coord_cartesian(ylim = c(0, 1)) +
  theme_bw() +
  labs(title = "Eigengene similarity to power 6",
       x = "Module", y = "|Correlation| of eigengenes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_cor

ggsave("~/projects/eco_genomics_2026/transcriptomics/images/eigengene_similarity_barplot.png",
       p_cor, width = 10, height = 5, dpi = 250)


###Number of Genes across thresholds

get_module_sizes <- function(colors, threshold) {
  as.data.frame(table(colors)) |>
    dplyr::rename(Module = colors, Count = Freq) |>
    dplyr::mutate(Threshold = threshold)
}

df <- rbind(
  get_module_sizes(bwnet_4$colors, "4"),
  get_module_sizes(bwnet_6$colors, "6"),
  get_module_sizes(bwnet_8$colors, "8")
)

# Order modules by size at power 6 so the x-axis is easier to read
mod_order <- df$Module[df$Threshold == "6"][order(-df$Count[df$Threshold == "6"])]
df$Module <- factor(df$Module, levels = unique(c(as.character(mod_order),
                                                 as.character(df$Module))))

p <- ggplot(df, aes(x = Module, y = Count, fill = Threshold)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.75) +
  theme_bw() +
  labs(
    title = "Module Sizes Across WGCNA Thresholds",
    x = "Module Color",
    y = "Gene Count",
    fill = "Soft power"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p

ggsave("~/projects/eco_genomics_2026/transcriptomics/images/module_sizes_barplot.png",
       p, width = 12, height = 6, dpi = 250)

###Visualizing module correlation
library(pheatmap)

#4 v 6
corME_46 <- cor(bwnet_4$MEs,
                bwnet_6$MEs)

pheatmap(corME_46,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE)

#4 v 8
corME_48 <- cor(bwnet_4$MEs,
                bwnet_8$MEs)

pheatmap(corME_48,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE)

##All together
ME4 <- bwnet_4$MEs
ME6 <- bwnet_6$MEs
ME8 <- bwnet_8$MEs

# Add prefixes so modules are identifiable
colnames(ME4) <- paste0("P4_", colnames(ME4))
colnames(ME6) <- paste0("P6_", colnames(ME6))
colnames(ME8) <- paste0("P8_", colnames(ME8))

# Combine eigengenes
allMEs <- cbind(ME4, ME6, ME8)

# Correlation matrix
MEcor <- cor(allMEs)

library(pheatmap)

pheatmap(
  MEcor,
  color = colorRampPalette(c("blue","white","red"))(100),
  border_color = NA,
  main = "Eigengene Correlations Across Soft Threshold Powers"
)

#Clustering it
pheatmap(
  MEcor,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  color = colorRampPalette(c("navy","white","firebrick3"))(100)
)

#Gene Reassignment
table46 <- table(
  Power4 = bwnet_4$colors,
  Power6 = bwnet_6$colors
)

table48 <- table(bwnet_4$colors,
                 bwnet_8$colors)

table68 <- table(bwnet_6$colors,
                 bwnet_8$colors)

pheatmap(log10(table46 + 1))
pheatmap(log10(table48 + 1))
pheatmap(log10(table68 + 1))

#Pulling Out Some Interesting genes
best_match <- apply(abs(MEcor), 1, which.max)

data.frame(
  Module = rownames(MEcor),
  BestMatch = colnames(MEcor)[best_match],
  Correlation = apply(abs(MEcor), 1, max)
)

#
MEcor2 <- MEcor

diag(MEcor2) <- NA

best_match <- apply(abs(MEcor2), 1, which.max)

results <- data.frame(
  Module = rownames(MEcor2),
  BestMatch = colnames(MEcor2)[best_match],
  Correlation = apply(abs(MEcor2), 1, max, na.rm = TRUE)
)

results[order(-results$Correlation), ]

#                           Module        BestMatch Correlation
#P4_MEturquoise     P4_MEturquoise   P8_MEturquoise   0.9996560
#P8_MEturquoise     P8_MEturquoise   P4_MEturquoise   0.9996560
#P6_MEturquoise     P6_MEturquoise   P8_MEturquoise   0.9995236


### 6A. GO Association - looking at turquoise--------------------------------------------------
# module trait associations
module <- "turquoise"

# WGCNA module assignments
moduleColors <- bwnet_6$colors

# names(moduleColors) should be gene IDs
table(moduleColors)

geneList <- factor(
  as.integer(moduleColors == module)
)

names(geneList) <- sub(
  "::g\\.[0-9]+::m\\.[0-9]+$",
  "",
  names(bwnet_6$colors)
)

table(geneList)

##GO Analysis
library(topGO)

mappingFile <- "trinotate_annotation_GOblastx_forTopGO.txt"

geneID2GO <- readMappings(
  file = mappingFile
)

cat("Genes in GO mapping:",
    length(geneID2GO),
    "\n")

#Filter GO terms by size:(we dont want very small categories with small counts that 
#could be error, noise, or too specific, as well as super large ones (e.g. "biology" category))
allGO <- table(unlist(geneID2GO))

keepTerms <- names(allGO)[
  allGO >= 5 &
    allGO <= 500
]

geneID2GO.filtered <- lapply(
  geneID2GO,
  function(x) intersect(x, keepTerms)
)

geneID2GO.filtered <- geneID2GO.filtered[
  lengths(geneID2GO.filtered) > 0
]

GOdata <- new(
  "topGOdata",
  ontology = "BP",
  allGenes = geneList,
  geneSelectionFun = function(x) x == 1,
  annot = annFUN.gene2GO,
  gene2GO = geneID2GO.filtered
)

resultWeight <- runTest(
  GOdata,
  algorithm = "weight01",
  statistic = "fisher"
)

GOresults <- GenTable(
  GOdata,
  weightFisher = resultWeight,
  orderBy = "weightFisher",
  topNodes = 100
)

head(GOresults)

#GO.ID                                        Term Annotated Significant Expected weightFisher
#1 GO:0006936                          muscle contraction       236         162    90.97      1.5e-26
#2 GO:0007527            adult somatic muscle development        85          74    32.77      1.2e-21
#3 GO:0071689               muscle thin filament assembly        81          70    31.22      3.9e-20
#4 GO:0043520 regulation of myosin II filament assembl...        74          65    28.52      1.3e-19
#5 GO:0030241 skeletal muscle myosin thick filament as...        86          72    33.15      5.1e-19
#6 GO:0060361                                      flight        83          68    31.99      4.9e-17

##Top Hub Genes:
# Calculate module membership
vsd <- vst(dds)
datExpr <- t(assay(vsd))
MM <- cor(datExpr, bwnet_6$MEs, use = "p")

kME <- signedKME(
  datExpr,
  bwnet_6$MEs
)

hub.df <- data.frame(
  Gene = colnames(datExpr),
  Module = bwnet_6$colors,
  kME = kME$kMEturquoise
)

hub.df <- subset(
  hub.df,
  Module == "turquoise"
)

hub.df <- hub.df[
  order(-abs(hub.df$kME)),
]

head(hub.df, 5)

#                                                                     Module       kME
#TRINITY_DN2332_c0_g1::TRINITY_DN2332_c0_g1_i9::g.11448::m.11448    turquoise 0.9668181
#TRINITY_DN2332_c0_g1::TRINITY_DN2332_c0_g1_i3::g.11447::m.11447    turquoise 0.9668181
#TRINITY_DN2215_c0_g1::TRINITY_DN2215_c0_g1_i1::g.11114::m.11114    turquoise 0.9636212
#TRINITY_DN12290_c0_g1::TRINITY_DN12290_c0_g1_i13::g.37170::m.37170 turquoise 0.9613168
#TRINITY_DN26848_c0_g1::TRINITY_DN26848_c0_g1_i3::g.50201::m.50201  turquoise 0.9588781
# Is the gene in the mapping at all?
grep("DN26848", names(geneID2GO), value = TRUE)

get_go <- function(id) {
  go <- geneID2GO[[id]]
  if (is.null(go)) return(data.frame(ID = id, GO.ID = NA, Term = NA, Ontology = NA))
  data.frame(ID = id, GO.ID = go, Term = Term(go), Ontology = Ontology(go))
}

top5 <- head(hub.df, 5)
top5$ID <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", rownames(top5))

top5$ID %in% names(geneID2GO)   # should all be TRUE now

top5_go <- do.call(rbind, lapply(top5$ID, get_go))
top5_go$kME <- top5$kME[match(top5_go$ID, top5$ID)]
top5_go


#                                                           Term Ontology       kME
#1                                                               axon       CC 0.9668181
#2                                                 dense core granule       CC 0.9668181
#3                                                   exocytic vesicle       CC 0.9668181
#4                                                               <NA>     <NA> 0.9668181
#5                                                               <NA>     <NA> 0.9668181
#6                                             neuromuscular junction       CC 0.9668181
#7                                                    plasma membrane       CC 0.9668181
#8                                                   synaptic vesicle       CC 0.9668181
#9                                          synaptic vesicle membrane       CC 0.9668181
#10                                                   terminal bouton       CC 0.9668181

library(topGO)
library(dplyr)
library(tidyr)
library(ggplot2)

run_go <- function(net, power, module = "turquoise") {
  geneList <- factor(as.integer(net$colors == module))
  names(geneList) <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", names(net$colors))
  
  GOdata <- new(
    "topGOdata",
    ontology = "BP",
    allGenes = geneList,
    geneSelectionFun = function(x) x == 1,
    annot = annFUN.gene2GO,
    gene2GO = geneID2GO.filtered
  )
  
  res <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  # topNodes = all tested terms, so every term has a p-value at every power
  tab <- GenTable(GOdata, weightFisher = res, orderBy = "weightFisher",
                  topNodes = length(usedGO(GOdata)))
  
  tab$p     <- as.numeric(gsub("<", "", tab$weightFisher))
  tab$Power <- power
  tab
}

go4 <- run_go(bwnet_4, "4")
go6 <- run_go(bwnet_6, "6")
go8 <- run_go(bwnet_8, "8")

# Top 10 terms at each power
head(go4[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)
head(go8[, c("GO.ID", "Term", "Significant", "weightFisher")], 10)

# 1. How many turquoise genes at each power, and how much do they overlap? ----
tq <- list(
  P4 = names(bwnet_4$colors)[bwnet_4$colors == "turquoise"],
  P6 = names(bwnet_6$colors)[bwnet_6$colors == "turquoise"],
  P8 = names(bwnet_8$colors)[bwnet_8$colors == "turquoise"]
)
sapply(tq, length)

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))
jaccard(tq$P4, tq$P6)
jaccard(tq$P4, tq$P8)
jaccard(tq$P6, tq$P8)

# 2. Overlap of the top-N enriched terms ---------------------------------------
N <- 50
top_terms <- list(
  P4 = head(go4$GO.ID, N),
  P6 = head(go6$GO.ID, N),
  P8 = head(go8$GO.ID, N)
)

length(intersect(top_terms$P4, top_terms$P6))
length(intersect(top_terms$P4, top_terms$P8))
length(intersect(top_terms$P6, top_terms$P8))
length(Reduce(intersect, top_terms))   # in the top N at all three powers

# 3. Side-by-side table of p-values --------------------------------------------
go_all <- bind_rows(go4, go6, go8)

go_wide <- go_all %>%
  dplyr::select(GO.ID, Term, Power, p) %>%
  pivot_wider(names_from = Power, values_from = p, names_prefix = "p_power") %>%
  mutate(best_p = pmin(p_power4, p_power6, p_power8, na.rm = TRUE)) %>%
  arrange(best_p)

head(go_wide, 20)

# Terms enriched at one power but not the others (p < 0.05 at only one power)
sig <- go_wide %>%
  dplyr::mutate(across(starts_with("p_power"), ~ .x < 0.05, .names = "sig_{.col}")) %>%
  dplyr::mutate(n_sig = sig_p_power4 + sig_p_power6 + sig_p_power8)

sig %>%
  dplyr::filter(n_sig == 1) %>%
  dplyr::select(GO.ID, Term, starts_with("p_power"))

# 4. Plot: top terms across powers ---------------------------------------------
top_union <- go_wide %>% slice_head(n = 15) %>% pull(GO.ID)

plot_df <- go_all %>%
  dplyr::filter(GO.ID %in% top_union) %>%
  dplyr::mutate(neglogp = -log10(pmax(p, 1e-30)),     # cap so tiny p's don't dominate
         Term = factor(Term, levels = rev(unique(go_wide$Term[go_wide$GO.ID %in% top_union]))))

p_go <- ggplot(plot_df, aes(x = neglogp, y = Term, fill = Power)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  theme_bw() +
  labs(title = "Top GO terms in the turquoise module across soft powers",
       x = "-log10(weightFisher p-value)", y = NULL, fill = "Soft power")

p_go

ggsave("~/projects/eco_genomics_2026/images/GO_comparison_powers.png",
       p_go, width = 11, height = 7, dpi = 250)

