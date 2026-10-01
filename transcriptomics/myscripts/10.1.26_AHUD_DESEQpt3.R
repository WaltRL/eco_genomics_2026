setwd ("~/projects/eco_genomics_2026/transcriptomics/mydata")
getwd()

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)

# Import the counts matrix
countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

# Import the counts matrix
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)

dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)
dim(dds)

# Filter
dds <- dds[rowSums(counts(dds) >= 15) >= 28,] #looking for gene coverage of over 15 instances in 75% of samples (28)
nrow(dds)

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

####################################################

### Check on the DE results from the DESeq 

####################################################

resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05) #alpha is p-value cut off
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),] #taking that output from above and ordering it by p adjusted value
head(res_OWAvsAM) #top significantly differentiated 
summary(res_OWAvsAM)
#ut of 25260 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2343, 9.3%
#LFC < 0 (down)     : 1575, 6.2%
#outliers [1]       : 22, 0.087%
#low counts [2]     : 979, 3.9%
#(mean count < 23)

res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),]
head(res_OWvsAM) 
summary(res_OWvsAM)
#ut of 25260 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2343, 9.3%
#LFC < 0 (down)     : 1575, 6.2%
#outliers [1]       : 22, 0.087%
#low counts [2]     : 979, 3.9%
#(mean count < 23)

res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM <- res_OAvsAM[order(res_OAvsAM$padj),]
head(res_OAvsAM)
summary(res_OAvsAM)
#Out of 25260 with nonzero total read count
#adjusted p-value < 0.05
#LFC > 0 (up)       : 2343, 9.3%
#LFC < 0 (down)     : 1575, 6.2%
#outliers [1]       : 22, 0.087%
#low counts [2]     : 979, 3.9%
#(mean count < 23)


#################################################################

#### Scatter plot to assess how correlated are responses to OWA vs OW?

#################################################################


# Create merged data frame - need to use rownames because differences in filtering
plot_OWA <- data.frame(
  gene = rownames(res_OWAvsAM),
  LFC_OWA = res_OWAvsAM$log2FoldChange,
  padj_OWA = res_OWAvsAM$padj
)

plot_OW <- data.frame(
  gene = rownames(res_OWvsAM),
  LFC_OW = res_OWvsAM$log2FoldChange,
  padj_OW = res_OWvsAM$padj
)

plot_df <- merge(plot_OWA,
                 plot_OW,
                 by = "gene")

# Remove genes with missing LFC values
plot_df <- plot_df %>%
  filter(!is.na(LFC_OWA),
         !is.na(LFC_OW))

# Classify significance
plot_df <- plot_df %>%
  mutate(
    SigGroup = case_when(
      padj_OWA < 0.05 & padj_OW < 0.05 ~ "Both",
      padj_OWA < 0.05 ~ "OWA only",
      padj_OW < 0.05 ~ "OW only",
      TRUE ~ "Neither"
    )
  )

# Correlation for noting on the plot 
r <- cor(plot_df$LFC_OWA,
         plot_df$LFC_OW,
         use = "complete.obs")

# Arrange the genes by significant to make the plotting easier/more interesting to see
# ggplot plots in the order of the df, so random

plot_df$SigGroup <- factor(
  plot_df$SigGroup,
  levels = c("Neither", "OWA only", "OW only", "Both")
)

plot_df <- plot_df %>%
  arrange(SigGroup)

# Now make the plot!

ggplot(plot_df,
       aes(x = LFC_OW,
           y = LFC_OWA,
           color = SigGroup)) +
  
  geom_point(alpha = 0.6, size = 1.5) +
  
  geom_abline(intercept = 0,
              slope = 1,
              linetype = "dashed",
              color = "black") +
  
  geom_hline(yintercept = 0,
             color = "grey70") +
  
  geom_vline(xintercept = 0,
             color = "grey70") +
  
  annotate("text",
           x = min(plot_df$LFC_OW, na.rm = TRUE),
           y = max(plot_df$LFC_OWA, na.rm = TRUE),
           hjust = 0,
           label = paste0("r = ", round(r, 3))) +
  
  scale_color_manual(values = c(
    "Both" = "purple",
    "OWA only" = "#CC3333",
    "OW only" = "#00A08A",
    "Neither" = "grey80"
  )) +
  
  coord_fixed() + # forces the same scaling on x and y axes
  
  labs(
    x = "Log2 Fold Change: OW vs AM",
    y = "Log2 Fold Change: OWA vs AM",
    color = "",
    title = "GE Responses to OW relative to OWA"
  ) +
  
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    legend.position = "right"
  )
