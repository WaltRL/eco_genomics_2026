setwd ("~/projects/eco_genomics_2026/transcriptomics/mydata/myresults")

library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)

countsTable <- read.table("salmon.isoform.counts.matrix.filteredAssembly", header=TRUE, row.names=1)
head(countsTable)
dim(countsTable)

countsTableRound <- round(countsTable) # bc DESeq2 doesn't like decimals (and Salmon outputs data with decimals)
head(countsTableRound)

#import the sample description table
conds <- read.delim("ahud_samples_R.txt", header=TRUE, stringsAsFactors = TRUE, row.names=1)
head(conds)


dds <- DESeqDataSetFromMatrix(countData = countsTableRound, colData=conds, 
                              design= ~ treatment)

dim(dds)
# [1] 130580     38

# Filter 
dds <- dds[rowSums(counts(dds) >= 15) >= 28,]
nrow(dds) 

# Subset the DESeqDataSet to the specific level of the "generation" factor
dds_F0 <- subset(dds, select = generation == 'F0')
dim(dds_F0)
# [1] 25260    12

# Perform DESeq2 analysis on the subset
dds_F0 <- DESeq(dds_F0)

resultsNames(dds_F0)
# [1] "Intercept"           "treatment_OA_vs_AM"  "treatment_OW_vs_AM"  "treatment_OWA_vs_AM"

res_OWAvsAM <- results(dds_F0, name="treatment_OWA_vs_AM", alpha=0.05)
res_OWAvsAM <- res_OWAvsAM[order(res_OWAvsAM$padj),]
head(res_OWAvsAM)  
summary(res_OWAvsAM)


res_OWvsAM <- results(dds_F0, name="treatment_OW_vs_AM", alpha=0.05)
res_OWvsAM <- res_OWvsAM[order(res_OWvsAM$padj),]
head(res_OWvsAM) 
summary(res_OWvsAM)

res_OAvsAM <- results(dds_F0, name="treatment_OA_vs_AM", alpha=0.05)
res_OAvsAM <- res_OAvsAM[order(res_OAvsAM$padj),]


### FOR TODAY:
setwd ("~/projects/eco_genomics_2026/transcriptomics/mydata")

res_OWAvsAM.df <- as.data.frame(res_OWAvsAM)
res_OWAvsAM.df$fullID <- rownames(res_OWAvsAM.df)

parts <- strsplit(res_OWAvsAM.df$fullID, "::")

res_OWAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWAvsAM.df,
  "myresults/F0_OWAvsAM_results.csv",
  row.names = FALSE
)

###
res_OWvsAM.df <- as.data.frame(res_OWvsAM)
res_OWvsAM.df$fullID <- rownames(res_OWvsAM.df)

parts <- strsplit(res_OWvsAM.df$fullID, "::")

res_OWvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OWvsAM.df,
  "myresults/F0_OWvsAM_results.csv",
  row.names = FALSE
)

###
res_OAvsAM.df <- as.data.frame(res_OAvsAM)
res_OAvsAM.df$fullID <- rownames(res_OAvsAM.df)

parts <- strsplit(res_OAvsAM.df$fullID, "::")

res_OAvsAM.df$shortID <- sapply(
  parts,
  function(x) paste(x[1:2], collapse="::")
)

write.csv(
  res_OAvsAM.df,
  "myresults/F0_OAvsAM_results.csv",
  row.names = FALSE
)

###
#copied over files from cl/biol3990

#prep to run GO analyses
library(topGO)

setwd ("~/projects/eco_genomics_2026/transcriptomics/mydata")

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

#Create a function to run a TopGO contrast
run_topGO_contrast <- function(
    infile,
    outfile,
    ontology = "BP",
    padj.cutoff = 0.05){
  
  deseq <- read.csv(
    infile,
    stringsAsFactors = FALSE
  )
  
  deseq <- subset(
    deseq,
    !is.na(padj)
  )
  
  geneList <- factor(
    as.integer(deseq$padj < padj.cutoff)
  )
  
  names(geneList) <- deseq$shortID
  
  cat("\nGenes tested:",
      length(geneList))
  
  cat("\nSignificant genes:",
      sum(geneList == 1),
      "\n")
  
  GOdata <- new(
    "topGOdata",
    ontology = ontology,
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
  
  write.csv(
    GOresults,
    outfile,
    row.names = FALSE
  )
  
  return(GOresults)
}

#Run our TopGO function for the three different contrasts
GO_OA <- run_topGO_contrast(
  "myresults/F0_OAvsAM_results.csv",
  "myresults/GO_F0_OAvsAM_BP.csv"
)

GO_OW <- run_topGO_contrast(
  "myresults/F0_OWvsAM_results.csv",
  "myresults/GO_F0_OWvsAM_BP.csv"
)

GO_OWA <- run_topGO_contrast(
  "myresults/F0_OWAvsAM_results.csv",
  "myresults/GO_F0_OWAvsAM_BP.csv"
)

#Visualizing GO

# Read TopGO results
setwd ("~/projects/eco_genomics_2026/transcriptomics/mydata/myresults")
list.files()



go <- read.csv(
  "GO_F0_OWAvsAM_BP.csv",
  stringsAsFactors = FALSE
)

# Convert p-values to numeric
#AI for this line:
# Force weightFisher into a plain character vector

go$weightFisher <- gsub("^<\\s*", "", go$weightFisher)

go$weightFisher <- as.numeric(go$weightFisher)

# Convert counts to numeric
go$Significant <- as.numeric(go$Significant)

# Create -log10(p)
go$minusLogP <- -log10(go$weightFisher)

head(go)
str(go)

# Keep top 10 GO terms
go_top10 <- go %>%
  arrange(weightFisher) %>%
  head(n=10)

# Order terms for plotting
go_top10$Term <- factor(
  go_top10$Term,
  levels = rev(go_top10$Term)
)

# Bubble plot
ggplot(
  go_top10,
  aes(
    x = minusLogP,
    y = Term
  )
) +
  geom_point(
    aes(
      size = Significant,
      color = minusLogP
    )
  ) +
  scale_color_viridis_c() +
  theme_bw(base_size = 14) +
  labs(
    title = "Top GO Terms: OWA vs AM",
    x = expression(-logp),
    y = "GO Term",
    color = expression(-logp),
    size = "Significant\nGenes"
  )

###DO FOR OTHER COMPARISIONS

