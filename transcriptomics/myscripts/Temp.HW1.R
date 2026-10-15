### Ecological Genomics - Homework #1 - 10/13/2026
### WGCNA: comparing soft-threshold powers 4, 6, 8

# 0. Setup ---------------------------------------------------------------------
setwd("~/projects/eco_genomics_2026/transcriptomics/mydata")
options(stringsAsFactors = FALSE)

# Load order matters: topGO/GO.db (AnnotationDbi) first, tidyverse last,
# so dplyr's select()/filter()/mutate() are not masked.
library(topGO)
library(GO.db)
library(WGCNA)
library(DESeq2)
library(gridExtra)
library(tidyverse)

img_dir <- "~/projects/eco_genomics_2026/transcriptomics/images"
dir.create(img_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Load data, filter, normalize ----------------------------------------------
countsTable      <- read.table("salmon.isoform.counts.matrix.filteredAssembly",
                               header = TRUE, row.names = 1)
countsTableRound <- round(countsTable)   # DESeq2 needs integers (Salmon gives decimals)
dim(countsTable)

sample_metadata <- read.table("Ahud_trait_data.txt", header = TRUE, row.names = 1)

# samples must be present and in the same order in both tables
all(rownames(sample_metadata) %in% colnames(countsTableRound))
all(rownames(sample_metadata) == colnames(countsTableRound))

dds <- DESeqDataSetFromMatrix(countData = countsTableRound,
                              colData   = sample_metadata,
                              design    = ~ 1)

# Keep genes with >= 15 reads in >= 28 samples (~75%), per the WGCNA RNA-seq FAQ
dds <- dds[rowSums(counts(dds) >= 15) >= 28, ]
nrow(dds)   # 25260

dds <- DESeq(dds)   # not strictly needed for vst(); kept so results are unchanged

# 2. QC: outlier genes and samples ---------------------------------------------
gsg <- goodSamplesGenes(t(countsTable))
summary(gsg)
gsg$allOK
table(gsg$goodGenes)
table(gsg$goodSamples)

# Sample clustering
plot(hclust(dist(t(countsTable)), method = "average"))

# PCA
pca             <- prcomp(t(countsTable))
pca.dat         <- as.data.frame(pca$x)
pca.var.percent <- round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0("PC1: ", pca.var.percent[1], " %"),
       y = paste0("PC2: ", pca.var.percent[2], " %"))
# Samples look good; poorly behaved genes should fall into the grey module.

# Variance-stabilized expression, samples in rows
norm.counts <- t(assay(vst(dds)))

# 3. Choosing soft-threshold powers --------------------------------------------
power <- c(1:10, seq(12, 50, by = 2))

sft <- pickSoftThreshold(norm.counts, powerVector = power,
                         networkType = "signed", verbose = 5)   # takes a few minutes
sft.data <- sft$fitIndices

a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = "red") +
  labs(x = "Power", y = "Scale free topology model fit, signed R^2") +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() + geom_text(nudge_y = 0.1) +
  labs(x = "Power", y = "Mean Connectivity") +
  theme_classic()

grid.arrange(a1, a2, nrow = 2)
# Pick the lowest power with R^2 > 0.8 and reasonable connectivity (here 6-8).

cor <- WGCNA::cor   # use WGCNA's cor()

# 4. Load networks (built earlier with blockwiseModules at powers 4, 6, 8) -------
bwnet_4 <- readRDS("bwnet.4.rds")
bwnet_6 <- readRDS("bwnet.6.rds")
bwnet_8 <- readRDS("bwnet.8.rds")

nets <- list("4" = bwnet_4, "6" = bwnet_6, "8" = bwnet_8)

# 5. Dendrograms ---------------------------------------------------------------
plot_dendro <- function(net, title, block = 1) {
  idx <- net$blockGenes[[block]]
  plotDendroAndColors(
    net$dendrograms[[block]],
    cbind(net$unmergedColors[idx], net$colors[idx]),
    groupLabels = c("Unmerged", "Merged"),
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE, guideHang = 0.05,
    main = paste("Gene Dendrogram and Module Colors -", title)
  )
}

# Save one PNG per power
for (p in names(nets)) {
  png(file.path(img_dir, paste0("dendrogram_Power_", p, ".png")),
      width = 2400, height = 1200, res = 200)
  plot_dendro(nets[[p]], paste("Power", p))
  dev.off()
}

# Optional: stack them into one image
library(magick)
imgs <- image_read(file.path(img_dir, paste0("dendrogram_Power_", names(nets), ".png")))
image_write(image_append(imgs, stack = TRUE),
            file.path(img_dir, "dendrograms_stacked.png"), format = "png")

# 6. How much do the networks differ across powers? ----------------------------
# Number of modules
sapply(nets, function(n) length(unique(n$colors)))   # 10, 11, 15

# Proportion of genes whose module label differs
mean(bwnet_4$colors != bwnet_6$colors)   # 0.195
mean(bwnet_4$colors != bwnet_8$colors)   # 0.281
mean(bwnet_6$colors != bwnet_8$colors)   # 0.303

# Module sizes
get_module_sizes <- function(colors, threshold) {
  as.data.frame(table(colors)) |>
    dplyr::rename(Module = colors, Count = Freq) |>
    dplyr::mutate(Threshold = threshold)
}

df <- rbind(get_module_sizes(bwnet_4$colors, "4"),
            get_module_sizes(bwnet_6$colors, "6"),
            get_module_sizes(bwnet_8$colors, "8"))

mod_order <- df$Module[df$Threshold == "6"][order(-df$Count[df$Threshold == "6"])]
df$Module <- factor(df$Module, levels = unique(c(as.character(mod_order),
                                                 as.character(df$Module))))

p_size <- ggplot(df, aes(x = Module, y = Count, fill = Threshold)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.75) +
  theme_bw() +
  labs(title = "Module Sizes Across WGCNA Thresholds",
       x = "Module Color", y = "Gene Count", fill = "Soft power") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_size
ggsave(file.path(img_dir, "module_sizes_barplot.png"), p_size,
       width = 12, height = 6, dpi = 250)

# Eigengene similarity of same-named modules (relative to power 6)
ME4 <- bwnet_4$MEs; ME6 <- bwnet_6$MEs; ME8 <- bwnet_8$MEs
shared <- setdiff(Reduce(intersect, list(colnames(ME4), colnames(ME6), colnames(ME8))),
                  "MEgrey")

cor_df <- rbind(
  data.frame(Module = shared, Comparison = "Power 4 vs 6",
             Correlation = sapply(shared, function(m) WGCNA::cor(ME4[[m]], ME6[[m]]))),
  data.frame(Module = shared, Comparison = "Power 8 vs 6",
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
ggsave(file.path(img_dir, "eigengene_similarity_barplot.png"), p_cor,
       width = 10, height = 5, dpi = 250)

# Best cross-power match for every module (ignores same-power modules)
allMEs <- cbind(setNames(ME4, paste0("P4_", colnames(ME4))),
                setNames(ME6, paste0("P6_", colnames(ME6))),
                setNames(ME8, paste0("P8_", colnames(ME8))))
MEcor <- WGCNA::cor(allMEs)

same_power <- outer(substr(rownames(MEcor), 1, 2),
                    substr(colnames(MEcor), 1, 2), "==")
MEcor_x <- abs(MEcor)
MEcor_x[same_power] <- NA

best <- apply(MEcor_x, 1, which.max)
results <- data.frame(Module      = rownames(MEcor_x),
                      BestMatch   = colnames(MEcor_x)[best],
                      Correlation = apply(MEcor_x, 1, max, na.rm = TRUE))
results[order(-results$Correlation), ]
# Turquoise is matched across powers with |r| ~ 0.9995+

# 7. GO enrichment setup -------------------------------------------------------
geneID2GO <- readMappings(file = "trinotate_annotation_GOblastx_forTopGO.txt")
cat("Genes in GO mapping:", length(geneID2GO), "\n")

# Drop very small (<5) and very large (>500) GO categories
allGO     <- table(unlist(geneID2GO))
keepTerms <- names(allGO)[allGO >= 5 & allGO <= 500]

geneID2GO.filtered <- lapply(geneID2GO, function(x) intersect(x, keepTerms))
geneID2GO.filtered <- geneID2GO.filtered[lengths(geneID2GO.filtered) > 0]

# Full GO table (all tested terms) for one module in one network
run_go <- function(net, power, module) {
  geneList <- factor(as.integer(net$colors == module))
  names(geneList) <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", names(net$colors))
  
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x) x == 1,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO.filtered)
  
  res <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  
  tab <- GenTable(GOdata, weightFisher = res, orderBy = "weightFisher",
                  topNodes = length(usedGO(GOdata)))
  tab$p     <- as.numeric(gsub("<", "", tab$weightFisher))
  tab$Power <- power
  tab
}

jaccard <- function(a, b) length(intersect(a, b)) / length(union(a, b))

# Everything for one module across powers: gene overlap, GO overlap, p-value table, plot
compare_module <- function(module, N = 50, n_top_plot = 15) {
  sizes <- vapply(nets, function(n) sum(n$colors == module), integer(1))
  if (any(sizes == 0))
    stop(module, " is absent at power(s): ", paste(names(sizes)[sizes == 0], collapse = ", "))
  
  go    <- setNames(lapply(names(nets), function(p) run_go(nets[[p]], p, module)), names(nets))
  genes <- lapply(nets, function(n) names(n$colors)[n$colors == module])
  
  prs       <- combn(names(nets), 2, simplify = FALSE)
  pair_name <- vapply(prs, function(x) paste0("P", x[1], " vs P", x[2]), character(1))
  
  # Gene-content overlap between powers (same color name != same module!)
  jac <- setNames(vapply(prs, function(x) jaccard(genes[[x[1]]], genes[[x[2]]]), numeric(1)),
                  pair_name)
  
  # Overlap of the top-N enriched GO terms
  top_terms    <- lapply(go, function(tab) head(tab$GO.ID, N))
  term_overlap <- setNames(vapply(prs, function(x)
    length(intersect(top_terms[[x[1]]], top_terms[[x[2]]])), integer(1)), pair_name)
  
  # Side-by-side p-values
  go_all  <- dplyr::bind_rows(go)
  pcols   <- paste0("p_power", names(nets))
  go_wide <- go_all |>
    dplyr::select(GO.ID, Term, Power, p) |>
    tidyr::pivot_wider(names_from = Power, values_from = p, names_prefix = "p_power")
  go_wide$best_p <- do.call(pmin, c(unname(as.list(go_wide[pcols])), na.rm = TRUE))
  go_wide        <- go_wide[order(go_wide$best_p), ]
  go_wide$n_sig  <- rowSums(go_wide[pcols] < 0.05, na.rm = TRUE)
  sig_one        <- go_wide[go_wide$n_sig == 1, c("GO.ID", "Term", pcols)]  # significant at one power only
  
  # Plot top terms
  top     <- head(go_wide, n_top_plot)
  plot_df <- go_all[go_all$GO.ID %in% top$GO.ID, ]
  plot_df$neglogp <- -log10(pmax(plot_df$p, 1e-30))   # cap tiny p-values
  plot_df$Term    <- factor(plot_df$Term, levels = rev(unique(top$Term)))
  
  p <- ggplot(plot_df, aes(x = neglogp, y = Term, fill = Power)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    theme_bw() +
    labs(title = paste("Top GO terms in the", module, "module across soft powers"),
         x = "-log10(weightFisher p-value)", y = NULL, fill = "Soft power")
  
  # Top 5 terms per power
  top5 <- do.call(rbind, lapply(go, function(tab) {
    x <- head(tab[, c("GO.ID", "Term", "Annotated", "Significant", "weightFisher")], 5)
    x$Module <- module
    x$Power  <- tab$Power[1]
    x$Rank   <- seq_len(nrow(x))
    x
  }))
  rownames(top5) <- NULL
  
  list(go = go, sizes = sizes, jaccard = jac, term_overlap = term_overlap,
       wide = go_wide, sig_one = sig_one, plot = p, top5 = top5)
}

# 8. Turquoise and yellow across powers ----------------------------------------
res_turq <- compare_module("turquoise")
res_yel  <- compare_module("yellow")
res_blue <- compare_module("blue")
res_brwn <- compare_module("brown")

#brown
res_brwn$sizes
res_brwn$jaccard
res_brwn$term_overlap
head(res_brwn$wide, 20)
res_brwn$plot

#blue
res_blue$sizes
res_blue$jaccard
res_blue$term_overlap
head(res_blue$wide, 20)
res_blue$plot

# Turquoise: stable (large module, muscle/flight signal)
res_turq$sizes
res_turq$jaccard
res_turq$term_overlap
head(res_turq$wide, 20)
res_turq$plot
ggsave(file.path(img_dir, "GO_comparison_turquoise.png"), res_turq$plot,
       width = 11, height = 7, dpi = 250)

# Yellow: small module, expected to be unstable across powers
res_yel$sizes
res_yel$jaccard         # low Jaccard = "yellow" is a different gene set at each power
res_yel$term_overlap
head(res_yel$wide, 20)
res_yel$plot
ggsave(file.path(img_dir, "GO_comparison_yellow.png"), res_yel$plot,
       width = 11, height = 7, dpi = 250)

# Top 5 GO terms per module and power
top5_all <- rbind(res_turq$top5, res_yel$top5)
top5_all[, c("Module", "Power", "Rank", "GO.ID", "Term", "Significant", "weightFisher")]
write.csv(top5_all, file.path(img_dir, "top5_GO_terms_by_module_power.csv"), row.names = FALSE)

# Where did the yellow genes from one power end up at another power?
fate <- function(from, to, module) {
  ids <- names(from$colors)[from$colors == module]
  sort(table(to$colors[ids]), decreasing = TRUE)
}
fate(bwnet_6, bwnet_4, "yellow")   # power 6 yellow -> power 4 modules
fate(bwnet_6, bwnet_8, "yellow")   # power 6 yellow -> power 8 modules
fate(bwnet_8, bwnet_6, "yellow")   # power 8 yellow -> power 6 modules





# 9. Hub genes (turquoise, power 6) --------------------------------------------
module_hubs <- function(net, module, n = 5) {
  stopifnot(identical(colnames(norm.counts), names(net$colors)))   # genes aligned?
  kME <- signedKME(norm.counts, net$MEs)
  hub <- data.frame(Module = net$colors, kME = kME[[paste0("kME", module)]])
  hub <- hub[hub$Module == module, ]
  head(hub[order(-abs(hub$kME)), ], n)
}

top5 <- module_hubs(bwnet_6, "turquoise")
top5$ID <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", rownames(top5))
top5

# Which hubs have a GO annotation? (DN12290 does not, so it has no GO terms)
top5$ID %in% names(geneID2GO)

get_go <- function(id) {
  go <- geneID2GO[[id]]
  if (is.null(go)) return(data.frame(ID = id, GO.ID = NA, Term = NA, Ontology = NA))
  data.frame(ID = id, GO.ID = go, Term = Term(go), Ontology = Ontology(go))
}

top5_go <- do.call(rbind, lapply(top5$ID, get_go))
top5_go$kME <- top5$kME[match(top5_go$ID, top5$ID)]
top5_go
# NA terms = GO IDs that are obsolete/not found in GO.db

### Module-trait and treatment/generation analysis ----------------------------
# Run after the streamlined script. Needs: bwnet_4/6/8, nets, sample_metadata,
# norm.counts, img_dir, and the tidyverse / WGCNA libraries.

# 1. Build a numeric trait table ------------------------------------------------
# (means only; the CI columns, Rep, and the duplicate generation label are dropped)
traits <- data.frame(
  Fitness    = sample_metadata$Fit_mean,
  EPR        = sample_metadata$EPR_mean,        # egg production rate
  HatchSucc  = sample_metadata$HS_mean,
  Survival   = sample_metadata$Surv_mean,
  DevTime    = sample_metadata$Dev_mean,
  Generation = sample_metadata$Generation,      # numeric: 0, 2, 4, 11
  row.names  = rownames(sample_metadata)
)

# One 0/1 column per treatment (AM, OA, OW, OWA)
trt <- model.matrix(~ 0 + treatment, data = sample_metadata)
colnames(trt) <- sub("^treatment", "trt_", colnames(trt))
traits <- cbind(traits, trt)

table(sample_metadata$treatment, sample_metadata$Generation)   # check the design
cor(traits[, 1:6])                                              # traits are correlated with each other

# 2. Module-trait correlations at each power ------------------------------------
module_trait <- function(net, power) {
  MEs <- net$MEs[, colnames(net$MEs) != "MEgrey", drop = FALSE]
  tr  <- traits[rownames(MEs), , drop = FALSE]          # align samples by name
  
  r <- WGCNA::cor(MEs, tr, use = "pairwise.complete.obs")
  p <- corPvalueStudent(r, nSamples = nrow(tr))
  
  out <- as.data.frame(as.table(r))
  names(out) <- c("Module", "Trait", "r")
  out$Module <- sub("^ME", "", as.character(out$Module))
  out$Trait  <- as.character(out$Trait)
  out$p      <- as.vector(p)
  out$Power  <- power
  out
}

mt_all <- do.call(rbind, lapply(names(nets), function(p) module_trait(nets[[p]], p)))

# BH-adjust within each power; ** = adjusted p < 0.05, * = raw p < 0.05
mt_all$padj  <- ave(mt_all$p, mt_all$Power, FUN = function(x) p.adjust(x, "BH"))
mt_all$label <- paste0(sprintf("%.2f", mt_all$r),
                       ifelse(mt_all$padj < 0.05, "\n**", ifelse(mt_all$p < 0.05, "\n*", "")))

plot_mt <- function(d, title) {
  d$Trait <- factor(d$Trait, levels = colnames(traits))
  ggplot(d, aes(Trait, Module, fill = r)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), size = 2.8, lineheight = 0.9) +
    scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", limits = c(-1, 1)) +
    theme_minimal() +
    labs(title = title, x = NULL, y = NULL, fill = "r") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

p_mt4 <- plot_mt(subset(mt_all, Power == "4"), "Module-trait correlations (power 4)")
p_mt4

p_mt6 <- plot_mt(subset(mt_all, Power == "6"), "Module-trait correlations (power 6)")
p_mt6

p_mt8 <- plot_mt(subset(mt_all, Power == "8"), "Module-trait correlations (power 8)")
p_mt8


# Is the turquoise-trait relationship robust to power? (turquoise is the one module
# that is the same module at every power; other colors are NOT comparable by name)
mt_all |>
  dplyr::filter(Module == "turquoise") |>
  dplyr::select(Trait, Power, r) |>
  tidyr::pivot_wider(names_from = Power, values_from = r, names_prefix = "P")

# 3. Eigengenes by treatment and generation (power 6) ---------------------------
MEs6 <- bwnet_6$MEs
meta <- sample_metadata[rownames(MEs6), ]

# Tests per module: treatment effect and generation effect (each adjusted for the other).
# Exploratory only: the design is unbalanced (F11 only in AM and OWA; AM F2 has 2 reps)
mod_tests <- do.call(rbind, lapply(setdiff(colnames(MEs6), "MEgrey"), function(m) {
  d  <- data.frame(ME = MEs6[[m]],
                   treatment = factor(meta$treatment),
                   gen = factor(meta$Generation))
  d1 <- drop1(lm(ME ~ treatment + gen, data = d), test = "F")
  data.frame(Module = sub("^ME", "", m),
             p_treatment  = d1["treatment", "Pr(>F)"],
             p_generation = d1["gen", "Pr(>F)"])
}))
mod_tests$padj_treatment  <- p.adjust(mod_tests$p_treatment,  "BH")
mod_tests$padj_generation <- p.adjust(mod_tests$p_generation, "BH")
mod_tests[order(mod_tests$p_treatment), ]

# Trajectories: eigengene across generations, one line per treatment
key_mods <- c("turquoise", "yellow", "blue", "brown")   # edit to whichever modules stand out in step 2

me_long <- data.frame(treatment = meta$treatment, Generation = meta$Generation, MEs6) |>
  tidyr::pivot_longer(dplyr::starts_with("ME"), names_to = "Module", values_to = "Eigengene") |>
  dplyr::mutate(Module = sub("^ME", "", Module)) |>
  dplyr::filter(Module %in% key_mods)

me_sum <- me_long |>
  dplyr::group_by(Module, treatment, Generation) |>
  dplyr::summarise(mean = mean(Eigengene),
                   se   = sd(Eigengene) / sqrt(dplyr::n()),
                   .groups = "drop")

p_traj <- ggplot(me_sum, aes(Generation, mean, color = treatment)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(data = me_long, aes(y = Eigengene), alpha = 0.3, size = 1) +
  scale_x_continuous(breaks = c(0, 2, 4, 11)) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw() +
  labs(title = "Module eigengenes across generations by treatment (power 6)",
       y = "Eigengene (mean +/- SE, points = replicates)")
p_traj
ggsave(file.path(img_dir, "eigengene_trajectories.png"), p_traj, width = 11, height = 7, dpi = 250)

##across different powers now:
# 3. Eigengenes by treatment and generation (power 6) ---------------------------
MEs4 <- bwnet_4$MEs
meta <- sample_metadata[rownames(MEs4), ]

# Tests per module: treatment effect and generation effect (each adjusted for the other).
# Exploratory only: the design is unbalanced (F11 only in AM and OWA; AM F2 has 2 reps)
mod_tests <- do.call(rbind, lapply(setdiff(colnames(MEs4), "MEgrey"), function(m) {
  d  <- data.frame(ME = MEs4[[m]],
                   treatment = factor(meta$treatment),
                   gen = factor(meta$Generation))
  d1 <- drop1(lm(ME ~ treatment + gen, data = d), test = "F")
  data.frame(Module = sub("^ME", "", m),
             p_treatment  = d1["treatment", "Pr(>F)"],
             p_generation = d1["gen", "Pr(>F)"])
}))
mod_tests$padj_treatment  <- p.adjust(mod_tests$p_treatment,  "BH")
mod_tests$padj_generation <- p.adjust(mod_tests$p_generation, "BH")
mod_tests[order(mod_tests$p_treatment), ]

# Trajectories: eigengene across generations, one line per treatment
key_mods <- c("turquoise", "yellow", "blue", "brown")   # edit to whichever modules stand out in step 2

me_long <- data.frame(treatment = meta$treatment, Generation = meta$Generation, MEs4) |>
  tidyr::pivot_longer(dplyr::starts_with("ME"), names_to = "Module", values_to = "Eigengene") |>
  dplyr::mutate(Module = sub("^ME", "", Module)) |>
  dplyr::filter(Module %in% key_mods)

me_sum <- me_long |>
  dplyr::group_by(Module, treatment, Generation) |>
  dplyr::summarise(mean = mean(Eigengene),
                   se   = sd(Eigengene) / sqrt(dplyr::n()),
                   .groups = "drop")

p_traj <- ggplot(me_sum, aes(Generation, mean, color = treatment)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(data = me_long, aes(y = Eigengene), alpha = 0.3, size = 1) +
  scale_x_continuous(breaks = c(0, 2, 4, 11)) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw() +
  labs(title = "Module eigengenes across generations by treatment (power 4)",
       y = "Eigengene (mean +/- SE, points = replicates)")
p_traj
ggsave(file.path(img_dir, "eigengene_trajectories.png"), p_traj, width = 11, height = 7, dpi = 250)

##
MEs8 <- bwnet_8$MEs
meta <- sample_metadata[rownames(MEs8), ]

# Tests per module: treatment effect and generation effect (each adjusted for the other).
# Exploratory only: the design is unbalanced (F11 only in AM and OWA; AM F2 has 2 reps)
mod_tests <- do.call(rbind, lapply(setdiff(colnames(MEs8), "MEgrey"), function(m) {
  d  <- data.frame(ME = MEs8[[m]],
                   treatment = factor(meta$treatment),
                   gen = factor(meta$Generation))
  d1 <- drop1(lm(ME ~ treatment + gen, data = d), test = "F")
  data.frame(Module = sub("^ME", "", m),
             p_treatment  = d1["treatment", "Pr(>F)"],
             p_generation = d1["gen", "Pr(>F)"])
}))
mod_tests$padj_treatment  <- p.adjust(mod_tests$p_treatment,  "BH")
mod_tests$padj_generation <- p.adjust(mod_tests$p_generation, "BH")
mod_tests[order(mod_tests$p_treatment), ]

# Trajectories: eigengene across generations, one line per treatment
key_mods <- c("turquoise", "yellow", "blue", "brown")   # edit to whichever modules stand out in step 2

me_long <- data.frame(treatment = meta$treatment, Generation = meta$Generation, MEs8) |>
  tidyr::pivot_longer(dplyr::starts_with("ME"), names_to = "Module", values_to = "Eigengene") |>
  dplyr::mutate(Module = sub("^ME", "", Module)) |>
  dplyr::filter(Module %in% key_mods)

me_sum <- me_long |>
  dplyr::group_by(Module, treatment, Generation) |>
  dplyr::summarise(mean = mean(Eigengene),
                   se   = sd(Eigengene) / sqrt(dplyr::n()),
                   .groups = "drop")

p_traj <- ggplot(me_sum, aes(Generation, mean, color = treatment)) +
  geom_line() +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
  geom_point(data = me_long, aes(y = Eigengene), alpha = 0.3, size = 1) +
  scale_x_continuous(breaks = c(0, 2, 4, 11)) +
  facet_wrap(~ Module, scales = "free_y") +
  theme_bw() +
  labs(title = "Module eigengenes across generations by treatment (power 8)",
       y = "Eigengene (mean +/- SE, points = replicates)")
p_traj



# 4. Candidate genes: module membership vs gene-trait correlation ----------------
focus_mod   <- "turquoise"   # edit after looking at the heatmap
focus_trait <- "Fitness"

stopifnot(identical(colnames(norm.counts), names(bwnet_6$colors)))
kME <- signedKME(norm.counts, bwnet_6$MEs)[[paste0("kME", focus_mod)]]
GS  <- as.vector(WGCNA::cor(norm.counts, traits[rownames(norm.counts), focus_trait]))

gene_df <- data.frame(Gene = colnames(norm.counts), Module = bwnet_6$colors, kME = kME, GS = GS)
gene_df <- gene_df[gene_df$Module == focus_mod, ]

cor.test(gene_df$kME, abs(gene_df$GS))   # do the best-connected genes track the trait most?

p_mmgs <- ggplot(gene_df, aes(kME, abs(GS))) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(title = paste(focus_mod, "module: membership vs", focus_trait, "correlation"),
       x = paste0("kME (", focus_mod, ")"), y = paste("|gene-trait correlation|:", focus_trait))
p_mmgs
ggsave(file.path(img_dir, "MM_vs_GS.png"), p_mmgs, width = 7, height = 5, dpi = 250)

# Genes high in both: candidates to annotate with the Trinotate BLAST hits
head(gene_df[order(-abs(gene_df$kME) * abs(gene_df$GS)), ], 10)

# Does yellow separate warming from non-warming at F0 and at F4?
for (g in c(0, 4)) {
  keep <- meta$Generation == g
  warm <- meta$treatment[keep] %in% c("OW", "OWA")
  cat("\nGeneration", g, "\n")
  print(t.test(bwnet_6$MEs$MEyellow[keep] ~ warm))
}

# 3. What is power-6 yellow?
head(run_go(bwnet_6, "6", "yellow")[, c("Term", "Annotated", "Significant", "weightFisher")], 15)
module_hubs(bwnet_6, "yellow", n = 10)