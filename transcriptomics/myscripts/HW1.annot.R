### Ecological Genomics - Homework #1 - 10/13/2026
### WGCNA: comparing soft-threshold powers 4, 6, 8
###
### ############################################################################
### PIPELINE MAP  (run top to bottom; later stages depend on objects made earlier)
### ############################################################################
###
### STAGE A | Data prep
###   0  Setup ............... libraries, paths, img_dir
###   1  Load/filter/normalize . counts + metadata -> dds, norm.counts
###   2  QC .................... outlier genes/samples (diagnostic only)
###
### STAGE B | Network inputs
###   3  Power selection ....... pickSoftThreshold -> scale-free fit / connectivity plots
###   4  Load networks ......... bwnet_4/6/8 from .rds -> nets (named list "4","6","8")
###
### STAGE C | Compare networks across powers
###   5  Dendrograms ........... PNGs + stacked image
###   6  Network differences ... module counts, relabeling, sizes, eigengene similarity
###
### STAGE D | Function (GO) and hub genes
###   7  GO setup .............. geneID2GO.filtered, run_go(), jaccard(), compare_module()
###   8  GO across powers ...... turquoise / yellow / blue / brown
###   9  Hub genes ............. module_hubs(), GO terms of top turquoise hubs
###
### STAGE E | Trait integration
###  10  Trait table ........... traits (numeric traits + treatment dummies)
###  11  Module-trait cor ...... mt_all, heatmaps per power
###  12  Treatment/generation .. eigengene tests + trajectories (powers 6, 4, 8)
###  13  Candidate genes ....... kME vs gene-trait correlation (MM vs GS)
###  14  Yellow follow-up ...... warm vs non-warm tests, GO, hubs
###  15  Methods summary ....... versions, sample/gene counts, parameters (for the report)
###
### EXTERNAL INPUTS (all in ~/projects/eco_genomics_2026/transcriptomics/mydata):
###   salmon.isoform.counts.matrix.filteredAssembly   raw isoform counts (Salmon)
###   Ahud_trait_data.txt                             sample metadata + traits
###   bwnet.4.rds, bwnet.6.rds, bwnet.8.rds           blockwiseModules results (BUILT ELSEWHERE)
###   trinotate_annotation_GOblastx_forTopGO.txt      gene -> GO mapping (topGO format)
###
### OUTPUTS (all written to img_dir):
###   dendrogram_Power_{4,6,8}.png, dendrograms_stacked.png
###   module_sizes_barplot.png, eigengene_similarity_barplot.png
###   GO_comparison_turquoise.png, GO_comparison_yellow.png
###   top5_GO_terms_by_module_power.csv
###   eigengene_trajectories.png (see NOTE in section 12), MM_vs_GS.png
###   methods_summary.txt, sessionInfo.txt, network_summary_by_power.csv (section 15)
###
### SLOW STEPS: DESeq(dds) (section 1), pickSoftThreshold (3), each run_go call
###   and therefore each compare_module() call (8), signedKME on all genes (9, 13)
###
### ADD-ON SCRIPTS that slot in after this one (they reuse its objects):
###   turquoise_children_P6_to_P8.R ... after section 8 (needs run_go, geneID2GO.filtered)
###   kME_across_powers.R ............. after section 7 (needs nets, norm.counts, jaccard)
###   trait_correlations_across_powers.R  after section 11 (needs mt_all, plot_mt)
### ############################################################################


### ############################################################################
### PACKAGES USED (what each one does in this script)
### ############################################################################
### topGO      Gene Ontology (GO) enrichment testing.
###            readMappings, annFUN.gene2GO, new("topGOdata"), runTest (weight01 + Fisher),
###            GenTable, usedGO. Sections 7-9 and 14.
### GO.db      Offline GO database (term names, ontologies). Term(), Ontology() in section 9.
###            topGO also needs it to build the GO graph.
### WGCNA      Weighted gene co-expression network analysis.
###            goodSamplesGenes (2), pickSoftThreshold (3), plotDendroAndColors (5),
###            signedKME (9, 13), corPvalueStudent (11), WGCNA::cor, a fast correlation
###            that handles NAs (6, 11, 13). The bwnet_*.rds objects were made offline with
###            blockwiseModules (not in this script).
### DESeq2     Count-data handling and normalization.
###            DESeqDataSetFromMatrix, counts, DESeq, vst (section 1). assay() comes from
###            SummarizedExperiment, which DESeq2 attaches.
### gridExtra  grid.arrange(): puts the two soft-power plots in one figure (section 3).
### tidyverse  Meta-package. Parts used: ggplot2 (all plots), dplyr (bind_rows, select, filter,
###            mutate, group_by, summarise, rename), tidyr (pivot_wider, pivot_longer).
###            Loaded LAST so its verbs are not masked by AnnotationDbi (via topGO/GO.db).
### magick     ImageMagick wrapper: image_read / image_append / image_write stack the
###            dendrogram PNGs (section 5). Needs the ImageMagick system library; this
###            step is optional if it will not install.
### Base R / stats (no library call): read.table, round, t, hclust, dist, prcomp, cor, lm,
###            drop1, t.test, p.adjust, model.matrix, ave, outer, combn, table, split,
###            apply/sapply/lapply. The |> pipe is base R (R >= 4.1), not magrittr's %>%.
### Masking to remember: after `cor <- WGCNA::cor` (section 3), cor() is WGCNA's version.
###            If tidyverse were loaded BEFORE topGO, AnnotationDbi::select would hide dplyr::select.
### ############################################################################

### ############################################################################
### DATA TRANSFORMATION TRACE (how the data changes shape and scale)
### ############################################################################
### Raw Salmon counts (genes x samples, decimals, isoform-level) ........ countsTable
###   -> round() ........... integer counts, because DESeq2 requires integers
###   -> DESeqDataSetFromMatrix() counts + sample metadata in one object (design ~ 1) ... dds
###   -> row filter ........ keep genes with >= 15 reads in >= 28 samples (25,260 remain)
###   -> DESeq() ........... size factors (library-size normalization) + dispersions
###   -> vst() ............. variance-stabilizing transformation: roughly log2 scale, corrects
###                          for library size, removes the mean-variance dependence of counts
###   -> assay() ........... extract the numeric matrix (genes x samples)
###   -> t() ............... transpose to samples x genes, the orientation WGCNA expects
###                          ........................................ norm.counts
###
### Network construction (done offline by blockwiseModules; makes the .rds files)
###   gene-gene correlations
###   -> adjacency: |cor|^power (unsigned) or ((1 + cor)/2)^power (signed), depending on the
###      networkType that was used. Higher power shrinks weak correlations toward 0
###   -> topological overlap matrix (TOM): how much two genes share the same neighbours
###   -> 1 - TOM distance -> hierarchical clustering -> dynamic tree cut -> modules
###   -> merge modules whose eigengenes are too similar (mergeCutHeight)
###   Result: $colors (gene -> module), $MEs (eigengenes), $dendrograms, $unmergedColors
###
### Quantities derived downstream
###   module eigengene (ME)   first principal component of a module's genes: one expression
###                           profile per module across samples ($MEs, samples x modules;
###                           MEgrey = genes assigned to no module)
###   kME (module membership) correlation of each gene's expression with a module eigengene
###   gene significance (GS)  correlation of each gene's expression with a trait (section 13)
###   Jaccard index           |A and B| / |A or B| for two gene sets; module overlap that does
###                           not depend on color names
###   -log10(p)               GO p-values on a log scale for plotting (floored at 1e-30)
### ############################################################################


### ############################################################################
### REPORT GUIDE: using this script for the "Bioinformatics Pipeline" section
### ############################################################################
### Model: Example2_TranscriptomicsHW.pdf (Background, Bioinformatics Pipeline, Results,
### Conclusion, Figures, References). In that example's rubric Methods is worth 30 of 100
### points, the largest single item, and its pipeline paragraphs name every tool, version,
### threshold, and sample count. Each section below has a "REPORT NOTES" block with:
###   REPORT   which part of the paper the section feeds
###   METHODS  a draft sentence; [square brackets] = a value to fill in. Section 15 prints
###            the values that can come from R
###   FIGURE   saved files that could become figures (and plots worth saving)
###
### Suggested order of the pipeline paragraphs (mirrors the example):
###   1 read cleaning + QC ........ NOT in this script (e.g. fastp / MultiQC)
###   2 mapping + quantification .. NOT in this script (Salmon)
###   3 import, filter, normalize . section 1
###   4 sample QC ................. section 2
###   5 soft power + network ...... sections 3-4
###   6 comparing powers .......... sections 5-6 and 8
###   7 functional enrichment ..... sections 7-9
###   8 traits and treatment ...... sections 10-14
###
### INFORMATION THAT LIVES OUTSIDE THIS SCRIPT (collect it before writing):
###   - species, collection site, rearing design, treatment definitions (AM/OA/OW/OWA?),
###     generations, replicates per group. The file names say "ahud", which suggests
###     Acartia hudsonica (the example paper used A. tonsa); confirm
###   - sequencing: platform, read length, library prep, provider
###   - read cleaning: tool, version and parameters; QC summary
###   - reference: how the transcriptome was built and filtered ("filteredAssembly" in the
###     counts file name), its size, and how it was annotated (Trinotate/BLASTx/GO sources
###     and versions)
###   - Salmon: version, index settings, mapping rate (per sample, in Salmon's
###     aux_info/meta_info.json -> percent_mapped, or its log)
###   - the blockwiseModules call: power, networkType, TOMType, minModuleSize,
###     mergeCutHeight, maxBlockSize, random seed (the .rds files do not store these)
###   - repository link (GitHub is a rubric item in the example)
### ############################################################################


# ==============================================================================
# 0. SETUP
# PURPOSE : Set the working directory, load packages, create the output folder.
# OUTPUTS : img_dir (character path); packages attached.
# NOTES   : Load order matters. topGO/GO.db (AnnotationDbi) go first and tidyverse
#           last so dplyr's select()/filter()/mutate() are not masked. Do NOT load
#           Rmisc/plyr after tidyverse: plyr::mutate() breaks across().
#           Even so, dplyr:: is used explicitly in later code for safety.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : topGO, GO.db, WGCNA, DESeq2, gridExtra, tidyverse (see the list at the top).
# TRANSFORMATIONS : none. setwd() sets the folder that relative file names resolve against;
#                   dir.create() makes the output folder if it is missing.
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, software
# METHODS : Analyses were run in R [R version] using DESeq2 [version], WGCNA [version], topGO
#           [version], GO.db [version], and tidyverse packages [versions]. Section 15 prints
#           every version.
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 1. LOAD DATA, FILTER, NORMALIZE
# PURPOSE : Read counts and sample metadata, build a DESeq2 object, remove
#           low-expression genes, and produce variance-stabilized expression.
# INPUTS  : salmon.isoform.counts.matrix.filteredAssembly, Ahud_trait_data.txt
# OUTPUTS : countsTable, countsTableRound, sample_metadata, dds,
#           norm.counts (samples x genes, VST; used by sections 3, 9, 13, add-ons)
# NOTES   : Counts are isoform-level (not gene-level), so isoforms of one gene can
#           both appear in a module and inflate module size / enrichment.
#           The two all(...) checks should both print TRUE; if not, samples are
#           missing or ordered differently between counts and metadata.
#           Filter: >= 15 reads in >= 28 samples (~75% of 38 samples).
#           DESeq(dds) is slow and not required for vst(); kept so results match.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : base R (read.table, round); DESeq2 (DESeqDataSetFromMatrix, counts,
#                   DESeq, vst); SummarizedExperiment via DESeq2 (assay).
# TRANSFORMATIONS (in order):
#   read.table(row.names = 1)   first column (gene/isoform ID) becomes row names; columns = samples
#   round()                     decimal Salmon estimates -> integers (DESeq2 requires integers)
#   DESeqDataSetFromMatrix      bundles counts with sample metadata; design = ~ 1 means "no model"
#                               (only normalization is needed here, not differential expression)
#   rowSums(counts(dds) >= 15) >= 28
#                               TRUE/FALSE per gene: ">= 15 reads in >= 28 samples". Failing rows
#                               are dropped (low counts are noisy and distort correlations)
#   DESeq()                     size-factor normalization + dispersion estimates (not needed by vst here)
#   vst()                       variance-stabilizing transformation (blind = TRUE by default, so it
#                               ignores the design). Output is on a roughly log2 scale, so highly
#                               expressed genes no longer dominate the variance
#   t(assay(...))               genes x samples -> samples x genes (WGCNA wants samples in rows)
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, data import and filtering; first sentences of Results (sequencing depth)
# METHODS : Salmon isoform-level counts ([n] transcripts x [38] samples) were rounded to
#           integers and imported into DESeq2 [version]. Transcripts with fewer than 15 reads
#           in at least 28 of the 38 samples (73.7%) were removed, leaving [25,260] transcripts.
#           Counts were variance-stabilized with DESeq2's vst() and transposed to samples x
#           transcripts for WGCNA.
# NOTE    : the code comment says "75%", but 28/38 = 73.7%; quote the exact figure.
# RESULTS : mean counts per sample, counts per transcript, and group sizes come from section 15
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 2. QC: OUTLIER GENES AND SAMPLES
# PURPOSE : Check for bad genes/samples before network construction.
# INPUTS  : countsTable (unfiltered raw counts)
# OUTPUTS : gsg (goodSamplesGenes result); diagnostic plots only, nothing saved.
# NOTES   : Diagnostic section; no object here is used downstream. Look for
#           outlier samples in the hclust tree and the PCA. Replicate-line
#           outliers can distort eigengenes (see the red/brown outliers noted
#           in the trajectory plots, section 12).
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : WGCNA (goodSamplesGenes); stats (hclust, dist, prcomp); ggplot2 via tidyverse.
# TRANSFORMATIONS:
#   t(countsTable)              samples in rows. NOTE: this is the RAW, unfiltered, untransformed
#                               count table, so high-count genes dominate the distances and the PCA.
#                               Fine for spotting gross outliers; results would look different on
#                               norm.counts
#   dist() + hclust(average)    Euclidean distance between samples, then average-linkage clustering tree
#   prcomp()                    principal components (centered, not scaled).
#                               sdev^2 / sum(sdev^2) * 100 = percent of variance per PC (axis labels)
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, QC; Results, sample clustering
# METHODS : Sample- and gene-level quality was checked with WGCNA::goodSamplesGenes,
#           average-linkage hierarchical clustering of Euclidean distances, and PCA (prcomp)
#           of the raw counts. [No samples were removed / state any removed; check gsg$allOK.]
# FIGURE  : the PCA plot (PC1 and PC2 with % variance) can serve as a Figure 1, like the
#           example's PCA. It is not saved: add ggsave(file.path(img_dir, "PCA_samples.png"), ...)
#           and color/shape the points by treatment and generation (join sample_metadata to pca.dat)
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 3. CHOOSING SOFT-THRESHOLD POWERS
# PURPOSE : Justify the powers compared in this project (4, 6, 8) using the
#           scale-free fit (R^2) and mean connectivity at each candidate power.
# INPUTS  : norm.counts
# OUTPUTS : sft, sft.data, plots a1 (R^2 vs power) and a2 (connectivity vs power)
# NOTES   : Slow (a few minutes). Rule of thumb: lowest power with R^2 > 0.8
#           and connectivity that has not collapsed (here 6 to 8). Power 4 is
#           included deliberately as a below-threshold comparison.
#           pickSoftThreshold uses networkType = "signed". The network-building
#           step (section 4 input) must use the SAME networkType, or the chosen
#           powers do not apply to it.
#           cor <- WGCNA::cor masks stats::cor for the rest of the session.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : WGCNA (pickSoftThreshold); ggplot2; gridExtra (grid.arrange).
# TRANSFORMATIONS:
#   c(1:10, seq(12, 50, by = 2))  the grid of candidate powers
#   pickSoftThreshold           for each power: correlations -> network connections (signed:
#                               ((1 + cor)/2)^power) -> each gene's connectivity k -> fit of
#                               log10(p(k)) against log10(k). SFT.R.sq = how straight that line is
#                               (scale-free fit); mean.k. = average connectivity. Higher power
#                               shrinks weak correlations toward 0, so connectivity falls
#   cor <- WGCNA::cor           replaces stats::cor with WGCNA's faster, NA-aware version (session-wide)
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, soft-threshold selection; figure candidate
# METHODS : Signed networks were evaluated over soft-thresholding powers 1-10 and 12-50
#           (step 2) with pickSoftThreshold, using the scale-free topology fit (signed R^2)
#           and mean connectivity. Powers 6-8 [reached R^2 > 0.8; confirm]. Powers 4, 6 and 8
#           were carried forward [state why these three were compared].
# FILL IN : R^2 and mean connectivity at 4, 6, 8 -> methods_info$soft_threshold_fit (section 15)
# FIGURE  : a1/a2 (R^2 and mean connectivity vs power) are drawn but not saved; to save them:
#           ggsave(file.path(img_dir, "soft_power_selection.png"), arrangeGrob(a1, a2, nrow = 2))
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 4. LOAD NETWORKS
# PURPOSE : Load the precomputed WGCNA networks for the three powers and gather
#           them in one named list used by almost every later section.
# INPUTS  : bwnet.4.rds, bwnet.6.rds, bwnet.8.rds (output of blockwiseModules)
# OUTPUTS : bwnet_4, bwnet_6, bwnet_8, nets (list named "4", "6", "8")
# NOTES   : PIPELINE GAP: the blockwiseModules() call that makes these files is
#           not in this script. For reproducibility, record here (or add the
#           call): power, networkType (signed/unsigned), TOMType,
#           minModuleSize, mergeCutHeight, maxBlockSize. The kME boxplot
#           (many genes with negative kME to their own module) suggests the
#           networks may have been built as UNSIGNED, which is the
#           blockwiseModules default; confirm this.
#           Each object has: $colors (named by gene), $MEs (samples x modules,
#           includes MEgrey), $dendrograms, $blockGenes, $unmergedColors.
#           Module COLOR NAMES are reassigned by module size at each power, so
#           only turquoise reliably means the same module across powers.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : base R (readRDS, list). No WGCNA function is called here; the objects were
#                   made by WGCNA::blockwiseModules offline.
# TRANSFORMATIONS : none; loading and grouping only. nets is a named list so later code can loop over
#                   powers by name ("4", "6", "8") instead of repeating code three times.
#                   Contents: $colors = module name per gene, $MEs = module eigengenes.
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, network construction
# METHODS : Co-expression networks were built with WGCNA::blockwiseModules [WGCNA version] at
#           powers 4, 6 and 8 (networkType = [signed/unsigned], TOMType = [ ], minModuleSize =
#           [ ], mergeCutHeight = [ ], maxBlockSize = [ ], [number of blocks]). Modules were
#           labeled by color; genes in no module are assigned to "grey".
# FILL IN : the blockwiseModules arguments (the .rds files do not store them); the number of
#           blocks is in methods_info$network_summary
# NOTE    : state networkType explicitly. pickSoftThreshold above used "signed"
# ----------------------------------------------------------------------------
bwnet_4 <- readRDS("bwnet.4.rds")
bwnet_6 <- readRDS("bwnet.6.rds")
bwnet_8 <- readRDS("bwnet.8.rds")

nets <- list("4" = bwnet_4, "6" = bwnet_6, "8" = bwnet_8)

# ==============================================================================
# 5. DENDROGRAMS
# PURPOSE : Visualize each gene tree with unmerged vs merged module colors.
# INPUTS  : nets
# OUTPUTS : dendrogram_Power_{4,6,8}.png, dendrograms_stacked.png (in img_dir)
# NOTES   : Each power has its own gene tree, so the dendrograms cannot be
#           overlaid. plot_dendro() indexes colors by block (net$blockGenes) so it
#           stays correct if the network was built in multiple blocks.
#           The magick step is optional; stacked images are small, so view the
#           individual PNGs for detail.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : WGCNA (plotDendroAndColors); grDevices (png, dev.off); magick (image_*).
# TRANSFORMATIONS:
#   cbind(unmergedColors[idx], colors[idx])
#                               two color rows drawn under the tree: modules before and after merging
#                               similar modules. [idx] selects the genes of this block
#   png(...) ... dev.off()      sends the plot to a file; dev.off() closes it and writes the file
#   image_append(stack = TRUE)  glues the three PNGs vertically into one image
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Results, figure
# METHODS : Gene dendrograms with module colors before and after merging were drawn with
#           plotDendroAndColors.
# FIGURE  : dendrogram_Power_{4,6,8}.png or dendrograms_stacked.png. The stacked image is hard
#           to read; use one power as a main figure and the others as supplementary figures
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 6. HOW MUCH DO THE NETWORKS DIFFER ACROSS POWERS?
# PURPOSE : Core homework comparison: module counts, gene relabeling, module
#           sizes, and eigengene similarity between powers.
# INPUTS  : nets, bwnet_4/6/8
# OUTPUTS : df (module sizes), p_size, p_cor, cor_df, ME4/ME6/ME8, allMEs,
#           MEcor, results; module_sizes_barplot.png, eigengene_similarity_barplot.png
# NOTES   : "Proportion of genes whose module label differs" compares COLOR
#           NAMES, which are reassigned between powers, so it overstates true
#           gene movement for non-turquoise modules. For a name-independent
#           measure, use the Jaccard / best-match approach (sections 7 and the
#           add-on scripts).
#           Eigengene similarity (cor_df) and the best-match table only make
#           sense for modules matched by correlation; a same-named module is
#           not guaranteed to be the same module.
#           MEcor_x masks same-power pairs so best matches are cross-power only.
#           Turquoise matches across powers at |r| ~ 0.9995+.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : base R (table, factor, cor, outer, apply); WGCNA (cor); dplyr (rename, mutate);
#                   ggplot2.
# TRANSFORMATIONS:
#   table(colors) -> as.data.frame   genes per module; Freq renamed to Count
#   rbind(...)                       stacks the three powers into one long table (module x power rows)
#   factor(levels = ...)             fixes the x-axis order (modules sorted by size at power 6)
#   cbind(setNames(...prefix))       puts all eigengenes in one table with P4_/P6_/P8_ prefixes so
#                                    columns stay distinguishable
#   cor(allMEs)                      eigengene-vs-eigengene correlation matrix
#   abs(), MEcor_x[same_power] <- NA, apply(which.max)
#                                    drop the sign, blank out same-power pairs, then take each row's
#                                    strongest remaining match (so matches are cross-power only)
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, network comparison metrics; Results (the core homework comparison)
# METHODS : Networks were compared across powers by the number of modules, module sizes, the
#           proportion of genes whose module label differed, and Pearson correlations between
#           module eigengenes (WGCNA::cor). A module's best cross-power match was the module
#           whose eigengene correlated most strongly with it at another power.
# RESULTS : numbers to report: labels per power (10, 11, 15; unique labels, including grey if
#           present), proportion relabeled (0.195, 0.281, 0.303), turquoise eigengene |r| ~ 0.9995.
#           The relabel proportion compares color names, which are reassigned by module size, so
#           present it as an upper bound on gene movement (the Jaccard / fate() results in
#           section 8 do not depend on names)
# FIGURE  : module_sizes_barplot.png, eigengene_similarity_barplot.png
# ----------------------------------------------------------------------------
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

get_module_sizes(bwent_4$colors, "4")

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

# ==============================================================================
# 7. GO ENRICHMENT SETUP
# PURPOSE : Build the annotation set and the reusable enrichment functions used
#           by every GO analysis (sections 8, 9, 14 and the add-on scripts).
# INPUTS  : trinotate_annotation_GOblastx_forTopGO.txt
# OUTPUTS : geneID2GO, geneID2GO.filtered, run_go(), jaccard(), compare_module()
# NOTES   : GO terms with < 5 or > 500 annotated genes are dropped to avoid
#           noisy tiny categories and uninformative huge ones. The same
#           filtered annotation is used at every power, so differences between
#           powers come only from module membership.
#           run_go(net, power, module): topGO weight01 + Fisher, ontology BP,
#           returns ALL tested terms (not just the top ones). p-values are NOT
#           corrected across terms; treat as a ranking ("enriched"), not
#           "significant". Very small p-values print as "< 1e-30"; the "<" is stripped.
#           run_go errors if the module has zero genes at that power.
#           compare_module() stops if the module is absent at any power. Same
#           color name != same genes; read $jaccard before comparing GO across
#           powers (turquoise high; yellow near zero).
#           Slow: each compare_module() call runs topGO three times.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : topGO (readMappings, topGOdata, runTest, GenTable, usedGO); dplyr/tidyr
#                   (bind_rows, select, pivot_wider); ggplot2; base R (intersect, union, combn, vapply).
# TRANSFORMATIONS:
#   readMappings                gene ID -> list of GO IDs (topGO's tab-delimited format)
#   table(unlist(geneID2GO))    number of genes directly annotated to each GO term. Terms with < 5 or
#                               > 500 are removed (intersect()), then genes left with no terms are dropped
#   factor(as.integer(colors == module))
#                               0/1 per gene: 1 = in the module, 0 = everything else. topGO tests whether
#                               the 1s are enriched for each GO term
#   sub("::g\\.[0-9]+::m\\.[0-9]+$", "", names(...))
#                               strips the ORF suffix (e.g. ::g.11447::m.11447) so gene IDs match the
#                               names in the mapping file
#   topGOdata                   builds the GO graph and pushes annotations up to parent terms (true-path
#                               rule). weight01 down-weights genes already explained by a more specific
#                               child term, so results are less redundant
#   Fisher test                 2x2 table per term: (in module vs not) x (annotated vs not)
#   as.numeric(gsub("<", "", weightFisher))
#                               GenTable returns p-values as text (e.g. "< 1e-30"); this makes them
#                               numeric (the bound itself becomes the value)
#   pivot_wider                 one row per GO term, one p-value column per power
#   pmin(..., na.rm = TRUE)     best (smallest) p-value across powers, used for ranking
#   -log10(pmax(p, 1e-30))      plotting scale; floored at 1e-30 so tiny p-values do not stretch the axis
#   jaccard()                   |intersection| / |union| of two gene sets
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, functional enrichment
# METHODS : GO enrichment of module genes used topGO [version] (ontology: biological process;
#           algorithm: weight01; statistic: Fisher's exact test). Annotation came from
#           [Trinotate/BLASTx; version] for [n] transcripts. GO terms with fewer than 5 or more
#           than 500 annotated transcripts were removed ([n] terms kept; [n] transcripts retained
#           at least one term). For each module, its transcripts were tested against all annotated
#           transcripts. P-values were not corrected across terms.
# FILL IN : the counts -> methods_info$go (section 15); the annotation source and version
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 8. GO ACROSS POWERS: TURQUOISE, YELLOW, BLUE, BROWN
# PURPOSE : Compare the functional enrichment of selected modules at each power.
# INPUTS  : compare_module() and nets from earlier sections
# OUTPUTS : res_turq, res_yel, res_blue, res_brwn (each a list: $go, $sizes,
#           $jaccard, $term_overlap, $wide, $sig_one, $plot, $top5);
#           GO_comparison_turquoise.png, GO_comparison_yellow.png,
#           top5_GO_terms_by_module_power.csv, fate()
# NOTES   : Read $jaccard first. Turquoise is the same module at all powers
#           (stable; muscle/flight enrichment). Yellow, blue, and brown share
#           little gene content across powers, so a GO comparison by color name
#           mostly compares different gene sets.
#           blue and brown plots are not saved, and top5_all includes only
#           turquoise and yellow; add res_blue/res_brwn to rbind() if needed.
#           fate(from, to, module): where the genes of one module (at one power)
#           land at another power. Use it to find the real counterpart of a module.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : base R (table, sort, write.csv); ggplot2 (plots returned by compare_module).
# TRANSFORMATIONS:
#   compare_module()            returns a list of tables and a plot (see section 7)
#   rbind(res_turq$top5, ...)   stacks the top-5 term tables into one
#   fate(): to$colors[ids]      looks up, for the genes of one module at power A, their module label at
#                               power B; table() counts them and sort() orders them by size
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods (module overlap); Results; table and figures
# METHODS : Gene-set overlap between modules at different powers was measured with the Jaccard
#           index. GO results were compared across powers by the overlap of the top 50 terms and
#           by -log10(p) of the top 15 terms.
# RESULTS : report the Jaccard values per module (res_*$jaccard). Yellow was 0.16 / 0 / 0.002
#           for 4v6 / 4v8 / 6v8, i.e. "yellow" is a different gene set at each power, while
#           turquoise is the module preserved across powers
# TABLE   : top5_GO_terms_by_module_power.csv is ready to become a results table
# FIGURE  : GO_comparison_turquoise.png, GO_comparison_yellow.png
# ----------------------------------------------------------------------------
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

# Run after section 8 of the script (needs top5_all and img_dir).
library(gridExtra)   # tableGrob, arrangeGrob, ttheme_default
library(grid)        # textGrob, gpar, grid.draw

# Top 5 GO terms per module and power (same table as before, plus the CSV)
top5_all <- rbind(res_turq$top5, res_yel$top5)
write.csv(top5_all, file.path(img_dir, "top5_GO_terms_by_module_power.csv"), row.names = FALSE)


# Display table with report-ready column names
tab <- data.frame(
  Module = top5_all$Module,
  Power  = top5_all$Power,
  Rank   = top5_all$Rank,
  `GO term (biological process)` = top5_all$Term,
  `Genes (module / annotated)`   = paste0(top5_all$Significant, " / ", top5_all$Annotated),
  `p (weight01)`                 = top5_all$weightFisher,
  check.names = FALSE
)

# Each row is lightly shaded in its module's color; text is left-aligned
theme_go <- ttheme_default(
  base_size = 11,
  core    = list(bg_params = list(fill = adjustcolor(tab$Module, alpha.f = 0.25), col = "white"),
                 fg_params = list(hjust = 0, x = 0.05)),
  colhead = list(bg_params = list(fill = "grey85", col = "white"),
                 fg_params = list(hjust = 0, x = 0.05, fontface = "bold"))
)

g <- tableGrob(tab, rows = NULL, theme = theme_go)

out <- arrangeGrob(
  g,
  top    = textGrob("Top 5 enriched GO terms per module and soft power",
                    gp = gpar(fontsize = 14, fontface = "bold"), x = 0, hjust = 0),
  bottom = textGrob("topGO, biological process, weight01 + Fisher's exact test; p-values are not corrected across terms.",
                    gp = gpar(fontsize = 9, fontface = "italic"), x = 0, hjust = 0)
)

# Height grows with the number of rows
png(file.path(img_dir, "top5_GO_terms_table.png"),
    width = 11, height = 0.32 * nrow(tab) + 1.6, units = "in", res = 250)
grid.newpage()
grid.draw(out)
dev.off()

# Where did the yellow genes from one power end up at another power?
fate <- function(from, to, module) {
  ids <- names(from$colors)[from$colors == module]
  sort(table(to$colors[ids]), decreasing = TRUE)
}
fate(bwnet_6, bwnet_4, "yellow")   # power 6 yellow -> power 4 modules
fate(bwnet_6, bwnet_8, "yellow")   # power 6 yellow -> power 8 modules
fate(bwnet_8, bwnet_6, "yellow")   # power 8 yellow -> power 6 modules

# ==============================================================================
# 9. HUB GENES (TURQUOISE, POWER 6)
# PURPOSE : Identify the most connected genes in a module (highest |kME|) and
#           look up their GO annotations.
# INPUTS  : norm.counts, net$MEs, geneID2GO
# OUTPUTS : module_hubs() (reused in section 14 and the add-ons), top5 (hubs),
#           top5_go (GO terms per hub gene, with kME)
# NOTES   : kME = correlation of a gene's expression with the module eigengene.
#           Gene order of norm.counts must match names(net$colors); the
#           stopifnot() checks this.
#           The ::g.#::m.# suffix is removed so IDs match the GO mapping file.
#           A hub with no GO annotation (e.g. DN12290) returns NA. For protein
#           names, use the Trinotate report (sprot_Top_BLASTX_hit, Pfam).
#           The top two hubs are isoforms of one gene (identical kME).
#           NA terms in top5_go = GO IDs obsolete or missing in GO.db.
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : WGCNA (signedKME); GO.db (Term, Ontology); base R.
# TRANSFORMATIONS:
#   signedKME(norm.counts, MEs) correlation of every gene with every eigengene
#                               (genes x modules; columns are named kME<color>)
#   order(-abs(kME))            ranks genes by strength of membership, regardless of sign
#   sub("::g...")               strips the ORF suffix so IDs match the geneID2GO names
#   geneID2GO[[id]] -> Term(go) GO IDs -> readable names (NA if the ID is obsolete in GO.db)
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods (hub genes); Results table
# METHODS : Module membership (kME) was the correlation of each transcript's variance-stabilized
#           expression with its module eigengene (WGCNA::signedKME). Hub genes were the
#           transcripts with the largest |kME| in the module; their GO annotations were looked up
#           in the same Trinotate-derived mapping.
# FILL IN : for a results table: transcript ID, kME, annotation (add the Trinotate BLAST hit
#           for protein names, since GO terms alone are uninformative here)
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 10. TRAIT TABLE
# PURPOSE : Turn the metadata into a numeric table that can be correlated with
#           module eigengenes.
# INPUTS  : sample_metadata
# OUTPUTS : traits (samples x: Fitness, EPR, HatchSucc, Survival, DevTime,
#           Generation, trt_AM, trt_OA, trt_OW, trt_OWA)
# NOTES   : Uses the mean columns only (the CI columns, Rep, and the duplicate
#           generation label are dropped).
#           Design is unbalanced: F11 exists only for AM and OWA, AM F2 has 2
#           replicates, and every cell has about 3. Traits are collinear with each
#           other and with treatment (development time is essentially a warming
#           indicator), so module-trait hits cannot say which factor drives them.
#           Confirm treatment labels (AM/OA/OW/OWA) against the methods.
# ==============================================================================
# (means only; the CI columns, Rep, and the duplicate generation label are dropped)
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : base R only (data.frame, model.matrix, cbind, table, cor).
# TRANSFORMATIONS:
#   data.frame(...)             keeps only the mean trait columns, renamed to short labels
#   model.matrix(~ 0 + treatment)
#                               "dummy coding": one 0/1 column per treatment (no intercept, so all four
#                               treatments get a column). Lets a categorical variable be correlated
#   cbind                       appends the dummy columns to the trait table
#   Generation is kept numeric (0, 2, 4, 11) here, which assumes a linear trend across generations
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, phenotype data
# METHODS : Per-sample phenotypes (fitness, egg production rate, hatching success, survival,
#           development time; means across [what was averaged]) and experimental factors
#           (treatment, generation) came from Ahud_trait_data.txt. Treatment was dummy coded;
#           generation was numeric (0, 2, 4, 11).
# FILL IN : how each phenotype was measured, and the definition of each treatment
# NOTE    : the design is unbalanced (F11 only in AM and OWA; AM F2 has 2 replicates). Report
#           the sample size per group using the table() output in this section
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 11. MODULE-TRAIT CORRELATIONS AT EACH POWER
# PURPOSE : Correlate every module eigengene (excluding grey) with every trait.
# INPUTS  : nets, traits
# OUTPUTS : mt_all (long table: Module, Trait, r, p, Power, padj, label),
#           plot_mt(), p_mt4/p_mt6/p_mt8; also the turquoise-across-powers table
# NOTES   : ~100 tests per power with only ~38 samples. BH-adjusted p < 0.05 is
#           marked "**", raw p < 0.05 is marked "*". Treat single-star hits as
#           suggestive. Replicate lines within a group are not independent
#           observations of a treatment effect, so p-values are somewhat optimistic.
#           Heatmaps are shown but not saved; add ggsave() if needed.
#           Only turquoise is the same module at every power; for other modules
#           use the gene-matched comparison (trait_correlations_across_powers.R).
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : WGCNA (cor, corPvalueStudent); dplyr/tidyr (filter, select, pivot_wider); ggplot2.
# TRANSFORMATIONS:
#   net$MEs[, != "MEgrey"]      drops the unassigned-gene "module"
#   traits[rownames(MEs), ]     reorders the trait rows to match the eigengene rows by sample name
#   cor(MEs, traits)            module x trait correlation matrix (Pearson)
#   corPvalueStudent            turns each r into a two-sided p-value through the t statistic
#                               t = r * sqrt((n - 2) / (1 - r^2)), n = number of samples
#   as.table -> as.data.frame   wide matrix -> long table (one row per module x trait)
#   ave(p, Power, p.adjust BH)  Benjamini-Hochberg multiple-testing correction within each power
#   sprintf/paste0/ifelse       builds the cell labels: r to 2 decimals plus * or **
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, module-trait analysis; Results figure
# METHODS : Module eigengenes (excluding grey) were correlated with each trait (Pearson).
#           P-values came from the Student t distribution (corPvalueStudent) and were
#           Benjamini-Hochberg adjusted within each power ([n] tests per power).
# FIGURE  : the p_mt4/p_mt6/p_mt8 heatmaps are drawn but not saved; to save one:
#           ggsave(file.path(img_dir, "module_trait_heatmap_P6.png"), p_mt6, width = 10, height = 6)
#           Cells marked * (raw p < 0.05) and ** (adjusted p < 0.05); report both honestly
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 12. EIGENGENES BY TREATMENT AND GENERATION (POWER 6, THEN 4, THEN 8)
# PURPOSE : Test whether module eigengenes differ by treatment or generation and
#           plot their trajectories across generations, one line per treatment.
# INPUTS  : bwnet_6/4/8 $MEs, sample_metadata
# OUTPUTS : MEs6/MEs4/MEs8, meta, mod_tests, me_long, me_sum, p_traj,
#           eigengene_trajectories.png
# NOTES   : Three copies of the same block, one per power (power 6 first).
#           PIPELINE WARNING: the blocks reuse the same object names (meta,
#           mod_tests, me_long, me_sum, p_traj) so each block OVERWRITES the last,
#           leaving the power-8 results in memory. The power-6 and power-4 blocks
#           both ggsave() to the same file (eigengene_trajectories.png), so power 4
#           overwrites power 6 on disk and the power-8 plot is never saved.
#           Give each power its own file name, or wrap the block in a function of
#           (net, power) (offer: I can write it).
#           mod_tests: per-module F-tests for treatment and generation, each
#           adjusted for the other with drop1(). Exploratory only: the design is
#           unbalanced, and a treatment x generation interaction is not tested.
#           Check single-replicate outliers (e.g. red OWA F4, brown OA F2).
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : stats (lm, drop1, p.adjust); dplyr (filter, mutate, group_by, summarise);
#                   tidyr (pivot_longer); ggplot2.
# TRANSFORMATIONS:
#   sample_metadata[rownames(MEs), ]
#                               puts the metadata in the same sample order as the eigengenes
#   factor(treatment), factor(generation)
#                               generation is treated as categories (0, 2, 4, 11) here, not a trend,
#                               which differs from section 10
#   lm(ME ~ treatment + gen) + drop1(test = "F")
#                               additive model; drop1 gives an F-test for each term with the other
#                               term already in the model
#   pivot_longer(starts_with("ME"))
#                               wide (one column per module) -> long (one row per sample x module),
#                               the shape ggplot/facet_wrap needs
#   sub("^ME", "")              removes the "ME" prefix so labels are plain color names
#   group_by + summarise        mean and standard error (sd / sqrt(n)) per module x treatment x generation
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, treatment x generation analysis; Results figure
# METHODS : For each module, the eigengene was modeled as eigengene ~ treatment + generation
#           (both as factors) with linear models, and each term was tested with marginal F-tests
#           (drop1), BH-adjusted across modules. Mean +/- SE trajectories were plotted per
#           treatment. The design is unbalanced, so tests are exploratory, and no treatment x
#           generation interaction was tested.
# FIGURE  : eigengene_trajectories.png (see the overwrite warning in the notes above: save one
#           file per power)
# ----------------------------------------------------------------------------
# --- Power 6 -------------------------------------------------------------------
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
# --- Power 4 (same analysis; overwrites meta, mod_tests, me_long, me_sum, p_traj
#     and the saved eigengene_trajectories.png; see NOTE above) -----------------
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

# --- Power 8 (same analysis; overwrites the objects again; plot is NOT saved) ---
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



# ==============================================================================
# 13. CANDIDATE GENES: MODULE MEMBERSHIP VS GENE-TRAIT CORRELATION
# PURPOSE : For one module and one trait, find genes that are both well connected
#           to the module (kME) and strongly correlated with the trait (GS).
# INPUTS  : norm.counts, bwnet_6, traits; edit focus_mod and focus_trait
# OUTPUTS : gene_df, p_mmgs, MM_vs_GS.png, a short list of candidate genes
# NOTES   : Uses power 6 only. A positive kME-vs-|GS| relationship means the
#           module's hub genes track the trait most, which supports the module as
#           trait-relevant.
#           Candidates need functional annotation (Trinotate BLAST hits) before
#           any biological claim; GO terms alone are limited here.
#           Fitness is strongly tied to treatment/generation in this design.
# ==============================================================================
# 4. Candidate genes: module membership vs gene-trait correlation ----------------
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : WGCNA (signedKME, cor); ggplot2 (geom_smooth with lm); stats (cor.test).
# TRANSFORMATIONS:
#   kME                         gene vs module-eigengene correlation (membership)
#   GS = cor(norm.counts, trait)
#                               gene vs trait correlation across samples (one value per gene)
#   abs(GS), -abs(kME) * abs(GS)
#                               ranks genes that are strong on both; sign is ignored so up- and
#                               down-regulated genes both count
#   cor.test(kME, abs(GS))      do hub genes (high kME) also track the trait more strongly?
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, candidate genes; Results
# METHODS : Gene significance was the Pearson correlation of each transcript's expression with
#           [fitness]; kME and |GS| were compared within [turquoise] (cor.test), and transcripts
#           high on both were listed as candidates.
# ----------------------------------------------------------------------------
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

# ==============================================================================
# 14. YELLOW (POWER 6) FOLLOW-UP: WARMING vs NON-WARMING
# PURPOSE : Test whether the power-6 yellow eigengene separates warming (OW, OWA)
#           from non-warming (AM, OA) samples at F0 and at F4, then characterize
#           the module (GO terms, hub genes).
# INPUTS  : bwnet_6$MEs, meta, run_go(), module_hubs()
# OUTPUTS : Welch t-tests printed for F0 and F4; top GO terms; top 10 hubs
# NOTES   : Observed: strong separation at F0 (p ~ 5e-5), none at F4 (p ~ 0.79),
#           consistent with a transient early response. Not proof of absence at F4
#           (wide CI), and "plastic vs adaptive" cannot be separated with this data.
#           This module was picked after looking at the data. Check that it is
#           special by running the same test on every module with BH adjustment, and
#           test the warm x generation interaction formally.
#           Power-6 yellow is NOT the power-8 yellow (Jaccard ~ 0.002), so do not
#           link their GO results. Top terms were cell-cycle/mitosis, with muscle
#           regulation, salinity, and collagen terms based on few genes.
#           Uses `meta`, which holds the sample metadata in the same order for all
#           powers (the last block of section 12 set it from MEs8; row order is identical).
# ==============================================================================
# ---- Packages and transformations in this section ---------------------------
# PACKAGES        : stats (t.test); topGO and WGCNA through run_go() and module_hubs().
# TRANSFORMATIONS:
#   meta$Generation == g        selects the samples of one generation
#   %in% c("OW", "OWA")         turns treatment labels into a TRUE/FALSE "warming" grouping
#   t.test(ME ~ warm)           Welch two-sample t-test (unequal variances) comparing the yellow
#                               eigengene between the two groups (6 vs 6 samples at F0 and at F4)
# ----------------------------------------------------------------------------
# ---- REPORT NOTES (feeds the Bioinformatics Pipeline section) -----------------
# REPORT  : Methods, planned contrast; Results
# METHODS : The yellow (power 6) eigengene was compared between warming (OW, OWA) and
#           non-warming (AM, OA) lines at F0 and at F4 with Welch two-sample t-tests (6 vs 6
#           samples per generation).
# NOTE    : state that yellow was selected after exploring the data (post hoc), and say whether
#           the all-modules test and the warm x generation interaction were run
# ----------------------------------------------------------------------------
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


# ==============================================================================
# 15. METHODS SUMMARY
# PURPOSE : Collect every number and software version the Bioinformatics Pipeline section
#           needs in one place, so the report is not written from memory.
# INPUTS  : objects from earlier sections (countsTable, dds, sft.data, nets, geneID2GO,
#           geneID2GO.filtered, keepTerms, traits, sample_metadata); img_dir
# OUTPUTS : methods_info (list); methods_summary.txt, sessionInfo.txt,
#           network_summary_by_power.csv (in img_dir)
# NOTES   : Run last. It reads objects only and changes no analysis result.
#           Hard-coded filter values (15 and 28) must match section 1; the TRUE/FALSE check
#           below confirms the filtered data satisfy them.
#           Things R cannot supply (Salmon mapping rate, fastp settings, blockwiseModules
#           arguments, annotation source) are listed in the REPORT GUIDE at the top.
# ==============================================================================
methods_info <- list()

# Software ----------------------------------------------------------------------
methods_info$R_version <- R.version.string
pkgs <- c("DESeq2", "WGCNA", "topGO", "GO.db", "ggplot2", "dplyr", "tidyr",
          "gridExtra", "magick")
methods_info$package_versions <- sapply(pkgs, function(p) as.character(packageVersion(p)))

# Data and filtering --------------------------------------------------------------
methods_info$n_samples                <- ncol(countsTable)
methods_info$n_transcripts_raw        <- nrow(countsTable)
methods_info$n_transcripts_after_filter <- nrow(dds)
methods_info$filter_rule              <- ">= 15 reads in >= 28 samples"
methods_info$filter_samples_pct       <- round(100 * 28 / ncol(countsTable), 1)
methods_info$filter_rule_satisfied    <- all(rowSums(counts(dds) >= 15) >= 28)
methods_info$total_counts_per_sample  <- round(summary(colSums(countsTableRound)))
methods_info$total_counts_per_transcript <- round(summary(rowSums(countsTableRound)))
methods_info$samples_per_group        <- table(Treatment = sample_metadata$treatment,
                                               Generation = sample_metadata$Generation)

# Soft threshold ------------------------------------------------------------------
methods_info$soft_threshold_fit <- sft.data[sft.data$Power %in% c(4, 6, 8),
                                            c("Power", "SFT.R.sq", "mean.k.")]

# Networks ------------------------------------------------------------------------
methods_info$network_summary <- do.call(rbind, lapply(names(nets), function(p) {
  cl <- nets[[p]]$colors
  sz <- table(cl[cl != "grey"])
  data.frame(Power                = p,
             n_genes              = length(cl),
             n_blocks             = length(nets[[p]]$dendrograms),
             n_modules_excl_grey  = length(sz),
             n_grey               = sum(cl == "grey"),
             pct_grey             = round(100 * mean(cl == "grey"), 1),
             largest_module       = names(sz)[which.max(sz)],
             largest_module_size  = max(sz),
             smallest_module_size = min(sz))
}))

# Proportion of genes with a different module LABEL (color names are reassigned by size)
methods_info$prop_genes_relabeled <- c(
  "4_vs_6" = mean(bwnet_4$colors != bwnet_6$colors),
  "4_vs_8" = mean(bwnet_4$colors != bwnet_8$colors),
  "6_vs_8" = mean(bwnet_6$colors != bwnet_8$colors))

# GO enrichment -------------------------------------------------------------------
methods_info$go <- list(ontology                  = "BP",
                        algorithm                 = "weight01",
                        statistic                 = "Fisher",
                        term_size_filter          = "5-500 annotated genes",
                        n_genes_in_mapping        = length(geneID2GO),
                        n_genes_after_term_filter = length(geneID2GO.filtered),
                        n_terms_kept              = length(keepTerms))

# Traits and statistical tests ----------------------------------------------------
methods_info$traits_used <- colnames(traits)
methods_info$n_module_trait_tests_per_power <- (ncol(bwnet_6$MEs) - 1) * ncol(traits)  # grey excluded; power 6 shown
methods_info$tests <- c(
  module_trait    = "Pearson correlation; Student-t p-values (corPvalueStudent); BH within power",
  eigengene_model = "lm(ME ~ treatment + generation), factors; marginal F-tests (drop1); BH across modules",
  yellow_contrast = "Welch two-sample t-test, warming (OW, OWA) vs non-warming (AM, OA), at F0 and F4",
  go              = "topGO weight01 + Fisher's exact test; no correction across terms")

# Write out -----------------------------------------------------------------------
writeLines(capture.output(print(methods_info)), file.path(img_dir, "methods_summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(img_dir, "sessionInfo.txt"))
write.csv(methods_info$network_summary, file.path(img_dir, "network_summary_by_power.csv"),
          row.names = FALSE)

methods_info