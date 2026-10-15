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
key_mods <- c("turquoise", "yellow", "red", "brown")   # edit to whichever modules stand out in step 2

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




### Fate of the power-6 turquoise module at power 8 ---------------------------

min_genes <- 30   # a power-8 module counts as a "child" if it holds >= this many turquoise-6 genes

# 1. Where did the power-6 turquoise genes go at power 8? ------------------------
tq6      <- names(bwnet_6$colors)[bwnet_6$colors == "turquoise"]
child_of <- bwnet_8$colors[tq6]
length(tq6)

fate_tab <- as.data.frame(table(Child = child_of), stringsAsFactors = FALSE)
names(fate_tab)[names(fate_tab) == "Freq"] <- "n_from_tq6"

sizes8 <- table(bwnet_8$colors)
fate_tab$pct_of_tq6           <- round(100 * fate_tab$n_from_tq6 / length(tq6), 1)
fate_tab$size_P8              <- as.integer(sizes8[fate_tab$Child])
fate_tab$pct_of_child_from_tq6 <- round(100 * fate_tab$n_from_tq6 / fate_tab$size_P8, 1)
fate_tab <- fate_tab[order(-fate_tab$n_from_tq6), ]
fate_tab
# pct_of_tq6:             how the old module was divided
# pct_of_child_from_tq6:  how "pure" each child is (low = it also absorbed non-turquoise genes)

children <- fate_tab$Child[fate_tab$n_from_tq6 >= min_genes & fate_tab$Child != "grey"]
children

p_fate <- ggplot(fate_tab, aes(x = reorder(Child, -n_from_tq6), y = n_from_tq6,
                               fill = Child)) +
  geom_col() +
  scale_fill_identity() +
  theme_bw() +
  labs(title = "Fate of power-6 turquoise genes at power 8",
       x = "Power-8 module", y = "Number of turquoise-6 genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_fate
ggsave(file.path(img_dir, "turquoise_fate_P6_to_P8.png"), p_fate,
       width = 8, height = 5, dpi = 250)

# 2. Do the children have distinct expression patterns? --------------------------
# (eigengene signs are arbitrary, so focus on |r|)
ME_cmp <- cbind(P6_turquoise = bwnet_6$MEs$MEturquoise,
                setNames(bwnet_8$MEs[, paste0("ME", children), drop = FALSE],
                         paste0("P8_", children)))
cm      <- WGCNA::cor(ME_cmp)
cm_long <- as.data.frame(as.table(cm))

p_cm <- ggplot(cm_long, aes(Var1, Var2, fill = Freq)) +
  geom_tile() +
  geom_text(aes(label = sprintf("%.2f", Freq)), size = 3) +
  scale_fill_gradient2(low = "navy", mid = "white", high = "firebrick3", limits = c(-1, 1)) +
  theme_minimal() +
  labs(title = "Eigengene correlations: turquoise (P6) and its power-8 children",
       x = NULL, y = NULL, fill = "r") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_cm
ggsave(file.path(img_dir, "turquoise_children_eigengene_cor.png"), p_cm,
       width = 8, height = 6, dpi = 250)
# Children that correlate weakly with each other have different expression patterns.

# 3. GO enrichment of each child ---------------------------------------------------
# A) Child vs the REST of turquoise-6 (background = turquoise-6 genes only).
#    This asks what distinguishes each child from its siblings.
run_go_subset <- function(selected, universe, label, nodeSize = 5) {
  geneList <- factor(as.integer(universe %in% selected))
  names(geneList) <- sub("::g\\.[0-9]+::m\\.[0-9]+$", "", universe)
  
  GOdata <- new("topGOdata",
                ontology = "BP",
                allGenes = geneList,
                geneSelectionFun = function(x) x == 1,
                annot = annFUN.gene2GO,
                gene2GO = geneID2GO.filtered,
                nodeSize = nodeSize)
  
  res <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  tab <- GenTable(GOdata, weightFisher = res, orderBy = "weightFisher",
                  topNodes = length(usedGO(GOdata)))
  tab$p     <- as.numeric(gsub("<", "", tab$weightFisher))
  tab$Child <- label
  tab
}

go_within <- setNames(lapply(children, function(ch) {
  genes_ch <- tq6[child_of == ch]
  run_go_subset(genes_ch, universe = tq6, label = ch)
}), children)

# B) Whole power-8 child module vs all genes (includes any non-turquoise genes it absorbed)
go_whole <- setNames(lapply(children, function(ch) {
  tab <- run_go(bwnet_8, "8", module = ch)
  tab$Child <- ch
  tab
}), children)

top_n <- function(go_list, n = 5) {
  out <- do.call(rbind, lapply(go_list, function(tab) {
    x <- head(tab[, c("Child", "GO.ID", "Term", "Annotated", "Significant", "weightFisher")], n)
    x$Rank <- seq_len(nrow(x))
    x
  }))
  rownames(out) <- NULL
  out
}

top_within <- top_n(go_within)
top_whole  <- top_n(go_whole)

top_within   # what distinguishes each child from the rest of turquoise
top_whole    # what each whole power-8 child is enriched for

write.csv(top_within, file.path(img_dir, "turquoise_children_GO_within.csv"), row.names = FALSE)
write.csv(top_whole,  file.path(img_dir, "turquoise_children_GO_whole.csv"),  row.names = FALSE)

# 4. Plot top terms per child (child vs rest of turquoise) ---------------------------
plot_df <- do.call(rbind, lapply(go_within, function(tab) head(tab, 5)))
plot_df$neglogp <- -log10(pmax(plot_df$p, 1e-30))
plot_df$key     <- paste(plot_df$Child, plot_df$Term, sep = "___")
plot_df$key     <- factor(plot_df$key, levels = plot_df$key[order(plot_df$Child, plot_df$neglogp)])

p_go <- ggplot(plot_df, aes(x = neglogp, y = key, fill = Child)) +
  geom_col() +
  scale_fill_identity() +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  facet_wrap(~ Child, ncol = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Top GO terms for each power-8 child (vs rest of turquoise-6)",
       x = "-log10(weightFisher p-value)", y = NULL)
p_go
ggsave(file.path(img_dir, "turquoise_children_GO.png"), p_go,
       width = 10, height = 3 * length(children), dpi = 200, limitsize = FALSE)
