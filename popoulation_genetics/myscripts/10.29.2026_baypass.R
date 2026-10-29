library(tidyverse)

OUTDIR="~/projects/eco_genomics_2026/popoulation_genetics/myresults/baypass"

setwd("~/projects/eco_genomics_2026/popoulation_genetics/mydata")

pops <- readLines("baypass/Lottia.baypass.pops.txt")
npop <- length(pops)

env <- read_csv("pop_level_biooracle.csv")
pop_order<-env$pop[order(env$lat)]

omega <- as.matrix(read.table("baypass/Lottia_bp_mat_omega.out"))
dimnames(omega) <- list(pops, pops)

omega_cor <- cov2cor(omega)
round(omega_cor, 2)

omega_tree <- hclust(as.dist(1 - omega_cor), method = "average")

### we are not using ggplot so it is saved in a different way:
png(paste0(OUTDIR,"/baypass_omega_tree.png"), width = 1800, height = 1200, res = 300)
plot(omega_tree, main = "Populations clustered by BayPass Omega", xlab = "", sub = "")
dev.off()

##Omega Matrix

#Outliers to table
omega_df <- as.data.frame(as.table(omega_cor)) %>%
  rename(pop1 = Var1, pop2 = Var2, correlation = Freq) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))

###if you get any errors for visuals, try dev.off()
ggplot(omega_df, aes(x = pop1, y = pop2, fill = correlation)) +
  geom_tile(colour = "white") +
  scale_fill_gradient2(low = "cadetblue", mid = "white", high = "red4",
                       midpoint = 0, limits = c(-1, 1)) +
  labs(x = NULL, y = NULL,
       title = "BayPass Omega as correlations (south to north)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUTDIR,"/baypass_omega.png"), width = 7, height = 6)

##XtX*: SNPs that are more differentiated than expected
snps <- read_tsv("baypass/Lottia.baypass.snps.txt")

snps$MRK<-rownames(snps)
snps$MRK<-as.numeric(snps$MRK)

npop <- 19

xtx <- read.table("baypass/Lottia_bp_summary_pi_xtx.out", header = TRUE, check.names = FALSE) %>%
  left_join(snps, by = "MRK") %>%
  mutate(lp_high   = pchisq(XtXst, npop, lower.tail = FALSE, log.p = TRUE),   #
         lp_low    = pchisq(XtXst, npop, lower.tail = TRUE,  log.p = TRUE),   
         direction = if_else(lp_high < lp_low, "high", "low"),
         log10_p   = -(pmin(lp_high, lp_low) + log(2)) / log(10),             
         pval      = 10^(-log10_p),
         qval      = p.adjust(pval, method = "BH"))

#Before trusting the p-values, look at their distribution:
# 
ggplot(xtx, aes(x = pval)) +
  geom_histogram(breaks = seq(0, 1, by = 0.025), fill = "grey60", colour = "white") +
  labs(x = "p-value", y = "no. of SNPs") +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_xtx_pvalues.png"), width = 5, height = 4)

##
##same helper tool as earlier with FST
add_genome_pos <- function(df) {
  scaffolds <- df %>% 
    group_by(CHROM) %>% ##group by CHROM
    summarise(length = max(POS)) %>% ## length = max position within a chromosome
    arrange(match(CHROM, unique(df$CHROM))) %>% ##make sure each is unique, arranged
    mutate(offset = cumsum(as.numeric(length)) - length, ###cumulatively number (helps with labels)
           shade  = factor(row_number() %% 2))    ##value for alternating colour 
  df %>%
    left_join(scaffolds, by = "CHROM") %>% ##join this data by given data with CHROM
    mutate(genome_pos = offset + POS) ###cumalative genome pos 
}


xtx <- xtx %>%
  add_genome_pos()     ##recall above if reloaded R

sum(xtx$outlier)                           

threshold <- quantile(xtx$qval, 0.001)   # top 0.1% (smallest q)

ggplot(xtx, aes(genome_pos / 1e6, -log10(qval))) +
  geom_point(aes(colour = shade), size = 0.4) +
  geom_hline(yintercept = -log10(threshold), linetype = "dashed", colour = "red") +
  scale_colour_manual(values = c("grey20", "grey80"), guide = "none") +
  labs(x = "Genome pos (Mb)", y = "-log10(q)") +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_xtx_scan.png"), width = 10, height = 5)

##Genotype–environment association

covariates <- readLines("baypass/Lottia.baypass.covariates.txt")

bf <- read.table("baypass/Lottia_bp_summary_betai_reg.out", header = TRUE, check.names = FALSE) %>%
  rename(BF_dB = `BF(dB)`) %>%
  mutate(covariate = covariates[COVARIABLE]) %>%
  left_join(snps, by = "MRK") %>%
  add_genome_pos()

bf %>%
  group_by(covariate) %>%
  summarise(n_decisive = sum(BF_dB > 20))

ggplot(bf, aes(x = genome_pos / 1e6, y = BF_dB)) +
  geom_point(aes(colour = shade), size = 0.6) +
  geom_point(data = filter(bf, BF_dB > 20), colour = "red", size = 1) +
  geom_hline(yintercept = 20, linetype = "dashed") +
  facet_wrap(~ covariate, ncol = 1) +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Genome pos (Mb)",
       y = "Bayes factor (dB)") +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_bf_scan.png"), width = 10, height = 6)

#What does an association actually look like? Let’s plot the allele frequencies of the SNP with the strongest association with sea surface temperature:

geno <- as.matrix(read.table("baypass/Lottia.baypass.geno"))

top_snp <- bf %>%
  filter(covariate == "thetao_mean") %>%
  slice_max(BF_dB, n = 1, with_ties = FALSE)
top_snp

counts <- geno[top_snp$MRK, ]
ref <- counts[seq(1, length(counts), by = 2)]
alt <- counts[seq(2, length(counts), by = 2)]

top_freq <- tibble(pop = pops, alt_freq = alt / (ref + alt)) %>%
  left_join(env, by = "pop")

ggplot(top_freq, aes(x = thetao_mean, y = alt_freq, label = pop)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  ylim(0, 1) +
  labs(x = "Mean sea surface temperature (C)", y = "Allele Frequency",
       title = paste0("Top SST-associated SNP: ", top_snp$CHROM, ":", top_snp$POS)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_top_sst_snp.png"), width = 5, height = 4)

##lets repeat with lat

top_snp <- bf %>%
  filter(covariate == "lat") %>%
  slice_max(BF_dB, n = 1, with_ties = FALSE)
top_snp

counts <- geno[top_snp$MRK, ]
ref <- counts[seq(1, length(counts), by = 2)]
alt <- counts[seq(2, length(counts), by = 2)]

top_freq <- tibble(pop = pops, alt_freq = alt / (ref + alt)) %>%
  left_join(env, by = "pop")

ggplot(top_freq, aes(x = lat, y = alt_freq, label = pop)) +
  geom_point(size = 3) +
  geom_text(vjust = -1, size = 3) +
  ylim(0, 1) +
  labs(x = "Latitude", y = "Allele Frequency",
       title = paste0("Top Latitude-associated SNP: ", top_snp$CHROM, ":", top_snp$POS)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/baypass_top_lat_snp.png"), width = 5, height = 4)