library(tidyverse)

setwd("~/projects/eco_genomics_2026/popoulation_genetics/mydata")
###point to your results folder
OUTDIR="~/projects/eco_genomics_2026/popoulation_genetics/myresults"

env <- read_csv("sample_metadata.csv")
pw  <- read_tsv("fst/pairwise_fst_summary.txt")

# order the populations from south to north
env <- env %>%
  distinct(Population, Latitude, Longitude)
pop_order <- env$Population[order(env$Latitude)]

# add each pair in both directions so the heatmap is symmetrical
pw_both <- bind_rows(pw, select(pw, pop1 = pop2, pop2 = pop1, mean_fst, weighted_fst)) %>%
  mutate(pop1 = factor(pop1, levels = pop_order),
         pop2 = factor(pop2, levels = pop_order))

ggplot(pw_both, aes(x = pop1, y = pop2, fill = weighted_fst)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(weighted_fst, 3)), size = 2.5) +
  scale_fill_gradient(low = "white", high = "red4", name = "Weighted\nFST") +
  labs(x = NULL, y = NULL,
       title = "Genome-wide pairwise FST (populations ordered south to north)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(OUTDIR,"/pairwise_fst_heatmap.png"), width = 7, height = 6)

##Does differentiation increase with the geographic distance between populations?

pw_dist <- pw %>%
  left_join(select(env, Population, lat1 = Latitude), by = c("pop1" = "Population")) %>%
  left_join(select(env, Population, lat2 = Latitude), by = c("pop2" = "Population")) %>%
  mutate(lat_diff = abs(lat1 - lat2))

ggplot(pw_dist, aes(x = lat_diff, y = weighted_fst)) +
  geom_point(size = 2) +
  geom_smooth(method = "lm", se = FALSE, colour = "grey40") +
  labs(x = "Difference in latitude (degrees)", y = "Weighted FST",
       title = "Isolation by distance?") +
  theme_classic()

ggsave(paste0(OUTDIR,"/fst_vs_distance.png"), width = 5, height = 4)

##

pop1 <- "ASU"   #the same two populations as in the bash code
pop2 <- "KRU"

fst_win <- read_tsv(paste0("fst/", pop1, "_vs_", pop2, "_FST_10kb.windowed.weir.fst"))

# Helper: turn (scaffold, position) into a position along the whole genome by laying the scaffolds end to end, in the order in which they appear in the data.
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

fst_win <- fst_win %>%
  mutate(POS = (BIN_START + BIN_END) / 2) %>% 
  add_genome_pos() ## this is our new function

# empirical outlier threshold: the top 1% of windows
fst_cut <- quantile(fst_win$WEIGHTED_FST, probs = 0.99)
fst_cut

fst_win <- mutate(fst_win, outlier = WEIGHTED_FST >= fst_cut)

ggplot(fst_win, aes(x = genome_pos / 1e6, y = WEIGHTED_FST)) +
  geom_point(aes(colour = shade), size = 0.8) +
  geom_point(data = filter(fst_win, outlier), colour = "red", size = 1.2) +
  geom_hline(yintercept = fst_cut, linetype = "dashed") +
  scale_colour_manual(values = c("grey35", "grey70"), guide = "none") +
  labs(x = "Position along the genome (Mb; scaffolds end to end)",
       y = "Weighted FST (10-kb windows)",
       title = paste(pop1, "vs", pop2)) +
  theme_classic()

ggsave(paste0(OUTDIR,"/fst_scan_", pop1, "_vs_", pop2, ".png"), width = 10, height = 4)