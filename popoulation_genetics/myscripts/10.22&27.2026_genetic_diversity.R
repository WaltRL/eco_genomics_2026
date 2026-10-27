library(tidyverse)

setwd("/gpfs1/home/w/r/wreganlo/projects/eco_genomics_2026/popoulation_genetics/mydata")

##we want to output results to myresults
OUTDIR="/gpfs1/home/w/r/wreganlo/projects/eco_genomics_2026/popoulation_genetics/myresults"

meta <- read_csv("sample_metadata.csv")

# one latitude per site, so we can order the sites north to south
sites <- meta %>%
  group_by(Population) %>%
  summarise(Latitude = mean(Latitude))

# the files to read: pi/BOD.windowed.pi, pi/PGR.windowed.pi, ...
pi_files <- list.files("pi", pattern = "windowed.pi$", full.names = TRUE)

pi <- pi_files %>%
  # name each file by its site, e.g. "pi/BOD.windowed.pi" -> "BOD"
  set_names(str_remove(basename(pi_files), fixed(".windowed.pi"))) %>%
  # read every file and stack them, putting the names in a new column, Site
  map_dfr(read_tsv, .id = "Population", show_col_types = FALSE) %>%
  left_join(sites, by = "Population") %>%
  mutate(Site = fct_reorder(Population, Latitude))

# summarise pi for each site
pi_summary <- pi %>%
  group_by(Population, Latitude) %>%
  summarise(n_windows = n(),
            mean_pi   = mean(PI),
            median_pi = median(PI),
            .groups = "drop")
pi_summary

dir.create(paste0(OUTDIR,"/pi"))

pi_plot <- ggplot(pi, aes(Site, PI, fill = Site)) +
  geom_violin(scale = "width", linewidth = 0.5, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3) +
  scale_fill_viridis_d(end = 0.9) +
  scale_y_log10() +
  labs(x = "Population", y = "log10(Nucleotide diversity)") +
  theme_classic() +
  theme(legend.position = "none") + 
  coord_flip()

ggsave(paste0(OUTDIR, "/pi/pi_by_site.png"), plot = pi_plot, width = 6, height = 6)

# does diversity change with latitude?
pi_lat_plot <- ggplot(pi_summary, aes(Latitude, mean_pi, fill = Latitude)) +
  geom_point(size = 4, shape =21, stroke =1.5) +
  labs(x = "Latitude (°N)", y = "Nucleotide diversity") +
  theme_classic()+
  theme(legend.position = "none")+
  scale_fill_viridis_c()+
  coord_flip()

ggsave(paste0(OUTDIR, "/pi/pi_vs_latitude.png"), plot = pi_lat_plot, width = 7, height = 7)

###skipping over Tajima's D-----------------------------------------------------

#####Diversity along the genome-------------------------------------------------
top_scaffold <- pi %>%
  count(CHROM, sort = TRUE) %>%
  slice(1) %>%
  pull(CHROM)

scan_plot <- pi %>%
  filter(CHROM == top_scaffold) %>%
  ggplot(aes(BIN_START / 1e6, PI, colour = Site)) +
  geom_line(alpha=0.3) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = paste(top_scaffold, "position (Mb)"), y = "Nucleotide diversity",
       colour = "Site") +
  theme_bw() +
  facet_wrap(~Site)

ggsave(paste0(OUTDIR, "/pi/pi_along_scaffold_facet.png"), plot = scan_plot, width = 12, height = 6)

###Tajima's D-------------------------------------------------------------------

taj_files <- list.files("TajD", pattern = "Tajima.D$", full.names = TRUE)

tajd <- taj_files %>%
  set_names(str_remove(basename(taj_files), fixed(".Tajima.D"))) %>%
  # vcftools writes "nan" when D can't be calculated, so treat that as missing (NA)
  map_dfr(read_tsv, .id = "Population", na = c("", "NA", "nan", "-nan"),
          show_col_types = FALSE) %>%
  filter(N_SNPS >= 5) %>%   # drop windows with too few SNPs to be reliable
  left_join(sites, by = "Population") %>%
  mutate(Site = fct_reorder(Population, Latitude))

tajd %>%
  group_by(Site) %>%
  summarise(n_windows = n(),
            mean_D    = mean(TajimaD, na.rm = TRUE),
            median_D  = median(TajimaD, na.rm = TRUE))

tajd_plot <- ggplot(tajd, aes(Site, TajimaD, fill = Site)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_violin() +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  scale_fill_viridis_d(end = 0.9) +
  labs(x = "Site", y = "Tajima's D") +
  theme_classic() +
  theme(legend.position = "none")+
  coord_flip()

ggsave(paste0(OUTDIR, "/TajD/tajimasD_by_site.png"), plot = tajd_plot, width = 4, height = 7)

###Observed Heterozygosity

#dir.create(paste0(OUTDIR, "/het"))

het_files <- list.files("het", pattern = "het$", full.names = TRUE)

het_all <- read_tsv(het_files, id = "file", show_col_types = FALSE) %>%
  rename(ID = INDV, O_HOM = `O(HOM)`, E_HOM = `E(HOM)`) %>%
  mutate(
    Pop = tools::file_path_sans_ext(basename(file)),  # e.g. het/SiteA.het -> SiteA
    Ho  = (N_SITES - O_HOM) / N_SITES
  ) %>%
  select(-file) %>%
  left_join(meta, by = "ID") %>%
  mutate(Site = fct_reorder(Site, -Latitude))


##lets make one of my favuorite types of plots (raincloud)
library(ggdist)
#install.packages("ggdist")


ho_plot <- ggplot(het_all, aes(Site, Ho, fill = Site)) +
  geom_jitter(width = 0.05, alpha = 0.6, shape = 21, stroke =1.5) +
  geom_violin(alpha = 0.6) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, alpha =0.5) +
  scale_fill_viridis_d(end = 0.9) +
  labs(x = "", y = "Observed heterozygosity") +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust =1))+
  theme(legend.position = "none")


ggsave(paste0(OUTDIR, "/het/Ho_by_site.png"), plot = ho_plot, width = 4, height = 12)


ho_plot <- ggplot(het_all, aes(x = Site, y = Ho, fill = Site)) +
  geom_jitter(width = 0.05, height = 0, shape = 21, stroke = 0.8, alpha = 0.6, size = 1.4) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.5,
               position = position_nudge(x = 0.25)) +
  ggdist::stat_halfeye(
    adjust = 0.5, width = 0.6, justification = -0.5,
    .width = 0, point_colour = NA, alpha = 0.7
  ) +
  scale_fill_viridis_d(end = 0.9) +
  coord_flip() +
  labs(x = "", y = "Observed heterozygosity") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(paste0(OUTDIR, "/het/Ho_by_site.png"), plot = ho_plot, width = 5, height = 8)


F_plot <- ggplot(het_all, aes(x = Site, y = F, fill = Site)) +
  geom_jitter(width = 0.05, height = 0, shape = 21, stroke = 0.8, alpha = 0.6, size = 1.4) +
  geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, alpha = 0.5,
               position = position_nudge(x = 0.25)) +
  ggdist::stat_halfeye(
    adjust = 0.5, width = 0.6, justification = -0.5,
    .width = 0, point_colour = NA, alpha = 0.7
  ) +
  scale_fill_viridis_d(end = 0.9) +
  coord_flip() +
  labs(x = "", y = "F Coefficient") +
  theme_classic() +
  theme(legend.position = "none")

ggsave(paste0(OUTDIR, "/het/F_by_site.png"), plot = F_plot, width = 5, height = 8)