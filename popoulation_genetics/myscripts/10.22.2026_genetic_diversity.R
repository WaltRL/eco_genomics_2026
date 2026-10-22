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
pi_lat_plot <- ggplot(pi_summary, aes(Latitude, mean_pi)) +
  geom_point(size = 4, shape =21, fill = "cadetblue2", stroke =1.5) +
  labs(x = "Latitude (°N)", y = "Nucleotide diversity", colour = "Site") +
  theme_classic()+
  theme(legend.position = "none")
coord_flip()

ggsave(paste0(OUTDIR, "/pi/pi_vs_latitude.png"), plot = pi_lat_plot, width = 7, height = 7)