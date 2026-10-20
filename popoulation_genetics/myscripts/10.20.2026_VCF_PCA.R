library(tidyverse)

setwd("/users/w/r/wreganlo/projects/eco_genomics_2026/popoulation_genetics/mydata")

#####Plotting the PCA in R-------------------------------------------------------------

###How Much Variation does each PC explain?

# the .eigenvec file has no column names, so col_names = FALSE (they become X1, X2, ...)
pca <- read_table("Lg_LD_0.2.eigenvec", col_names = FALSE)

# the .eigenval file is just a list of numbers, so scan() reads it into a vector
eigenval <- scan("Lg_LD_0.2.eigenval")

# columns 1 and 2 are the family and individual IDs, which are identical
# thanks to --double-id, so drop column 2 (remember: [rows, columns])
pca <- pca[, -2]

# rename the columns: the first is the sample ID, the rest are PC1, PC2, ...
# paste0() joins text and numbers, e.g. paste0("PC", 1:3) gives "PC1" "PC2" "PC3"
names(pca)[1] <- "ID"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca) - 1))

head(pca)

# convert the eigenvalues to a percentage variance explained (pve)
pve <- data.frame(PC = seq_along(eigenval),
                  pve = eigenval / sum(eigenval) * 100)

# make a folder for our figures

pve_plot <- ggplot(pve, aes(PC, pve)) +
  geom_col() +
  labs(x = "Principal component", y = "Percentage variance explained") +
  theme_bw()
pve_plot

OUTDIR="~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure"

ggsave(pve_plot, file=paste0(OUTDIR, "/variance_0.2.png"), width = 6, height = 4)

###Making the PCA---------------------------------------------------------------

meta <- read_csv("sample_metadata.csv")
pca <- left_join(pca, meta, by = "ID")

# check every individual matched a row in the metadata: this should print 0
sum(is.na(pca$Site))

# order the sites from north to south, so the legend follows the coast
pca <- pca %>% mutate(Site = fct_reorder(Site, -Latitude))

# axis labels with the % variance explained, e.g. "PC1 (4.2%)"
pc1_lab <- paste0("PC1 (", round(pve$pve[1], 1), "%)")
pc2_lab <- paste0("PC2 (", round(pve$pve[2], 1), "%)")

pca_plot <- ggplot(pca, aes(PC1, PC2, colour = Site)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = pc1_lab, y = pc2_lab, colour = "Site (north to south)") +
  theme_bw()

OUTDIR="~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure"

ggsave(pca_plot, file = paste0(OUTDIR, "/Lg_0.2_pca.png"), width = 8, height = 6)

###Is there a pattern with latitude?--------------------------------------------

lat_plot <- pca %>%
  select(ID, Site, Latitude, PC1, PC2) %>% ##Select columnns to use
  pivot_longer(c(PC1, PC2), names_to = "PC", values_to = "score") %>%
  ggplot(aes(Latitude, score, colour = Site)) +
  geom_point(size = 2, alpha = 0.6) +
  facet_wrap(~PC, ncol = 1, scales = "free_y") + ##plots PC1 and PC2 seperately
  scale_colour_viridis_d(end = 0.9) + 
  labs(x = "Latitude (°N)", y = "PC score", colour = "Site (north to south)") +
  theme_bw()+
  coord_flip()

OUTDIR="~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure"

ggsave(lat_plot, file = paste0(OUTDIR, "/Lg_0.2_pca_latitude.png"), width = 8, height = 6)

###BONUS: 
library(patchwork) ##combined plots

world <- map_data("world")
sites <- pca %>% distinct(Site, Latitude, Longitude)

map_plot <- ggplot() +
  geom_polygon(data = world, aes(long, lat, group = group),
               fill = "grey90", colour = "grey60") +
  geom_point(data = sites, aes(Longitude, Latitude, colour = Site), size = 3) +
  scale_colour_viridis_d(end = 0.9) +
  # zoom to our sites, with a 2-degree margin all round
  coord_quickmap(xlim = range(sites$Longitude) + c(-2, 2),
                 ylim = range(sites$Latitude) + c(-2, 2)) +
  labs(x = "Longitude", y = "Latitude", colour = "Site (north to south)") +
  theme_bw()

# with patchwork, + puts plots side by side
# we hide the map's legend, as the PCA's legend already shows the sites
map_pca <- map_plot + theme(legend.position = "none") + pca_plot
map_pca

OUTDIR="~/projects/eco_genomics_2026/popoulation_genetics/myresults/population_structure"

ggsave(map_pca, file = paste0(OUTDIR, "/Lg_0.2_pca_map.png"), width = 8, height = 6)
