







# 100kb WINDOW PHYLOGENY — L. infantum-only
# PopAnnotated version: adds ADMIXTURE K=13 population layer to tip points
# Reproducible script — only change BASE_DIR to run on another machine




# 1. LOAD LIBRARIES ------------------------------------------------------------

library(tidyverse)
library(ggtree)
library(treeio)
library(phytools)
library(ape)
library(ggtext)
library(ggnewscale)




# 2. PROJECT DIRECTORY CONFIG --------------------------------------------------

# Change ONLY this line when moving machines
BASE_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

# Define reusable directories
METADATA_DIR <- file.path(BASE_DIR, "Metadata")
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/Phylogenetics Work/100kb Windows/L.infantum-Only")
OUTPUTS_DIR  <- file.path(BASE_DIR, "Phylogeny Outputs/100kb Windows/Linfantum_Only")

# Create output folders if missing
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUTS_DIR, "sessionInfo.txt")
)

# Paul Tol palette for K=13 populations (LI1-LI13)
paul_tol_palette_k13 <- c(
  LI1  = "#4477AA", LI2  = "#EE6677", LI3  = "#228833", LI4  = "#CCBB44",
  LI5  = "#66CCEE", LI6  = "#AA3377", LI7  = "#BBBBBB", LI8  = "#EE7733",
  LI9  = "#009988", LI10 = "#882255", LI11 = "#332288", LI12 = "#AA7744",
  LI13 = "#44BB99"
)




# 3. FILE INPUTS ---------------------------------------------------------------

TREE_FILE <- file.path(
  BASE_DIR,
  "Data/Phylogenetic_Data/100kb_windows/diploid.LinJ.31.100kb.filt_vcf.dir/Phylogeny_analysis.dir",
  "diploid.LinJ.31.100kb_Fast_bootstrap_tree.treefile"
)

METADATA_FILE <- file.path(
  METADATA_DIR,
  "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"
)

stopifnot(file.exists(TREE_FILE))
stopifnot(file.exists(METADATA_FILE))

# Read tree
tree <- read.newick(TREE_FILE)

# Read metadata
metadata <- read_tsv(METADATA_FILE) %>%
  rename(label = sra_run_accession) %>%
  select(label, who_id, source_geographic_region, source_country,
         source_state_or_region, host_species, coverage_rounded) %>%
  mutate(
    host_species = case_when(
      host_species %in% c("Human", "Homo sapiens") ~ "Human",
      host_species %in% c("Dog", "Canis familiaris") ~ "Dog",
      host_species == "Nyctomys procyonoides" ~ "Rodent",
      TRUE ~ host_species
    ),
    highlight_zero  = coverage_rounded == 0,
    label_fontface  = ifelse(coverage_rounded == 0, "bold", "plain"),
    label_color     = ifelse(coverage_rounded == 0, "tomato3", "dodgerblue3"),
    host_msl_label  = case_when(
      !is.na(host_species) & !is.na(coverage_rounded) ~ paste0(host_species, " (", coverage_rounded, ")"),
      is.na(host_species)  & !is.na(coverage_rounded) ~ paste0("NA (", coverage_rounded, ")"),
      !is.na(host_species) & is.na(coverage_rounded)  ~ paste0(host_species, " (NA)"),
      TRUE ~ "NA (NA)"
    )
  )

# Load K=13 population assignments and join to metadata
pop_membership <- read_tsv(
  file.path(BASE_DIR, "Data/population_structure/population_membership_table_infantum_only_k13_LI_notation.tsv")
) %>%
  rename(label = sample) %>%
  select(label, population)

metadata <- metadata %>%
  select(-any_of("population")) %>%
  left_join(pop_membership, by = "label") %>%
  mutate(population = factor(population, levels = paste0("LI", 1:13)))




# 4. ROOT TREE -----------------------------------------------------------------

tree_ape    <- as.phylo(tree)
rooted_tree <- midpoint.root(tree_ape)




# 5. CIRCULAR TREES ------------------------------------------------------------

circular_tree <- ggtree(rooted_tree, layout = "circular", branch.length = "none") %<+% metadata

circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree$node.label


# Circular tree coloured by country, with population tip points
linfantum_100kb_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_country),
             size = 0,
             show.legend = TRUE) +
  scale_color_discrete(name = "Country", na.translate = FALSE, guide = guide_legend(ncol = 2)) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = ifelse(highlight_zero, "bold", "plain"),
        color    = ifelse(highlight_zero, "tomato3", "dodgerblue3")),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide  = "none") +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k13,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 5)),
    na.translate = FALSE) +
  theme_tree() +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 24),
    legend.text       = element_text(size = 20),
    legend.key.size   = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y  = unit(1, "cm"),
    plot.title        = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 0)),
    plot.margin       = margin(t = 2, r = 20, b = 2, l = 2)
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree by Country (100kb MSL Window)\nwith Host Origin, CNV Coverage, and ADMIXTURE Population")

# Display the plot
linfantum_100kb_circular_tree_by_country_with_host_and_cnv

# Save the circular country-coloured tree as jpeg and pdf images 
ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_Only_100kb_circular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_100kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_Only_100kb_circular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_100kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)



# Circular tree coloured by region, with population tip points
linfantum_100kb_circular_tree_by_region_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_geographic_region), size = 1) +
  scale_color_discrete(
    name         = "Geographic Region",
    na.translate = FALSE,
    guide        = guide_legend(ncol = 1),
    labels       = c(
      "Central.East.Asia" = "Central East Asia",
      "East.Africa"       = "East Africa",
      "Middle.East"       = "Middle East",
      "North.Africa"      = "North Africa"
    )
  ) +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_geographic_region),
             size = 0,
             show.legend = TRUE) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = label_fontface,
        color    = label_color),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide  = "none") +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k13,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 5)),
    na.translate = FALSE) +
  theme_tree() +
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 24),
    legend.text       = element_text(size = 20),
    legend.key.size   = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y  = unit(1, "cm"),
    plot.title        = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 0)),
    plot.margin       = margin(t = 2, r = 20, b = 2, l = 2)
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree by Region (100kb MSL Window)\nwith Host Origin, CNV Coverage, and ADMIXTURE Population")

# Display the plot
linfantum_100kb_circular_tree_by_region_with_host_and_cnv

# Save the circular region-coloured tree as jpeg and pdf images 
ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_Only_100kb_circular_tree_by_region_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_100kb_circular_tree_by_region_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)

ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_Only_100kb_circular_tree_by_region_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_100kb_circular_tree_by_region_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)




# 6. RECTANGULAR TREES ---------------------------------------------------------

rectangular_tree <- ggtree(rooted_tree, layout = "rectangular", branch.length = "none") %<+% metadata

rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree$node.label


# Rectangular tree coloured by country, with population tip points
linfantum_100kb_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree +
  geom_tree(aes(color = source_country), size = 1.2, na.rm = TRUE) +
  scale_color_discrete(name = "Country", na.translate = FALSE, guide = guide_legend(ncol = 2)) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = ifelse(coverage_rounded == 0, "bold", "plain"),
        color    = ifelse(coverage_rounded == 0, "tomato3", "dodgerblue3")),
    size = 3,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide  = "none") +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k13,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 4)),
    na.translate = FALSE) +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    legend.position      = "right",
    legend.direction     = "vertical",
    legend.box           = "vertical",
    legend.key.size      = unit(2, "cm"),
    legend.spacing.y     = unit(0.5, "cm"),
    legend.spacing.x     = unit(0.1, "cm"),
    legend.text          = element_text(size = 20, margin = margin(r = 10)),
    legend.title         = element_text(size = 26, margin = margin(b = 10)),
    plot.title           = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.margin          = margin(t = 10, r = 10, b = 10, l = 120)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree by Country (100kb MSL Window)\nwith Host Origin, CNV Coverage, and ADMIXTURE Population")

# Display the plot
linfantum_100kb_rectangular_tree_by_country_with_host_and_cnv

# Save the rectangular country-coloured tree as jpeg and pdf images
ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_Only_100kb_rectangular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_100kb_rectangular_tree_by_country_with_host_and_cnv,
       width = 30, height = 80, units = "in", dpi = 300, limitsize = FALSE)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_Only_100kb_rectangular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_100kb_rectangular_tree_by_country_with_host_and_cnv,
       width = 30, height = 80, units = "in", dpi = 300, limitsize = FALSE)



# Rectangular tree coloured by region, with population tip points
linfantum_100kb_rectangular_tree_by_region_with_host_and_cnv <- rectangular_tree +
  geom_tree(aes(color = source_geographic_region), size = 1.2, na.rm = TRUE) +
  scale_color_discrete(
    name         = "Geographic Region",
    na.translate = FALSE,
    guide        = guide_legend(ncol = 1),
    labels       = c(
      "Central.East.Asia" = "Central East Asia",
      "East.Africa"       = "East Africa",
      "Middle.East"       = "Middle East",
      "North.Africa"      = "North Africa"
    )
  ) +
  new_scale_color() +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_geographic_region),
             size = 0,
             show.legend = TRUE) +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = ifelse(highlight_zero, "bold", "plain"),
        color    = ifelse(highlight_zero, "tomato3", "dodgerblue3")),
    size = 3,
    show.legend = FALSE) +
  scale_color_identity() +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k13,
    name         = "Population",
    guide = guide_legend(ncol = 2, override.aes = list(size = 4)),
    na.translate = FALSE) +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +
  theme(
    axis.text.x          = element_blank(),
    axis.ticks.x         = element_blank(),
    axis.line.x          = element_blank(),
    axis.title.x         = element_blank(),
    legend.position      = c(0.08, 0.5),
    legend.justification = c(0, 0.5),
    legend.direction     = "vertical",
    legend.box           = "vertical",
    legend.key.size      = unit(2, "cm"),
    legend.key.width     = unit(3, "cm"),
    legend.key.height    = unit(2, "cm"),
    legend.spacing.y     = unit(0.5, "cm"),
    legend.spacing.x     = unit(0.1, "cm"),
    legend.text          = element_text(size = 20, margin = margin(r = 10)),
    legend.title         = element_text(size = 26, margin = margin(b = 10)),
    plot.title           = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.margin          = margin(t = 10, r = 10, b = 10, l = 50)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree by Region (100kb MSL Window)\nwith Host Origin, CNV Coverage, and ADMIXTURE Population")

# Display the plot
linfantum_100kb_rectangular_tree_by_region_with_host_and_cnv

# Save the rectangular region-coloured tree as jpeg and pdf images
ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_Only_100kb_rectangular_tree_by_region_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_100kb_rectangular_tree_by_region_with_host_and_cnv,
       width = 30, height = 80, units = "in", dpi = 300, limitsize = FALSE)

ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_Only_100kb_rectangular_tree_by_region_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_100kb_rectangular_tree_by_region_with_host_and_cnv,
       width = 30, height = 80, units = "in", dpi = 300, limitsize = FALSE)



