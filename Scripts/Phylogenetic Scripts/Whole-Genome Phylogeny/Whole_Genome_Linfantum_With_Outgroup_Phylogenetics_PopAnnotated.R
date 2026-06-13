



# WHOLE-GENOME PHYLOGENY — L. infantum + L. donovani outgroup
# PopAnnotated version: adds ADMIXTURE K=12 population layer to tip points
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
DATA_DIR     <- file.path(BASE_DIR, "Data/phylogenetics")
TREE_DIR     <- file.path(BASE_DIR, "Trees/Whole-Genome")
FIGURES_DIR <- file.path(BASE_DIR, "Figures/Phylogenetics Work/Whole-Genome/L.infantum Plus Outgroup")
OUTPUTS_DIR  <- file.path(BASE_DIR, "Phylogeny Outputs/Whole-Genome/Linfantum_Plus_Outgroup")

# Create output folders if missing
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUTS_DIR, "sessionInfo.txt")
)

# Paul Tol palette for K=12 populations (LO1-LO12)
paul_tol_palette_k12 <- c(
  LO1  = "#4477AA", LO2  = "#EE6677", LO3  = "#228833", LO4  = "#CCBB44",
  LO5  = "#66CCEE", LO6  = "#AA3377", LO7  = "#BBBBBB", LO8  = "#EE7733",
  LO9  = "#009988", LO10 = "#882255", LO11 = "#332288", LO12 = "#AA7744"
)




# 3. FILE INPUTS ---------------------------------------------------------------

TREE_FILE <- file.path(
  BASE_DIR,
  "Phylogenetic Analyses/Subset Files/Linfantum_samples_with_outgroup",
  "Linfantum_samples_with_outgroup_2025-03-12_Fast_bootstrap_tree.treefile"
)

METADATA_FILE <- file.path(
  METADATA_DIR,
  "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"
)

# -------------------------------------------------------------------------


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

# Load K=12 population assignments and join to metadata
# NOTE: must happen after metadata is defined
pop_membership <- read_tsv(
  file.path(BASE_DIR, "Data/population_structure/population_membership_table_allsamples_k12_LO_notation.tsv")
) %>%
  rename(label = sample) %>%
  select(label, population)

metadata <- metadata %>%
  select(-any_of("population")) %>%
  left_join(pop_membership, by = "label") %>%
  mutate(population = factor(population, levels = paste0("LO", 1:12)))




# 4. ROOT TREE -----------------------------------------------------------------

outgroup_srr_ids <- c("ERR3550127", "ERR205758")
tree_ape    <- as.phylo(tree)
rooted_tree <- root(tree_ape, outgroup = outgroup_srr_ids, resolve.root = TRUE)
mrca_node   <- getMRCA(rooted_tree, outgroup_srr_ids)




# 5. CIRCULAR TREES ------------------------------------------------------------

# Base circular tree with metadata joined
circular_tree <- ggtree(rooted_tree, layout = "circular", branch.length = "none") %<+% metadata

# Add bootstrap values
circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree$node.label

mrca_node <- getMRCA(rooted_tree, outgroup_srr_ids)


# Circular tree coloured by country, with population tip points
linfantum_outgroup_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_country),
             size = 0,
             show.legend = TRUE) +
  scale_color_discrete(name = "Country", na.translate = FALSE) +
  new_scale_color() +
  geom_tiplab(
    aes(
      label    = host_msl_label,
      fontface = ifelse(highlight_zero, "bold", "plain"),
      color    = ifelse(highlight_zero, "tomato3", "dodgerblue3")
    ),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide  = "none") +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k12,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 5)),
    na.translate = FALSE) +
  geom_hilight(node = mrca_node, fill = "yellow", alpha = 0.3, extend = 1) +
  geom_cladelabel(node = mrca_node,
                  label    = "L. donovani outgroup",
                  offset   = 3.2,
                  barsize  = 0.5,
                  angle    = 0,
                  fontsize = 4) +
  theme_tree() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 24),
    legend.text      = element_text(size = 20),
    legend.key.size  = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title       = element_text(size = 26, face = "bold", hjust = 0.5, margin = margin(b = 0)),
    plot.margin      = margin(t = 2, r = 20, b = 2, l = 2)
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree (Whole-Genome) by Country\nwith L. donovani Outgroup, Host Origin, CNV Coverage, and ADMIXTURE Population")

linfantum_outgroup_circular_tree_by_country_with_host_and_cnv

# Save the circular country-coloured tree as jpeg and pdf images
ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_outgroup_circular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_outgroup_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_outgroup_circular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_outgroup_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)




# Circular tree coloured by region, with population tip points
linfantum_outgroup_circular_tree_by_region_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_geographic_region), size = 1) +
  scale_color_discrete(
    name   = "Geographic Region",
    na.translate = FALSE,
    labels = c(
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
    size = 3,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide  = "none") +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k12,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 5)),
    na.translate = FALSE) +
  geom_hilight(node = mrca_node, fill = "yellow", alpha = 0.3, extend = 1) +
  geom_cladelabel(node = mrca_node,
                  label    = "L. donovani outgroup",
                  offset   = 3.2,
                  barsize  = 0.5,
                  angle    = 0,
                  fontsize = 4) +
  theme_tree() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 24),
    legend.text      = element_text(size = 20),
    legend.key.size  = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title       = element_text(size = 36, face = "bold", hjust = 0.5, margin = margin(b = 5)),
    plot.margin      = margin(t = -100, r = 20, b = 2, l = 200),
    panel.spacing    = unit(0, "lines")
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree (Whole-Genome) by Region\nwith L. donovani Outgroup, Host Origin, CNV Coverage, and ADMIXTURE Population")

linfantum_outgroup_circular_tree_by_region_with_host_and_cnv

# Save the circular region-coloured tree as jpeg and pdf images
ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_outgroup_circular_tree_by_region_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_outgroup_circular_tree_by_region_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)

ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_outgroup_circular_tree_by_region_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_outgroup_circular_tree_by_region_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)




# 6. RECTANGULAR TREES ---------------------------------------------------------

# Clean country names for rectangular plots
metadata_rect <- metadata %>%
  mutate(source_country = str_replace_all(source_country, "\\.", " "))

rectangular_tree <- ggtree(rooted_tree, layout = "rectangular", branch.length = "none") %<+% metadata_rect

rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree$node.label

mrca_node <- getMRCA(rooted_tree, outgroup_srr_ids)


# Rectangular tree coloured by country, with population tip points
linfantum_outgroup_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree +
  geom_tree(aes(color = source_country), size = 1.2, na.rm = TRUE) +
  scale_color_discrete(name = "Country", na.translate = FALSE, guide = guide_legend(ncol = 2)) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = label_fontface,
        color    = label_color),
    size = 3,
    show.legend = FALSE) +
  scale_color_identity() +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k12,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 4)),
    na.translate = FALSE) +
  geom_hilight(node = mrca_node, fill = "yellow", alpha = 0.3, extend = 1) +
  geom_cladelabel(node = mrca_node,
                  label    = "L. donovani outgroup",
                  offset   = 3.2,
                  barsize  = 0.5,
                  angle    = 0,
                  fontsize = 4) +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 6), ncol = 1)) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.line.x       = element_blank(),
    axis.title.x      = element_blank(),
    legend.position   = c(0.05, 0.5),
    legend.justification = c(0, 0.5),
    legend.title      = element_text(size = 20, face = "plain"),
    legend.text       = element_text(size = 18),
    legend.key.height = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.8, "cm"),
    plot.title        = element_text(size = 40, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree (Whole-Genome) by Country\nwith L. donovani Outgroup, Host Origin, CNV Coverage, and ADMIXTURE Population")

linfantum_outgroup_rectangular_tree_by_country_with_host_and_cnv

# Save the rectangular country-coloured tree as jpeg and pdf images
ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_outgroup_rectangular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_outgroup_rectangular_tree_by_country_with_host_and_cnv,
       width = 25, height = 80, units = "in", dpi = 300, limitsize = FALSE)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Linfantum_outgroup_rectangular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_outgroup_rectangular_tree_by_country_with_host_and_cnv,
       width = 25, height = 80, units = "in", dpi = 300, limitsize = FALSE)




# Clean region names and rebuild tree for region plot
metadata_rect_region <- metadata_rect %>%
  mutate(source_geographic_region = str_replace_all(source_geographic_region, "\\.", " "))

rectangular_tree_region <- ggtree(rooted_tree, layout = "rectangular", branch.length = "none") %<+% metadata_rect_region


# Rectangular tree coloured by region, with population tip points
linfantum_outgroup_rectangular_tree_by_region_with_host_and_cnv <- rectangular_tree_region +
  geom_tree(aes(color = source_geographic_region), size = 1.2, na.rm = TRUE) +
  scale_color_discrete(name = "Geographic Region", na.translate = FALSE, guide = guide_legend(ncol = 2)) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = label_fontface,
        color    = label_color),
    size = 3,
    show.legend = FALSE) +
  scale_color_identity() +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k12,
    name         = "Population",
    guide = guide_legend(ncol = 4, override.aes = list(size = 4)),
    na.translate = FALSE) +
  geom_hilight(node = mrca_node, fill = "yellow", alpha = 0.3, extend = 1) +
  geom_cladelabel(node = mrca_node,
                  label    = "L. donovani outgroup",
                  offset   = 3.2,
                  barsize  = 0.5,
                  angle    = 0,
                  fontsize = 4) +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 6), ncol = 1)) +
  theme(
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.line.x       = element_blank(),
    axis.title.x      = element_blank(),
    legend.position   = c(0.1, 0.5),
    legend.justification = c(0, 0.5),
    legend.title      = element_text(size = 20, face = "plain"),
    legend.text       = element_text(size = 18),
    legend.key.height = unit(1.2, "cm"),
    legend.spacing.y  = unit(0.8, "cm"),
    plot.title        = element_text(size = 40, face = "bold", hjust = 0.5)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree (Whole-Genome) by Region\nwith L. donovani Outgroup, Host Origin, CNV Coverage, and ADMIXTURE Population")

linfantum_outgroup_rectangular_tree_by_region_with_host_and_cnv

# Save the rectangular region-coloured tree as jpeg and pdf images
ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_outgroup_rectangular_tree_by_region_with_host_and_cnv_popannotated.jpeg"),
       plot = linfantum_outgroup_rectangular_tree_by_region_with_host_and_cnv,
       width = 25, height = 80, units = "in", dpi = 300, limitsize = FALSE)

ggsave(file.path(FIGURES_DIR, "By Region",
                 "Linfantum_outgroup_rectangular_tree_by_region_with_host_and_cnv_popannotated.pdf"),
       plot = linfantum_outgroup_rectangular_tree_by_region_with_host_and_cnv,
       width = 25, height = 80, units = "in", dpi = 300, limitsize = FALSE)



