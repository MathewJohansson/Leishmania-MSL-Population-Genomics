



# WHOLE-GENOME PHYLOGENY — Americas-only
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
DATA_DIR     <- file.path(BASE_DIR, "Data/phylogenetics")
TREE_DIR     <- file.path(BASE_DIR, "Trees/Whole-Genome")
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/Phylogenetics Work/Whole-Genome/Americas-Only")
OUTPUTS_DIR  <- file.path(BASE_DIR, "Phylogeny Outputs/Whole-Genome/Americas_Only")

# Create output folders if missing
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)

# Save session info
writeLines(
  capture.output(sessionInfo()),
  file.path(OUTPUTS_DIR, "sessionInfo.txt")
)

# Paul Tol palette for K=13 populations (AM1-AM13)
paul_tol_palette_k13 <- c(
  AM1  = "#4477AA", AM2  = "#EE6677", AM3  = "#228833", AM4  = "#CCBB44",
  AM5  = "#66CCEE", AM6  = "#AA3377", AM7  = "#BBBBBB", AM8  = "#EE7733",
  AM9  = "#009988", AM10 = "#882255", AM11 = "#994F00", AM12 = "#F032E6",
  AM13 = "#44BB99"
)




# 3. FILE INPUTS ---------------------------------------------------------------

TREE_FILE <- file.path(
  BASE_DIR,
  "Phylogenetic Analyses/Subset Files/Americas_samples",
  "Americas_samples_n388_Fast_bootstrap_tree.treefile"
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
  file.path(BASE_DIR, "Data/population_structure/population_membership_table_Americas_only_k13_AM_notation.tsv")
) %>%
  rename(label = sample) %>%
  select(label, population)

metadata <- metadata %>%
  select(-any_of("population")) %>%
  left_join(pop_membership, by = "label") %>%
  mutate(population = factor(population, levels = paste0("AM", 1:13)))




# 4. ROOT TREE -----------------------------------------------------------------

tree_ape <- as.phylo(tree)

# Prune tree to Americas-only samples
americas_labels         <- metadata$label
rooted_tree_americas    <- drop.tip(tree_ape, setdiff(tree_ape$tip.label, americas_labels))

# Filter metadata to match pruned tree tips
metadata_americas <- metadata %>%
  filter(label %in% rooted_tree_americas$tip.label)




# 5. CIRCULAR TREE -------------------------------------------------------------

circular_tree <- ggtree(rooted_tree_americas, layout = "circular", branch.length = "none") %<+% metadata_americas

# Add bootstrap values
circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree_americas$node.label


# Circular tree coloured by country, with population tip points
americas_only_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_country),
             size = 0,
             show.legend = TRUE) +
  scale_color_discrete(name = "Country", na.translate = FALSE, guide = guide_legend(ncol = 1)) +
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
  ggtitle("Circular L. infantum Phylogenetic Tree (Whole-Genome) by Country\nfor American Strains Only, with Host Origin, CNV Coverage, and ADMIXTURE Population")

americas_only_circular_tree_by_country_with_host_and_cnv

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_circular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = americas_only_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_circular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = americas_only_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)




# 6. RECTANGULAR TREE ----------------------------------------------------------

rectangular_tree <- ggtree(rooted_tree_americas, layout = "rectangular", branch.length = "none") %<+% metadata_americas

# Add bootstrap values
rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree_americas$node.label


# Rectangular tree coloured by country, with population tip points
americas_only_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree +
  geom_tree(aes(color = source_country), size = 1) +
  scale_color_manual(
    values = c(
      "Brazil"   = "#F8766D",
      "Honduras" = "#A3A500",
      "Panama"   = "#00BF7D",
      "Paraguay" = "#00B0F6",
      "Uruguay"  = "#E76EF4"
    ),
    na.translate = FALSE,
    name         = "Country",
    guide        = guide_legend(ncol = 1)
  ) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = label_fontface,
        color    = label_color),
    size = 3,
    show.legend = FALSE,
    offset = 0) +
  scale_color_identity() +
  new_scale_color() +
  geom_tippoint(aes(subset = isTip, color = population), size = 2) +
  scale_color_manual(
    values       = paul_tol_palette_k13,
    name         = "Population",
    guide = guide_legend(ncol = 2, override.aes = list(size = 5)),
    na.translate = FALSE) +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +
  theme(
    legend.position      = c(0.85, 0.5),
    legend.justification = c(0, 0.5),
    legend.title         = element_text(size = 26),
    legend.text          = element_text(size = 24),
    legend.key.size      = unit(0.8, "cm"),
    legend.key.height    = unit(1.8, "cm"),
    legend.key.width     = unit(2.5, "cm"),
    legend.spacing.y     = unit(1.2, "cm"),
    legend.spacing.x     = unit(1.2, "cm"),
    plot.title           = element_text(size = 35, face = "bold", hjust = 0.5),
    plot.margin          = margin(t = 10, r = 200, b = 10, l = 200)
  ) +
  ggtitle("Standard Rectangular L. infantum Phylogenetic Tree (Whole-Genome) by Country\nfor American Strains Only, with Host Origin, CNV Coverage, and ADMIXTURE Population")

americas_only_rectangular_tree_by_country_with_host_and_cnv

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_rectangular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = americas_only_rectangular_tree_by_country_with_host_and_cnv,
       width = 30, height = 60, units = "in", dpi = 300, limitsize = FALSE)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_rectangular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = americas_only_rectangular_tree_by_country_with_host_and_cnv,
       width = 30, height = 120, units = "in", dpi = 300, limitsize = FALSE)


# Summary - the above plots are:
# 1. Circular tree coloured by country, with ADMIXTURE K=13 population tip points
# 2. Rectangular tree coloured by country, with ADMIXTURE K=13 population tip points


