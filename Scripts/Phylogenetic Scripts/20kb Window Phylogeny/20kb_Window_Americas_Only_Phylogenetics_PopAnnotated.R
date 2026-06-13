



# 20kb WINDOW PHYLOGENY — Americas-only
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
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/Phylogenetics Work/20kb Windows/Americas-Only")
OUTPUTS_DIR  <- file.path(BASE_DIR, "Phylogeny Outputs/20kb Windows/Americas_Only")

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
  AM9  = "#009988", AM10 = "#882255", AM11 = "#332288", AM12 = "#AA7744",
  AM13 = "#44BB99"
)

# Fixed country colours for Americas plots
americas_country_colours <- c(
  "Brazil"   = "#F8766D",
  "Honduras" = "#A3A500",
  "Panama"   = "#00BF7D",
  "Paraguay" = "#00B0F6",
  "Uruguay"  = "#E76EF4"
)




# 3. FILE INPUTS ---------------------------------------------------------------

TREE_FILE <- file.path(
  BASE_DIR,
  "Data/Phylogenetic_Data/20kb_windows/diploid.LinJ.31.20kb.filt_vcf.dir/Phylogeny_analysis.dir",
  "diploid.LinJ.31.20kb_Fast_bootstrap_tree.treefile"
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
         source_state_or_region, host_species, coverage_rounded)

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

# Clean tip labels and metadata labels
tree_ape$tip.label <- trimws(gsub('^\"|\"$', '', tree_ape$tip.label))
metadata$label     <- trimws(gsub('^\"|\"$', '', metadata$label))

# Normalise source_country
metadata <- metadata %>%
  mutate(source_country = ifelse(is.na(source_country) | source_country == "", NA_character_, trimws(source_country)))

# Define Americas countries and prune tree
americas_country_list <- c("Brazil", "Paraguay", "Panama", "Honduras", "Uruguay")
americas_labels       <- metadata %>% filter(source_country %in% americas_country_list) %>% pull(label)
rooted_tree_americas  <- drop.tip(tree_ape, setdiff(tree_ape$tip.label, americas_labels))

# Build tip-ordered metadata_americas
metadata_americas <- tibble(label = rooted_tree_americas$tip.label) %>%
  left_join(metadata, by = "label") %>%
  mutate(
    host_species = case_when(
      host_species %in% c("Human", "Homo sapiens") ~ "Human",
      host_species %in% c("Dog", "Canis familiaris") ~ "Dog",
      host_species == "Nyctomys procyonoides" ~ "Rodent",
      TRUE ~ host_species
    ),
    highlight_zero  = coverage_rounded == 0,
    label_fontface  = ifelse(highlight_zero, "bold", "plain"),
    tip_color       = ifelse(highlight_zero, "tomato3", "dodgerblue3"),
    host_msl_label  = case_when(
      !is.na(host_species) & !is.na(coverage_rounded) ~ paste0(host_species, " (", coverage_rounded, ")"),
      is.na(host_species)  & !is.na(coverage_rounded) ~ paste0("NA (", coverage_rounded, ")"),
      !is.na(host_species) & is.na(coverage_rounded)  ~ paste0(host_species, " (NA)"),
      TRUE ~ "NA (NA)"
    ),
    source_country = ifelse(is.na(source_country) | source_country == "", "Unknown", source_country)
  )




# 5. CIRCULAR TREE -------------------------------------------------------------

circular_tree <- ggtree(rooted_tree_americas, layout = "circular", branch.length = "none") %<+% metadata_americas

circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree_americas$node.label


# Circular tree coloured by country, with population tip points
americas_only_20kb_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_country),
             size = 0,
             show.legend = TRUE) +
  scale_color_manual(
    values       = americas_country_colours,
    na.translate = FALSE,
    name         = "Country",
    guide        = guide_legend(ncol = 1)
  ) +
  new_scale_color() +
  geom_tiplab(
    aes(label    = host_msl_label,
        fontface = label_fontface,
        color    = tip_color),
    size = 2,
    show.legend = FALSE) +
  scale_color_identity() +
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
  ggtitle("Circular L. infantum Phylogenetic Tree by Country (20kb MSL Window)\nfor American Strains Only, with Host Origin, CNV Coverage, and ADMIXTURE Population")

# Display the plot
americas_only_20kb_circular_tree_by_country_with_host_and_cnv

# Save the circular country-coloured tree as jpeg and pdf images 
ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_20kb_circular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = americas_only_20kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 600)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_20kb_circular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = americas_only_20kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, height = 25, units = "in", dpi = 300, limitsize = FALSE)




# 6. RECTANGULAR TREE ----------------------------------------------------------

rectangular_tree <- ggtree(rooted_tree_americas, layout = "rectangular", branch.length = "none") %<+% metadata_americas

rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree_americas$node.label


# Rectangular tree coloured by country, with population tip points
americas_only_20kb_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree +
  geom_tree(aes(color = source_country), size = 1.2, na.rm = TRUE) +
  scale_color_manual(
    values       = americas_country_colours,
    na.translate = FALSE,
    name         = "Country",
    guide        = guide_legend(ncol = 1)
  ) +
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
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.02))) +
  coord_cartesian(clip = "off") +
  theme(
    legend.position      = "right",
    legend.direction     = "vertical",
    legend.box           = "vertical",
    legend.key.size      = unit(2, "cm"),
    legend.key.width     = unit(3, "cm"),
    legend.spacing.y     = unit(0.5, "cm"),
    legend.spacing.x     = unit(0.1, "cm"),
    legend.box.margin    = margin(0, 0, 0, 0),
    legend.margin        = margin(0, 0, 0, 0),
    legend.text          = element_text(size = 20, margin = margin(r = 10)),
    legend.title         = element_text(size = 26, margin = margin(b = 10)),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.title.position = "plot",
    plot.margin          = margin(20, 10, 20, 10)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree by Country (20kb MSL Window)\nfor American Strains Only, with Host Origin, CNV Coverage, and ADMIXTURE Population")

# Display the plot
americas_only_20kb_rectangular_tree_by_country_with_host_and_cnv

# Save the rectangular country-coloured tree as jpeg and pdf images 
ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_20kb_rectangular_tree_by_country_with_host_and_cnv_popannotated.jpeg"),
       plot = americas_only_20kb_rectangular_tree_by_country_with_host_and_cnv,
       width = 30, height = 60, dpi = 150, limitsize = FALSE)

ggsave(file.path(FIGURES_DIR, "By Country",
                 "Americas_only_20kb_rectangular_tree_by_country_with_host_and_cnv_popannotated.pdf"),
       plot = americas_only_20kb_rectangular_tree_by_country_with_host_and_cnv,
       width = 33, height = 120, units = "in", dpi = 300, limitsize = FALSE)



