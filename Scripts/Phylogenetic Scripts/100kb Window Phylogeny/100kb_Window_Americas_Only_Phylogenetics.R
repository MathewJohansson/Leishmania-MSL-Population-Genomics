



# 100kb WINDOW PHYLOGENY — Americas-only
# Reproducible script — only change BASE_DIR to run on another machine




# LOAD LIBRARIES ---------------------------------------------------------------

# Load required libraries
library(tidyverse)
library(ggtree)
library(treeio)
library(phytools)
library(ape)
library(ggtext)
library(ggnewscale)





# AMERICAS-ONLY ----------------------------------------------------------------

# Read tree and metadata files
tree <- read.newick("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/100kb_windows/diploid.LinJ.31.100kb.filt_vcf.dir/Phylogeny_analysis.dir/diploid.LinJ.31.100kb_Fast_bootstrap_tree.treefile")

# Read metadata and prepare it to match tree tip labels (ERR IDs)
metadata <- read_tsv("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv") %>%
  rename(label = sra_run_accession) %>%
  select(label, who_id, source_geographic_region, source_country, source_state_or_region, host_species, coverage_rounded) %>%
  mutate(
    host_species = case_when(
      host_species %in% c("Human", "Homo sapiens") ~ "Human",
      host_species %in% c("Dog", "Canis familiaris") ~ "Dog",
      host_species == "Nyctomys procyonoides" ~ "Rodent",
      TRUE ~ host_species
    ),
    highlight_zero = coverage_rounded == 0,
    label_fontface = ifelse(highlight_zero, "bold", "plain"),
    label_color = ifelse(highlight_zero, "red", "blue"),
    label_fontface = ifelse(coverage_rounded == 0, "bold", "plain"),
    label_color = ifelse(coverage_rounded == 0, "tomato3", "dodgerblue3"),
    host_msl_label = case_when(
      !is.na(host_species) & !is.na(coverage_rounded) ~ paste0(host_species, " (", coverage_rounded, ")"),
      is.na(host_species) & !is.na(coverage_rounded) ~ paste0("NA (", coverage_rounded, ")"),
      !is.na(host_species) & is.na(coverage_rounded) ~ paste0(host_species, " (NA)"),
      TRUE ~ "NA (NA)"
    )
  )



# Filter metadata to Americas samples only
metadata_americas <- metadata %>% 
  filter(source_geographic_region == "Americas")

# Convert to ape format
tree_ape <- as.phylo(tree)

# Prune tree to match Americas-only tips
rooted_tree_americas <- drop.tip(tree_ape, setdiff(tree_ape$tip.label, metadata_americas$label))

# Ensure tip_color and label_fontface exist in metadata
metadata_americas <- metadata_americas %>%
  mutate(
    tip_color = ifelse(highlight_zero, "tomato3", "dodgerblue3"),
    label_fontface = ifelse(highlight_zero, "bold", "plain")
  )




# 1. CIRCULAR TREE -------------------------------------------------------------

# Create the circular tree with metadata
circular_tree <- ggtree(rooted_tree_americas, layout = "circular", branch.length = "none") %<+% metadata_americas

# Add bootstrap values
circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree_americas$node.label

# Create the circular tree coloured by country with tip labels
americas_only_100kb_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) +
  geom_point(data = . %>% filter(!isTip), aes(color = source_country), size = 0, show.legend = TRUE) +
  scale_color_manual(
    values = c(
      "Brazil"   = "#F8766D",
      "Honduras" = "#A3A500",
      "Panama"   = "#00BF7D",
      "Paraguay" = "#00B0F6",
      "Uruguay"  = "#E76EF4"
    ),
    na.translate = FALSE,
    name = "Country"
  ) +
  new_scale_color() + 
  geom_tiplab(
    aes(label = host_msl_label, fontface = label_fontface, color = tip_color),
    size = 2,
    show.legend = FALSE
  ) +
  scale_color_identity() + 
  theme_tree() +
  theme(
    legend.position = c(1.05, 0.5),
    legend.justification = c(1, 0.5),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 2, r = 200, b = 2, l = 2)
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree by Country (100kb MSL Window)\nfor American Strains Only, Host Origin, and CNV Coverage") +
  labs(color = "Country")

# Display the plot
americas_only_100kb_circular_tree_by_country_with_host_and_cnv

# Save the circular country-coloured tree as jpeg and pdf images 
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/100kb Windows/Americas-Only/By Country/Americas_only_100kb_circular_tree_by_country_with_host_and_cnv.jpeg",
       plot = americas_only_100kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/100kb Windows/Americas-Only/By Country/Americas_only_100kb_circular_tree_by_country_with_host_and_cnv.pdf",
       plot = americas_only_100kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)




# 2. STANDARD RECTANGULAR TREE -------------------------------------------------

# Create a standard rectangular tree
rectangular_tree <- ggtree(rooted_tree_americas, layout = "rectangular", branch.length = "none") %<+% metadata_americas

# Add bootstrap values
rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree_americas$node.label



# Create the standard rectangular tree coloured by country
americas_only_100kb_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree + 
  geom_tree(aes(color = source_country), size = 1.2, na.rm = TRUE) +
  scale_color_discrete(name = "Country", na.translate = FALSE) +
  new_scale_color() + 
  geom_tiplab(
    aes(
      label = host_msl_label,
      fontface = ifelse(coverage_rounded == 0, "bold", "plain"),
      color = ifelse(coverage_rounded == 0, "tomato3", "dodgerblue3")  
    ),
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide = "none" 
  ) +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0))) +
  guides(
    color = guide_legend(
      ncol = 1,
      byrow = FALSE,
      title = "Country",
      override.aes = list(size = 4, shape = 15)
    )
  ) +
  theme(
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.size = unit(2, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 20, margin = margin(r = 10)),
    legend.title = element_text(size = 26, margin = margin(b = 10))
  ) +
  ggtitle("Rectangular L. infantum Phylogeny (100kb MSL Window)\nAmerican Strains Only, Host Origin, and CNV Coverage") 

# Display the plot
americas_only_100kb_rectangular_tree_by_country_with_host_and_cnv


# Save the rectangular country-coloured tree as jpeg and pdf images
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/100kb Windows/Americas-Only/By Country/Americas_only_100kb_rectangular_tree_by_country_with_host_and_cnv.jpeg",
       plot = americas_only_100kb_rectangular_tree_by_country_with_host_and_cnv, 
       device = "jpeg",
       width = 25, 
       height = 80, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/100kb Windows/Americas-Only/By Country/Americas_only_100kb_rectangular_tree_by_country_with_host_and_cnv.pdf",
       plot = americas_only_100kb_rectangular_tree_by_country_with_host_and_cnv, 
       width = 25, 
       height = 80, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)



