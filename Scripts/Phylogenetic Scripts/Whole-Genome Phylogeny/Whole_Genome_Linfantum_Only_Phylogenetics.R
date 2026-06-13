



# WHOLE-GENOME PHYLOGENY — L. infantum-only
# Reproducible script — only change BASE_DIR to run on another machine




# LOAD LIBRARIES ---------------------------------------------------------------

library(tidyverse)
library(ggtree)
library(treeio)
library(phytools)
library(ape)
library(ggtext)
library(ggnewscale)





# L. INFANTUM ONLY -------------------------------------------------------------

# Read tree and metadata files
tree <- read.newick("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Phylogenetic Analyses/Subset Files/Linfantum_samples_only/infantum_only_n463_Fast_bootstrap_tree.treefile")

# Read metadata and prepare it
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
    label_color = ifelse(highlight_zero, "tomato3", "dodgerblue3"),
    host_msl_label = case_when(
      !is.na(host_species) & !is.na(coverage_rounded) ~ paste0(host_species, " (", coverage_rounded, ")"),
      is.na(host_species) & !is.na(coverage_rounded) ~ paste0("NA (", coverage_rounded, ")"),
      !is.na(host_species) & is.na(coverage_rounded) ~ paste0(host_species, " (NA)"),
      TRUE ~ "NA (NA)"
    )
  )

# Convert to ape format
tree_ape <- as.phylo(tree)

# Midpoint root 
rooted_tree <- midpoint.root(tree_ape)




# 1. CIRCULAR TREES ----------------------------------------------------------

circular_tree <- ggtree(rooted_tree, layout = "circular", branch.length = "none") %<+% metadata

# Add bootstrap values
circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree$node.label



# Create circular tree coloured by country
linfantum_only_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) + 
  scale_color_discrete(name = "Country", na.translate = FALSE) +
  new_scale_color() +
  geom_tiplab(
    aes(label = host_msl_label,
        fontface = label_fontface,
        color = label_color),
    size = 2,
    show.legend = FALSE
  ) +
  scale_color_manual(values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"), guide = "none") +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 24),
    legend.key.height = unit(1.5, "cm"),
    legend.spacing.x = unit(1, "cm"), 
    legend.text = element_text(size = 20),
    plot.title = element_text(size = 26, face = "bold", hjust = 0.5)) +
  guides(
    color = guide_legend(
    ncol = 2,
    byrow = TRUE,
    keyheight = unit(1.2, "cm")
    )
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree (Whole-Genome) by Country\nWith Host Origin and CNV Coverage")



# Save the circular country-coloured tree as jpeg and pdf images 
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Country/Linfantum_only_circular_tree_by_country_with_host_and_cnv.jpeg",
       plot = linfantum_only_circular_tree_by_country_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Country/Linfantum_only_circular_tree_by_country_with_host_and_cnv.pdf",
       plot = linfantum_only_circular_tree_by_country_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)



# Create circular tree coloured by region
linfantum_only_circular_tree_by_region_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_geographic_region), size = 1) +
  scale_color_discrete(
    name = "Geographic Region",
    na.translate = FALSE,
    labels = c(
      "Central.East.Asia" = "Central East Asia",
      "North.Africa"      = "North Africa",
      "Middle.East"       = "Middle East",
      "Europe"            = "Europe",
      "Americas"          = "Americas"
    )
  ) +
  geom_point(data = . %>% filter(!isTip),
             aes(color = source_geographic_region),
             size = 0,
             show.legend = TRUE) +
  new_scale_color() +
  geom_tiplab(
    aes(label = host_msl_label, 
        fontface = label_fontface,
        color = label_color),
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide = "none"
  ) +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(
      size = 36,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 5)
    ),
    plot.margin = margin(t = -100, r = 20, b = 2, l = 2),
    panel.spacing = unit(0, "lines")  
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree (Whole-Genome) by Region\nWith Host Origin and CNV Coverage") +
  labs(color = "Geographic Region")


# Save the circular region-coloured tree as jpeg and pdf images 
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Region/Linfantum_only_circular_tree_by_region_with_host_and_cnv.jpeg",
       plot = linfantum_only_circular_tree_by_region_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Region/Linfantum_only_circular_tree_by_region_with_host_and_cnv.pdf",
       plot = linfantum_only_circular_tree_by_region_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)




# 2. STANDARD RECTANGULAR TREES ----------------------------------------------

# Add tip order from circular tree
tip_order <- circular_tree$data %>%
  filter(isTip) %>%
  pull(label)
  

# Base rectangular tree with metadata
rectangular_tree <- ggtree(rooted_tree, layout = "rectangular", branch.length = "none") %<+% metadata

# Add bootstrap values
rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree$node.label



# Add a column for tip label colours in metadata
metadata <- metadata %>%
  mutate(
    tip_label_color = ifelse(highlight_zero, "tomato3", "dodgerblue3")
  )

# Rebuild rectangular tree with the updated metadata
linfantum_only_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree + 
  geom_tree(aes(color = source_country), size = 1) +
  scale_color_discrete(
    name = "Country",
    na.translate = FALSE
  ) +
  new_scale_color() +
  geom_tiplab(
    aes(label = host_msl_label, fontface = label_fontface, color = I(ifelse(highlight_zero, "tomato3", "dodgerblue3"))),
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide = "none"
  ) +
  theme_tree2() +
  theme(
    legend.position = c(1, 0.5), 
    legend.justification = c(1, 0.5), 
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.height = unit(1.5, "cm"), 
    legend.spacing.y = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 10, r = 50, b = 10, l = 10)
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.3))) +  
  ggtitle("Rectangular L. infantum Phylogenetic Tree (Whole-Genome) by Country\nWith Host Origin and CNV Coverage")


# Save the rectangular country-coloured tree as jpeg and pdf images
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Country/Linfantum_only_rectangular_tree_by_country_with_host_and_cnv.jpeg",
       plot = linfantum_only_rectangular_tree_by_country_with_host_and_cnv,
       width = 35, 
       height = 80, 
       units = "in", 
       dpi = 300, 
       limitsize = FALSE)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Country/Linfantum_only_rectangular_tree_by_country_with_host_and_cnv.pdf",
       plot = linfantum_only_rectangular_tree_by_country_with_host_and_cnv,
       width = 35, 
       height = 80, 
       units = "in", 
       dpi = 300, 
       limitsize = FALSE)



# Create rectangular tree by region
linfantum_only_rectangular_tree_by_region_with_host_and_cnv <- rectangular_tree +
  geom_tree(aes(color = source_geographic_region), size = 1) +
  geom_tiplab(
    aes(label = host_msl_label, fontface = label_fontface, color = I(ifelse(highlight_zero, "tomato3", "dodgerblue3"))),
    size = 3,
    show.legend = FALSE
  ) +
  geom_point(
    data = . %>% filter(!isTip),
    aes(color = source_geographic_region),
    size = 0,
    show.legend = TRUE
  ) +
  theme_tree2() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.height = unit(1.5, "cm"),
    legend.spacing.y = unit(0.8, "cm"),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.x = element_blank(),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.margin = margin(t = 10, r = 10, b = 10, l = 10) 
  ) +
  scale_color_discrete(
    name = "Geographic Region",
    na.translate = FALSE,
    labels = c(
      "Central.East.Asia" = "Central East Asia",
      "North.Africa"      = "North Africa",
      "Middle.East"       = "Middle East",
      "Europe"            = "Europe",
      "Americas"          = "Americas"
    )
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree (Whole-Genome) by Region\nWith Host Origin and CNV Coverage")


# Save the rectangular country-coloured tree as jpeg and pdf images
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Region/Linfantum_only_rectangular_tree_by_region_with_host_and_cnv.jpeg",
       plot = linfantum_only_rectangular_tree_by_region_with_host_and_cnv,
       width = 30, 
       height = 80, 
       units = "in", 
       dpi = 300, 
       limitsize = FALSE)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/Whole-Genome/L.infantum-Only/By Region/Linfantum_only_rectangular_tree_by_region_with_host_and_cnv.pdf",
       plot = linfantum_only_rectangular_tree_by_region_with_host_and_cnv,
       width = 30, 
       height = 80, 
       units = "in", 
       dpi = 300, 
       limitsize = FALSE)



