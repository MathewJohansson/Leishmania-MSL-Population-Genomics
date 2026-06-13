



# 20kb WINDOW PHYLOGENY — L. infantum-only
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





# L. INFANTUM-ONLY --------------------------------------------------

# Read tree and metadata files
tree <- read.newick("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/20kb_windows/diploid.LinJ.31.20kb.filt_vcf.dir/Phylogeny_analysis.dir/diploid.LinJ.31.20kb_Fast_bootstrap_tree.treefile")

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
    label_color = ifelse(highlight_zero, "red", "black"),
    label_fontface = ifelse(coverage_rounded == 0, "bold", "plain"),
    label_color = ifelse(coverage_rounded == 0, "tomato3", "dodgerblue3"),
    host_msl_label = case_when(
      !is.na(host_species) & !is.na(coverage_rounded) ~ paste0(host_species, " (", coverage_rounded, ")"),
      is.na(host_species) & !is.na(coverage_rounded) ~ paste0("NA (", coverage_rounded, ")"),
      !is.na(host_species) & is.na(coverage_rounded) ~ paste0(host_species, " (NA)"),
      TRUE ~ "NA (NA)"
    )
  )

# Convert tree to ape format if needed
tree_ape <- as.phylo(tree)

# Root the tree at midpoint (appropriate for intraspecific trees without outgroups)
rooted_tree <- midpoint.root(tree_ape)




# 1. CIRCULAR TREES ----------------------------------------------------------

# Create a circular tree
circular_tree <- ggtree(rooted_tree, layout = "circular", branch.length = "none") %<+% metadata

# Add bootstrap values
circular_tree$data$bootstrap <- NA
circular_tree$data$bootstrap[!circular_tree$data$isTip] <- rooted_tree$node.label




# Create the circular tree coloured by country
linfantum_20kb_circular_tree_by_country_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_country), size = 1) + 
  geom_point(data = . %>% filter(!isTip), 
             aes(color = source_country),
             size = 0,
             show.legend = TRUE) +
  scale_color_discrete(name = "Country", 
                       na.translate = FALSE) + 
  new_scale_color() +
  geom_tiplab(
    aes(label = host_msl_label, 
        fontface = label_fontface,
        color = label_color),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide = "none") +
  theme_tree() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(
      size = 26,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 0)
    ),
    plot.margin = margin(t = 2, r = 20, b = 2, l = 2)
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree by Country (20kb MSL Window)\nHost Origin and CNV Coverage") +
  labs(color = "Country")

# Display the plot
linfantum_20kb_circular_tree_by_country_with_host_and_cnv


# Save the circular country-coloured tree as jpeg and pdf images 
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Country/Linfantum_Only_20kb_circular_tree_by_country_with_host_and_cnv.jpeg",
       plot = linfantum_20kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Country/Linfantum_Only_20kb_circular_tree_by_country_with_host_and_cnv.pdf",
       plot = linfantum_20kb_circular_tree_by_country_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)




# Create the circular tree coloured by geographic region
linfantum_20kb_circular_tree_by_region_with_host_and_cnv <- circular_tree +
  geom_tree(aes(color = source_geographic_region), size = 1) +
  scale_color_discrete(
    name = "Geographic Region",
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
    aes(label = host_msl_label, 
        fontface = label_fontface,
        color = label_color),
    size = 2,
    show.legend = FALSE) +
  scale_color_manual(
    values = c("tomato3" = "tomato3", "dodgerblue3" = "dodgerblue3"),
    guide = "none") +
  theme_tree() +
  theme(
    legend.position = c(1.09, 0.5), 
    legend.justification = c(1, 0.5),  
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 20),
    legend.key.size = unit(2, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(1, "cm"),
    plot.title = element_text(
      size = 26,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 0)
    ),
    plot.margin = margin(t = 2, r = 200, b = 2, l = 2) 
  ) +
  ggtitle("Circular L. infantum Phylogenetic Tree by Region (20kb MSL Window)\nHost Origin and CNV Coverage") + 
  labs(color = "Geographic Region")


# Display the plot
linfantum_20kb_circular_tree_by_region_with_host_and_cnv


# Save the circular region-coloured tree as jpeg and pdf images
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Region/Linfantum_Only_20kb_circular_tree_by_region_with_host_and_cnv.jpeg", 
       plot = linfantum_20kb_circular_tree_by_region_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Region/Linfantum_Only_20kb_circular_tree_by_region_with_host_and_cnv.pdf",
       plot = linfantum_20kb_circular_tree_by_region_with_host_and_cnv,
       width = 25, 
       height = 25, 
       units = "in", 
       dpi = 600)




# 2. STANDARD RECTANGULAR TREES ----------------------------------------------

# Create a standard rectangular tree
rectangular_tree <- ggtree(rooted_tree, layout = "rectangular", branch.length = "none") %<+% metadata

# Add bootstrap values
rectangular_tree$data$bootstrap <- NA
rectangular_tree$data$bootstrap[!rectangular_tree$data$isTip] <- rooted_tree$node.label



# Create the standard rectangular tree coloured by country
linfantum_20kb_rectangular_tree_by_country_with_host_and_cnv <- rectangular_tree +
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
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(
    color = guide_legend(
      ncol = 1,
      byrow = FALSE,
      title = "Country",
      override.aes = list(size = 2, linetype = 1)
    )
  ) +
  theme(
    legend.position = "right", 
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.size = unit(2, "cm"),
    legend.key.width = unit(3, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.box.margin = margin(0, 0, 0, 0),
    legend.margin = margin(0, 0, 0, 0),
    legend.text = element_text(size = 20, margin = margin(r = 10)),
    legend.title = element_text(size = 26, margin = margin(b = 10)),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.margin = margin(10)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree by Country (20kb MSL Window)\nHost Origin and CNV Coverage")

# Display the plot
linfantum_20kb_rectangular_tree_by_country_with_host_and_cnv


# Save the rectangular country-coloured tree as jpeg and pdf images
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Country/Linfantum_Only_20kb_rectangular_tree_by_country_with_host_and_cnv.jpeg",
       plot = linfantum_20kb_rectangular_tree_by_country_with_host_and_cnv, 
       width = 35, 
       height = 80, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Country/Linfantum_Only_20kb_rectangular_tree_by_country_with_host_and_cnv.pdf",
       plot = linfantum_20kb_rectangular_tree_by_country_with_host_and_cnv, 
       width = 35, 
       height = 80, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)




# Create the standard rectangular tree coloured by geographic region
linfantum_20kb_rectangular_tree_by_region_with_host_and_cnv <- rectangular_tree + 
  geom_tree(aes(color = source_geographic_region), size = 1.2, na.rm = TRUE) +
  scale_color_discrete(
    name = "Geographic Region",
    na.translate = FALSE,
    labels = c(
      "Central.East.Asia" = "Central East Asia",
      "East.Africa"       = "East Africa",
      "Middle.East"       = "Middle East",
      "North.Africa"      = "North Africa"
    )
  ) +
  new_scale_color() +
  geom_point(data = . %>% filter(!isTip & !is.na(source_geographic_region)), 
             aes(color = source_geographic_region),
             size = 0,
             show.legend = TRUE) +
  geom_tiplab(
    aes(label = host_msl_label,
        fontface = ifelse(highlight_zero, "bold", "plain"),
        color = ifelse(highlight_zero, "tomato3", "dodgerblue3")),
    size = 3,
    show.legend = FALSE
  ) +
  scale_color_identity() +
  theme_tree2() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  guides(
    color = guide_legend(
      override.aes = list(
        size = 2,
        linetype = 1
      ),
      title = "Geographic Region",
      ncol = 1
    )
  ) +
  theme(
    legend.position = "right", 
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.key.size = unit(2, "cm"),
    legend.key.width = unit(3, "cm"),
    legend.key.height = unit(2, "cm"),
    legend.spacing.y = unit(0.5, "cm"),
    legend.spacing.x = unit(0.1, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 26),
    plot.title = element_text(size = 40, face = "bold", hjust = 0.5),
    plot.margin = margin(10)
  ) +
  ggtitle("Rectangular L. infantum Phylogenetic Tree by Region (20kb MSL Window)\nWith Host Origin and CNV Coverage")

# Display the plot
linfantum_20kb_rectangular_tree_by_region_with_host_and_cnv


# Save the rectangular region-coloured tree as jpeg and pdf images
ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Region/Linfantum_Only_20kb_rectangular_tree_by_region_with_host_and_cnv.jpeg", 
       plot = linfantum_20kb_rectangular_tree_by_region_with_host_and_cnv,
       width = 30, 
       height = 80, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)

ggsave("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/20kb Windows/L.infantum-Only/By Region/Linfantum_Only_20kb_rectangular_tree_by_region_with_host_and_cnv.pdf",
       plot = linfantum_20kb_rectangular_tree_by_region_with_host_and_cnv,
       width = 30, 
       height = 80, 
       units = "in", 
       dpi = 300,
       limitsize = FALSE)



