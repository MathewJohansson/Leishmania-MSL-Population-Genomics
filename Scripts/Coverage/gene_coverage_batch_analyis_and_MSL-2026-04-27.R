# MSL LOCUS COVERAGE PLOT — L. INFANTUM AND OUTGROUPS
# Adapted from: gene_coverage_batch_analyis_and_MSL-2025-03-22.R (Dan Jeffares)
# Loads from step5 RData checkpoint and generates final coverage barplot + boxplot.




# 1. SETUP ---------------------------------------------------------------------

rm(list = ls())
library(tidyverse)
library(ggpubr)

base_dir <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
setwd(base_dir)

# Create output directory if needed
# Output directory already exists: Figures/Read Coverage




# 2. LOAD DATA -----------------------------------------------------------------

load("Scripts/Coverage/gene_coverage_batch_analysis.step5.RData")




# 3. FILTER FOR MSL GENES ------------------------------------------------------

msl_genes <- c("LINF_310031200", "LINF_310031300", "LINF_310031400", "LINF_310031500")

diploid_coverage_genes_msl <- diploid_coverage_genes |>
  filter(gene %in% msl_genes)

dim(diploid_coverage_genes_msl)




# 4. CALCULATE MEDIAN COVERAGE PER STRAIN --------------------------------------

msl_median_coverage <- diploid_coverage_genes_msl |>
  summarise(across(7:ncol(diploid_coverage_genes_msl), median))

msl_median_coverage_df <- data.frame(
  strain   = names(msl_median_coverage),
  coverage = as.vector(t(msl_median_coverage))
)

# Round to nearest integer, cap at 4
msl_median_coverage_df <- msl_median_coverage_df |>
  mutate(coverage_rounded = as.factor(pmin(round(coverage), 4)))

# Sort strains by coverage_rounded for clean barplot ordering
msl_median_coverage_df <- msl_median_coverage_df |>
  arrange(coverage_rounded) |>
  mutate(strain = factor(strain, levels = strain))




# 5. DEFINE PAUL TOL PALETTE ---------------------------------------------------
# MSL0 (deletion, 0 copies) = red; MSL4 (wild-type, 4 copies) = blue
# MSL1/2/3 (artefactual partial deletions) = neutral greys

msl_colours <- c(
  "4" = "#364B9A",   # deep blue — wild-type
  "3" = "#98CAE1",   # light blue
  "2" = "#EAECCC",   # pale neutral
  "1" = "#F99858",   # orange
  "0" = "#A50026"    # deep red — deletion
)



# 6. BARPLOT — MEDIAN MSL COVERAGE PER STRAIN ----------------------------------

barplot_coverage <- ggplot(msl_median_coverage_df,
                           aes(x = strain, y = coverage, fill = coverage_rounded)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = msl_colours, name = "MSL copies",
                    labels = c("0 (deletion)", "1", "2", "3", "4 (wild-type)")) +
  theme_minimal() +
  labs(title = "MSL gene coverage (L. infantum and outgroups)",
       x = "Strain",
       y = "Median Coverage") +
  geom_hline(yintercept = c(1:4), linetype = "dashed", colour = "grey40") +
  theme(axis.text.x  = element_blank(),
        axis.ticks.x = element_blank())




# 7. BOXPLOT — COVERAGE BY MSL STATUS ------------------------------------------

boxplot_coverage <- ggplot(msl_median_coverage_df,
                           aes(x = coverage_rounded, y = coverage, fill = coverage_rounded)) +
  geom_boxplot() +
  scale_fill_manual(values = msl_colours, name = "MSL copies",
                    labels = c("0 (deletion)", "1", "2", "3", "4 (wild-type)")) +
  theme_minimal() +
  labs(title = "",
       x = "MSL copy number estimate",
       y = "MSL read coverage") +
  geom_hline(yintercept = c(1:4), linetype = "dashed", colour = "grey40")




# 8. COMBINE AND SAVE ----------------------------------------------------------

boxplot_barplot_coverage <- ggarrange(
  barplot_coverage, boxplot_coverage,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  common.legend = TRUE, legend = "right"
)

ggsave("Figures/Read Coverage/Linfantum_and_outgroups_msl_median_coverage_plot.jpeg",
       boxplot_barplot_coverage, width = 16, height = 12, dpi = 300)

ggsave("Figures/Read Coverage/Linfantum_and_outgroups_msl_median_coverage_plot.pdf",
       boxplot_barplot_coverage, width = 16, height = 12, dpi = 300)

print("MSL coverage plot saved to Figures/Read Coverage/")



