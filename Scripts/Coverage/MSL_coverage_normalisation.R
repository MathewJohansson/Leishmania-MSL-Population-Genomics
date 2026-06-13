# MSL COVERAGE NORMALISATION
# Task: Normalise read depth relative to other chromosomes.

# Set median read depth for all non-chr31 genes to 2.
# Plot 1 dot per gene: 4 genes before MSL + 4 MSL genes + 4 genes after MSL.

# Data: diploid_coverage_genes.2025-03-12.RData

# MSL gene positions on chr31:
#   LINF_310031200: 1,181,282 – 1,182,328
#   LINF_310031300: 1,184,205 – 1,185,341
#   LINF_310031400: 1,185,827 – 1,188,553
#   LINF_310031500: 1,191,357 – 1,192,406




# 1. SETUP ---------------------------------------------------------------------

library(tidyverse)

BASE_DIR <- "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
DATA_DIR <- file.path(BASE_DIR, "Data")
FIGURES_DIR <- file.path(BASE_DIR, "Figures/MSL Coverage Normalisation")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

load(file.path(DATA_DIR, "diploid_coverage_genes.2025-03-12.RData"))




# 2. IDENTIFY FLANKING GENES ---------------------------------------------------

msl_genes <- c("LINF_310031200", "LINF_310031300", "LINF_310031400", "LINF_310031500")

# Get ordered chr31 gene list and find flanking genes
chr31_genes <- diploid_coverage_genes %>%
  filter(chrom == "LinJ.31") %>%
  arrange(start)

msl_indices  <- which(chr31_genes$gene %in% msl_genes)
first_msl    <- min(msl_indices)
last_msl     <- max(msl_indices)

flanking_indices <- (first_msl - 4):(last_msl + 4)
genes_of_interest <- chr31_genes$gene[flanking_indices]

cat("Genes of interest (4 before + 4 MSL + 4 after):\n")
print(chr31_genes[flanking_indices, c("gene", "start", "end")])




# 3. NORMALISE COVERAGE -------------------------------------------------------

# Sample columns start at column 7
sample_cols <- colnames(diploid_coverage_genes)[7:ncol(diploid_coverage_genes)]

# For each sample: calculate median depth across all non-chr31 genes,
# then normalise so that median = 2
non_chr31 <- diploid_coverage_genes %>%
  filter(chrom != "LinJ.31") %>%
  select(all_of(sample_cols))

sample_medians <- apply(non_chr31, 2, median, na.rm = TRUE)

# Normalise entire dataset
normalised <- diploid_coverage_genes %>%
  mutate(across(all_of(sample_cols),
                ~ . / sample_medians[cur_column()] * 2))




# 4. EXTRACT GENES OF INTEREST ------------------------------------------------

plot_data <- normalised %>%
  filter(gene %in% genes_of_interest) %>%
  select(gene, start, all_of(sample_cols)) %>%
  pivot_longer(cols = all_of(sample_cols),
               names_to  = "sample",
               values_to = "normalised_depth") %>%
  mutate(
    gene_type = case_when(
      gene %in% msl_genes ~ "MSL gene",
      TRUE                ~ "Flanking gene"
    ),
    gene = factor(gene, levels = genes_of_interest)
  )




# 5. PLOT — ALL SAMPLES -------------------------------------------------------

reference_lines <- list(
  geom_hline(yintercept = 0, linetype = "solid",
             colour = "grey60", linewidth = 0.5),
  geom_hline(yintercept = 2, linetype = "dashed",
             colour = "grey50", linewidth = 0.5),
  geom_hline(yintercept = 4, linetype = "solid",
             colour = "grey60", linewidth = 0.5)
)

msl_coverage_normalisation_plot <- ggplot(plot_data,
                                          aes(x = gene,
                                              y = normalised_depth,
                                              colour = gene_type)) +
  reference_lines +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  stat_summary(fun = median, geom = "point",
               shape = 18, size = 4, colour = "black") +
  scale_colour_manual(values = c("MSL gene"      = "#BB5566",
                                 "Flanking gene"  = "#4477AA"),
                      name = "Gene type") +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 8),
    plot.title       = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.key.size  = unit(0.8, "cm")
  ) +
  labs(
    title = "Normalised MSL Locus Coverage — All Samples",
    x     = "Gene",
    y     = "Normalised read depth (median non-chr31 = 2)"
  )

msl_coverage_normalisation_plot




# 6. SAVE — ALL SAMPLES -------------------------------------------------------

ggsave(file.path(FIGURES_DIR, "MSL_coverage_normalisation_all_samples.jpeg"),
       plot  = msl_coverage_normalisation_plot,
       width = 10, height = 6, dpi = 300)

ggsave(file.path(FIGURES_DIR, "MSL_coverage_normalisation_all_samples.pdf"),
       plot  = msl_coverage_normalisation_plot,
       width = 10, height = 6)




# 7. STRATIFY BY MSL STATUS ---------------------------------------------------

# Load metadata to get coverage_rounded (MSL status) per sample
metadata <- read_tsv(
  file.path("~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Metadata",
            "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv")
)

msl_status <- metadata %>%
  select(sra_run_accession, coverage_rounded) %>%
  filter(coverage_rounded %in% c(0, 4)) %>%
  rename(sample = sra_run_accession,
         msl_status = coverage_rounded) %>%
  mutate(msl_label = ifelse(msl_status == 0, "MSL0 (deletion)", "MSL4 (wild-type)"))

# Join MSL status to plot data — keeps only MSL0 and MSL4 samples
plot_data_stratified <- plot_data %>%
  inner_join(msl_status, by = "sample")




# 8. PLOT — STRATIFIED BY MSL STATUS ------------------------------------------

msl_coverage_stratified_plot <- ggplot(plot_data_stratified,
                                       aes(x = gene,
                                           y = normalised_depth,
                                           colour = gene_type)) +
  reference_lines +
  geom_jitter(width = 0.15, alpha = 0.3, size = 0.8) +
  stat_summary(fun = median, geom = "point",
               shape = 18, size = 4, colour = "black") +
  scale_colour_manual(values = c("MSL gene"      = "#BB5566",
                                 "Flanking gene"  = "#4477AA"),
                      name = "Gene type") +
  facet_wrap(~ msl_label, ncol = 1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 10, by = 2)) +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1, size = 8),
    plot.title         = element_text(hjust = 0.5),
    panel.grid.major.x = element_blank(),
    legend.key.size    = unit(0.8, "cm"),
    strip.text         = element_text(size = 11, face = "bold"),
    panel.spacing      = unit(1.5, "cm")
  ) +
  labs(
    title = "Normalised MSL Locus Coverage — MSL0 vs MSL4",
    x     = "Gene",
    y     = "Normalised read depth (median non-chr31 = 2)"
  )

msl_coverage_stratified_plot




# 9. SAVE — STRATIFIED --------------------------------------------------------

ggsave(file.path(FIGURES_DIR, "MSL_coverage_normalisation_stratified.jpeg"),
       plot  = msl_coverage_stratified_plot,
       width = 10, height = 10, dpi = 300)

ggsave(file.path(FIGURES_DIR, "MSL_coverage_normalisation_stratified.pdf"),
       plot  = msl_coverage_stratified_plot,
       width = 10, height = 10)