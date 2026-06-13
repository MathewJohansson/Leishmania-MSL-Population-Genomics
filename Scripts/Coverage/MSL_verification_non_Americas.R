# MSL COVERAGE VERIFICATION
# Task: Verify all non-American strains with MSL copy number 2 or 3.




# 1. SETUP ---------------------------------------------------------------------

library(tidyverse)

BASE_DIR <- "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
DATA_DIR <- file.path(BASE_DIR, "Data")

load(file.path(DATA_DIR, "diploid_coverage_genes.2025-03-12.RData"))




# 2. IDENTIFY CHR 31 GENE WINDOW --------------------------------------------------

chr31_genes <- diploid_coverage_genes %>%
  filter(chrom == "LinJ.31") %>%
  arrange(start)

msl_genes   <- c("LINF_310031200", "LINF_310031300", "LINF_310031400", "LINF_310031500")
msl_indices <- which(chr31_genes$gene %in% msl_genes)
first_msl   <- min(msl_indices)
last_msl    <- max(msl_indices)

window_indices <- (first_msl - 10):(last_msl + 10)
cat("Window genes:\n")
print(chr31_genes[window_indices, c("gene", "start", "end")])




# 3. LOAD SAMPLE LIST ---------------------------------------------------------

non_americas <- read_tsv(file.path(DATA_DIR,
                                   "samples_outside_Americas_with_less_than_4_MSL_copies.tsv"))

genes_of_interest <- chr31_genes$gene[window_indices]




# 4. EXTRACT AND NORMALISE COVERAGE -------------------------------------------

sample_cols <- colnames(diploid_coverage_genes)[7:ncol(diploid_coverage_genes)]

# Normalise each sample so median non-chr31 depth = 2
non_chr31 <- diploid_coverage_genes %>%
  filter(chrom != "LinJ.31") %>%
  select(all_of(sample_cols))

sample_medians <- apply(non_chr31, 2, median, na.rm = TRUE)

normalised <- diploid_coverage_genes %>%
  mutate(across(all_of(sample_cols),
                ~ . / sample_medians[cur_column()] * 2))

# Filter to target samples and genes of interest
target_samples <- non_americas$sra_run_accession

plot_data <- normalised %>%
  filter(gene %in% genes_of_interest) %>%
  select(gene, all_of(target_samples)) %>%
  pivot_longer(cols = all_of(target_samples),
               names_to  = "sample",
               values_to = "normalised_depth") %>%
  left_join(non_americas %>%
              select(sra_run_accession, source_country, coverage_rounded),
            by = c("sample" = "sra_run_accession")) %>%
  mutate(
    gene = factor(gene, levels = genes_of_interest),
    gene_type = case_when(
      gene %in% msl_genes ~ "MSL gene",
      TRUE                ~ "Flanking gene"
    ),
    sample_label = paste0(source_country, " (", sample, ", MSL", coverage_rounded, ")")
  )




# 5. PLOT FUNCTION ------------------------------------------------------------

# Horizontal reference lines at expected copy number levels:
#   0 = full MSL deletion, 2 = diploid baseline, 4 = wild-type tetraploid
reference_lines <- list(
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4),
  geom_hline(yintercept = 2, linetype = "dashed", colour = "grey40", linewidth = 0.4),
  geom_hline(yintercept = 4, linetype = "dashed", colour = "grey40", linewidth = 0.4)
)

make_single_plot <- function(data, strain_label) {
  ggplot(data, aes(x = gene, y = normalised_depth,
                   colour = gene_type, group = 1)) +
    reference_lines +
    geom_point(size = 3) +
    geom_line(colour = "grey70", linewidth = 0.5) +
    scale_colour_manual(values = c("MSL gene"      = "#BB5566",
                                   "Flanking gene"  = "#4477AA"),
                        name = "Gene type") +
    scale_y_continuous(breaks = seq(0, 10, by = 2)) +
    theme_minimal() +
    theme(
      axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
      panel.grid.major.x = element_blank(),
      legend.key.size    = unit(0.8, "cm"),
      plot.title         = element_text(hjust = 0.5, size = 13)
    ) +
    labs(
      title = strain_label,
      x     = "Gene",
      y     = "Normalised read depth (median non-chr31 = 2)"
    )
}




# 6. GENERATE AND SAVE INDIVIDUAL PLOTS ---------------------------------------

FIGURES_DIR <- file.path(BASE_DIR, "Figures/MSL Verification Non-Americas")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)

sample_labels_ordered <- unique(plot_data$sample_label)

for (lbl in sample_labels_ordered) {
  strain_data <- plot_data %>% filter(sample_label == lbl)
  p <- make_single_plot(strain_data, lbl)
  
  # Clean filename from label
  safe_name <- gsub("[^a-zA-Z0-9_]", "_", lbl)
  safe_name <- gsub("_+", "_", safe_name)
  
  ggsave(file.path(FIGURES_DIR, paste0(safe_name, ".jpeg")),
         plot = p, width = 10, height = 6, dpi = 300)
  
  ggsave(file.path(FIGURES_DIR, paste0(safe_name, ".pdf")),
         plot = p, width = 10, height = 6)
  
  cat("Saved:", safe_name, "\n")
}