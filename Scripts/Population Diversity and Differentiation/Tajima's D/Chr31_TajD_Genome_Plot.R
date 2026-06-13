# Chr31 Tajima's D Genome-Wide Boxplot
# Compares Tajima's D across all chromosomes for MSL0 vs MSL4 strains
# Goal: test whether chr31 is an outlier in MSL0 strains specifically



# 1. Load packages -------------------------------------------------------------

library(tidyverse)
library(ggplot2)




# 2. Paths ---------------------------------------------------------------------

data_dir <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Population Diversity and Differentiation Analysis/Americas_wide"

output_dir <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Population Diversity and Differentiation"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)




# 3. Load data --------------------------------------------------------------------

# Non-chr31 chromosomes (diploid calls)
msl0_genome <- read_tsv(file.path(data_dir, "Americas_only_MSL0_samples.tsv.Tajima.D"),
                        show_col_types = FALSE) %>%
  mutate(MSL_status = "MSL0")

msl4_genome <- read_tsv(file.path(data_dir, "Americas_only_MSL4_samples.tsv.Tajima.D"),
                        show_col_types = FALSE) %>%
  mutate(MSL_status = "MSL4")

# Chr31 only (diploidised tetraploid calls)
msl0_chr31 <- read_tsv(file.path(data_dir, "Americas_only_MSL0_samples.diploidised.tsv.Tajima.D"),
                       show_col_types = FALSE) %>%
  mutate(MSL_status = "MSL0")

msl4_chr31 <- read_tsv(file.path(data_dir, "Americas_only_MSL4_samples.diploidised.tsv.Tajima.D"),
                       show_col_types = FALSE) %>%
  mutate(MSL_status = "MSL4")

# Combine all
all_data <- bind_rows(msl0_genome, msl4_genome, msl0_chr31, msl4_chr31)

cat("Total windows loaded:", nrow(all_data), "\n")
cat("Chromosomes present:", paste(sort(unique(all_data$CHROM)), collapse = ", "), "\n")




# 4. Prepare data -----------------------------------------------------------------

# Remove windows with very few SNPs (unreliable estimates) and NA chromosomes
all_data <- all_data %>%
  filter(N_SNPS >= 10, !is.na(CHROM))

# Clean chromosome names for display (LinJ.01 -> 01, etc.)
all_data <- all_data %>%
  mutate(
    chrom_clean = str_replace(CHROM, "LinJ\\.", ""),
    chrom_clean = factor(chrom_clean, levels = c(
      sprintf("%02d", 1:30), "31", "32", "33", "34", "35", "36"
    )),
    is_chr31 = CHROM == "LinJ.31",
    MSL_status = factor(MSL_status, levels = c("MSL0", "MSL4"))
  ) %>%
  filter(!is.na(chrom_clean))




# 5. Paul Tol Palette ------------------------------------------------------------

# MSL0 = red (deletion), MSL4 = blue (wild-type) — Paul Tol palette
msl_colours <- c("MSL0" = "#BB5566", "MSL4" = "#4477AA")

# Chr31 highlight colour
chr31_highlight <- "#FFDD44"




# 6. Plot ----------------------------------------------------------------------

p <- ggplot(all_data, aes(x = chrom_clean, y = TajimaD, fill = MSL_status)) +
  
  # Highlight chr31 panel background
  annotate("rect",
           xmin = which(levels(all_data$chrom_clean) == "31") - 0.5,
           xmax = which(levels(all_data$chrom_clean) == "31") + 0.5,
           ymin = -Inf, ymax = Inf,
           fill = chr31_highlight, alpha = 0.3) +
  
  # Box and whisker plots
  geom_boxplot(
    outlier.size = 0.5,
    outlier.alpha = 0.4,
    linewidth = 0.4,
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  
  # Colour scale
  scale_fill_manual(
    values = msl_colours,
    labels = c("MSL0 (deletion)", "MSL4 (wild-type)"),
    name = "MSL status"
  ) +
  
  # Reference line at 0
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  
  # Labels
  labs(
    title = "Tajima's D across all chromosomes: MSL0 vs MSL4",
    subtitle = "Chr31 highlighted — test for chromosome-specific selection signal in deletion strains",
    x = "Chromosome",
    y = "Tajima's D (50kb windows)",
    caption = "MSL0: complete deletion strains; MSL4: wild-type strains. Windows with <10 SNPs excluded.\nChr31 uses diploidised tetraploid calls; all other chromosomes use diploid calls."
  ) +
  
  # Theme
  theme_bw(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9, colour = "grey30"),
    plot.caption = element_text(size = 7, colour = "grey50")
  )

print(p)




# 7. Save ----------------------------------------------------------------------

ggsave(file.path(output_dir, "chr31_tajd_genome_boxplot.pdf"),
       plot = p, width = 14, height = 6, dpi = 300)

ggsave(file.path(output_dir, "chr31_tajd_genome_boxplot.jpeg"),
       plot = p, width = 14, height = 6, dpi = 300)

cat("Saved to:", output_dir, "\n")




# 8. Summary stats — median Tajima's D per chromosome per MSL status -----------

summary_stats <- all_data %>%
  group_by(CHROM, MSL_status) %>%
  summarise(
    n_windows = n(),
    median_TajimaD = median(TajimaD, na.rm = TRUE),
    mean_TajimaD = mean(TajimaD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(CHROM, MSL_status)

cat("\nSummary — median Tajima's D per chromosome:\n")
print(summary_stats, n = Inf)

