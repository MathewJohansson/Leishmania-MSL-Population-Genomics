# SLIDING WINDOW FST — CHR31
# Populations: AM2 (pop2) and AM5 (pop5) MSL0 vs MSL4
# Input: VCFtools 1kb windowed Fst files (already computed)
# Method: aggregate 1kb windows into 50kb bins (weighted by N_VARIANTS)
# Output: per-population panels + combined figure




## 1. LIBRARIES ----------------------------------------------------------------

library(tidyverse)
library(ggpubr)




## 2. PATHS --------------------------------------------------------------------

BASE_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
DATA_DIR <- file.path(BASE_DIR, "Data/Population Diversity and Differentiation")
OUT_DIR  <- file.path(BASE_DIR, "Figures/Population Diversity and Differentiation/Fst")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# MSL locus coordinates
MSL_START <- 1181282
MSL_END   <- 1192406

# Window size for aggregation
BIN_SIZE  <- 50000

# Neutral expectation — mean weighted Fst across non-chr31 chromosomes
# Replace NA with computed values once genome-wide Fst data downloaded from Viking
neutral_fst_am2 <- NA
neutral_fst_am5 <- NA




## 3. LOAD FST FILES -----------------------------------------------------------

fst_files <- list(
  AM2 = file.path(DATA_DIR, "fst_2__Americas.K13.pop2.MSL0.diploidised.sample_names.txt__Americas.K13.pop2.MSL4.diploidised.sample_names.txt__tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.diploidised.vcf.gz.windowed.weir.fst"),
  AM5 = file.path(DATA_DIR, "fst_5__Americas.K13.pop5.MSL0.diploidised.sample_names.txt__Americas.K13.pop5.MSL4.diploidised.sample_names.txt__tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.diploidised.vcf.gz.windowed.weir.fst")
)

load_fst <- function(path, pop_label) {
  read_tsv(path, show_col_types = FALSE) |>
    mutate(
      population = pop_label,
      WEIGHTED_FST = pmax(WEIGHTED_FST, 0, na.rm = TRUE)  # floor at 0
    )
}

fst_raw <- bind_rows(
  load_fst(fst_files$AM2, "AM2"),
  load_fst(fst_files$AM5, "AM5")
)

print(paste("Total 1kb windows loaded:", nrow(fst_raw)))




## 4. AGGREGATE INTO 50KB BINS -------------------------------------------------

aggregate_to_bins <- function(df, bin_size = 50000) {
  df |>
    mutate(bin = floor((BIN_START - 1) / bin_size) * bin_size + 1) |>
    group_by(CHROM, population, bin) |>
    summarise(
      bin_start   = min(BIN_START),
      bin_end     = max(BIN_END),
      bin_mid     = (min(BIN_START) + max(BIN_END)) / 2,
      n_variants  = sum(N_VARIANTS, na.rm = TRUE),
      weighted_fst = if (sum(N_VARIANTS, na.rm = TRUE) > 0)
        sum(WEIGHTED_FST * N_VARIANTS, na.rm = TRUE) / sum(N_VARIANTS, na.rm = TRUE)
      else NA_real_,
      .groups = "drop"
    ) |>
    filter(n_variants >= 3)  # require at least 3 SNPs per 50kb bin
}

fst_binned <- aggregate_to_bins(fst_raw, BIN_SIZE)

print(paste("Total 50kb bins after aggregation:", nrow(fst_binned)))
print(paste("AM2 bins:", sum(fst_binned$population == "AM2")))
print(paste("AM5 bins:", sum(fst_binned$population == "AM5")))




## 5. PLOT FUNCTION ------------------------------------------------------------

plot_fst_window <- function(data, title, colour, y_max = NULL, neutral_fst = NULL,
                            x_max = NULL, legend_pos = c(0.92, 0.88)) {
  
  if (is.null(y_max)) y_max <- max(data$weighted_fst, na.rm = TRUE) * 1.1
  if (is.null(x_max)) x_max <- max(data$bin_mid, na.rm = TRUE) / 1e6 + 0.1
  
  p <- ggplot(data, aes(x = bin_mid / 1e6, y = weighted_fst)) +
    # MSL locus highlight
    geom_vline(xintercept = (MSL_START + MSL_END) / 2 / 1e6,
               colour = "#BB5566", linetype = "dashed", linewidth = 0.8) +
    annotate("text",
             x = (MSL_START + MSL_END) / 2 / 1e6 + 0.03,
             y = y_max * 0.95,
             label = "MSL", colour = "#BB5566", size = 3.5, fontface = "bold",
             hjust = 0) +
    geom_line(colour = colour, linewidth = 0.6) +
    geom_point(aes(size = n_variants), colour = colour, alpha = 0.7) +
    scale_size_continuous(name = "SNPs in bin", range = c(2, 8)) +
    scale_x_continuous(name = "Chr31 position (Mb)",
                       breaks = seq(0, 2.5, by = 0.5),
                       limits = c(0, x_max)) +
    scale_y_continuous(name = "Weighted Fst",
                       limits = c(0, y_max)) +
    labs(title = title) +
    theme_classic(base_size = 12) +
    theme(
      plot.title   = element_text(face = "bold", size = 13),
      legend.position  = legend_pos,
      legend.background = element_rect(fill = alpha("white", 0.7), colour = NA),
      legend.key.size  = unit(0.4, "cm"),
      legend.title     = element_text(size = 9),
      legend.text      = element_text(size = 8)
    )
  
  # Add neutral expectation line once genome-wide Fst data available
  if (!is.null(neutral_fst) && !is.na(neutral_fst)) {
    p <- p +
      geom_hline(yintercept = neutral_fst, linetype = "dashed",
                 colour = "grey40", linewidth = 0.7) +
      annotate("text", x = 0.05, y = neutral_fst + y_max * 0.03,
               label = "Genome-wide mean Fst", colour = "grey40",
               size = 3, hjust = 0)
  }
  
  p
}




## 6. INDIVIDUAL POPULATION PLOTS ----------------------------------------------

# Shared y-axis max across both populations for comparability
y_max_shared <- max(fst_binned$weighted_fst, na.rm = TRUE) * 1.1

# Individual plots — x axis auto-fits to data extent
plot_am2 <- plot_fst_window(
  data        = fst_binned |> filter(population == "AM2"),
  title       = "AM2: MSL0 vs MSL4 — Fst across Chr31 (50kb windows)",
  colour      = "#EE6677",  # AM2 Paul Tol palette colour
  y_max       = y_max_shared,
  neutral_fst = neutral_fst_am2
)

plot_am5 <- plot_fst_window(
  data        = fst_binned |> filter(population == "AM5"),
  title       = "AM5: MSL0 vs MSL4 — Fst across Chr31 (50kb windows)",
  colour      = "#66CCEE",  # AM5 Paul Tol palette colour
  y_max       = y_max_shared,
  neutral_fst = neutral_fst_am5
)

# Combined panel versions — fixed x axis so both panels align
plot_am2_combined <- plot_fst_window(
  data        = fst_binned |> filter(population == "AM2"),
  title       = "AM2: MSL0 vs MSL4 — Fst across Chr31 (50kb windows)",
  colour      = "#EE6677",
  y_max       = y_max_shared,
  legend_pos  = c(0.92, 0.65),
  neutral_fst = neutral_fst_am2
)

plot_am5_combined <- plot_fst_window(
  data        = fst_binned |> filter(population == "AM5"),
  title       = "AM5: MSL0 vs MSL4 — Fst across Chr31 (50kb windows)",
  colour      = "#66CCEE",
  y_max       = y_max_shared,
  legend_pos  = c(0.92, 0.65),
  neutral_fst = neutral_fst_am5
)

ggsave(file.path(OUT_DIR, "fst_sliding_window_chr31_AM2.jpeg"),
       plot = plot_am2, width = 14, height = 5, dpi = 300)

ggsave(file.path(OUT_DIR, "fst_sliding_window_chr31_AM5.jpeg"),
       plot = plot_am5, width = 14, height = 5, dpi = 300)




## 7. NEUTRAL EXPECTATION — GENOME-WIDE FST BASELINE -------------------------
# TODO: Download genome-wide Fst files from Viking for all non-chr31 chromosomes.
# These should be the equivalent windowed Fst files for AM2 and AM5 across all
# other chromosomes (LinJ.01 to LinJ.30, LinJ.32 to LinJ.36).
# Once downloaded, compute the mean weighted Fst across all non-chr31 windows
# and add as a horizontal dashed line to both the AM2 and AM5 plots.
#
# Steps:
#   1. Load all non-chr31 .windowed.weir.fst files for AM2 and AM5
#   2. Compute mean weighted Fst (weighted by N_VARIANTS) across all windows
#      separately for AM2 and AM5
#   3. Add geom_hline(yintercept = neutral_fst_am2, linetype = "dashed") to plot_am2
#      and equivalent for plot_am5
#
# REQUIRES: Genome-wide Fst files on Viking —
#   /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/

# Placeholder values — replace once Viking data downloaded:
# neutral_fst_am2 and neutral_fst_am5 are defined in the config block (section 2).




## 8. COMBINED PANEL (AM2 + AM5) -----------------------------------------------

combined_panel <- ggarrange(
  plot_am2_combined, plot_am5_combined,
  ncol = 1, nrow = 2,
  labels = c("A", "B"),
  font.label = list(size = 14, face = "bold")
)

ggsave(file.path(OUT_DIR, "fst_sliding_window_chr31_combined.jpeg"),
       plot = combined_panel, width = 14, height = 10, dpi = 300)

print("All sliding window Fst figures saved.")