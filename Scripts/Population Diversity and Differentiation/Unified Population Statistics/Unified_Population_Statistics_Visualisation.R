
# UNIFIED POPULATION STATISTICS — VISUALISATION

# Reads the TSV outputs from:
#   Unified_Population_Statistics_Chr31.R    → within_pop_chr31.tsv

# CURRENT SCOPE: π only (Tajima's D, Dxy, Fst to be added once π is verified).

# Produces two figures:
#   A1. π per group — boxplot, one box per eligible group, log10 y-axis
#   A2. π along chr31 — windowed line plot, one line per group, MSL locus marked

# COLOUR CONVENTIONS (project-wide — do not alter):
#   MSL0 = #BB5566 (red)
#   MSL4 = #4477AA (blue)




# 1. LIBRARIES -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})




# 2. PATHS ---------------------------------------------------------------------

BASE_DIR   <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
DATA_DIR   <- file.path(BASE_DIR, "Data/Unified_Popgen")
FIG_DIR    <- file.path(BASE_DIR, "Figures/Unified_Popgen")

dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# MSL locus position (for windowed plots, added later)
MSL_START <- 1181282
MSL_END   <- 1192406
MSL_MID   <- (MSL_START + MSL_END) / 2e6   # in Mb, for ggplot x-axis




# 3. COLOUR AND THEME DEFINITIONS ----------------------------------------------

# Project colour conventions
COL_MSL0 <- "#BB5566"   # red — MSL0 and MSL0-MSL0 comparisons
COL_MSL4 <- "#4477AA"   # blue — MSL4 and MSL4-MSL4 comparisons

# Named colour vector for msl_type factor levels
msl_type_cols <- c("MSL0-MSL0" = COL_MSL0, "MSL4-MSL4" = COL_MSL4)

# Base ggplot theme used across all figures
theme_unified <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title       = element_text(face = "bold", size = base_size + 1),
      plot.subtitle    = element_text(size = base_size - 1, colour = "grey40"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      legend.position  = "bottom",
      legend.title     = element_text(face = "bold"),
      strip.background = element_blank(),
      strip.text       = element_text(face = "bold")
    )
}

# Helper: determine MSL type from group label (ends in _0 or _4)
msl_type_from_group <- function(g) {
  ifelse(grepl("_0$", g), "MSL0", "MSL4")
}




# 4. LOAD DATA -----------------------------------------------------------------

cat("Loading chr31 within-population data...\n")

within_chr31 <- read_tsv(
  file.path(DATA_DIR, "within_pop_chr31.tsv"),
  show_col_types = FALSE
) |>
  mutate(
    msl_type      = ifelse(grepl("_0$", group), "MSL0", "MSL4"),
    window_mid_mb = window_mid / 1e6
  )

cat(sprintf("Rows: %d | Groups: %d | Windows: %d\n",
            nrow(within_chr31),
            n_distinct(within_chr31$group),
            n_distinct(within_chr31$window_mid)))

# Group order: by population number, then MSL type within each number
group_order <- within_chr31 |>
  distinct(group, msl_type) |>
  mutate(num = as.integer(str_extract(group, "[0-9]+"))) |>
  arrange(num, msl_type) |>
  pull(group)

within_chr31 <- within_chr31 |>
  mutate(group = factor(group, levels = group_order))




# 5. WILCOXON TEST — MSL0 vs MSL4 (pooled across all groups and windows) -------
# Pooled test: all windowed π values from MSL0 groups vs all from MSL4 groups.

wilcox_pi <- wilcox.test(pi ~ msl_type, data = within_chr31, exact = FALSE)

cat(sprintf("\nWilcoxon π  MSL0 vs MSL4: W = %.0f, p = %.3e\n",
            wilcox_pi$statistic, wilcox_pi$p.value))

fmt_p <- function(p) {
  if (p < 0.001) sprintf("p = %.2e", p) else sprintf("p = %.3f", p)
}




# =============================================================================
# SECTION A — π FIGURES
# =============================================================================

# A1. π per group — boxplot ----------------------------------------------------
# One box per eligible group. Log10 y-axis to match Dan's L figures.
# Jitter removed — boxes alone are cleaner at this scale.
# Wilcoxon p annotated at top.

p_pi_boxplot <- ggplot(within_chr31, aes(x = group, y = pi, fill = msl_type)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4,
               alpha = 0.8, width = 0.6) +
  scale_fill_manual(
    values = c("MSL0" = COL_MSL0, "MSL4" = COL_MSL4),
    name   = "MSL type"
  ) +
  scale_y_log10(
    labels = scales::label_scientific()
  ) +
  annotate(
    "text",
    x     = length(group_order) / 2,
    y     = max(within_chr31$pi, na.rm = TRUE) * 1.5,
    label = fmt_p(wilcox_pi$p.value),
    size  = 3.5,
    colour = "grey30"
  ) +
  labs(
    title    = "Nucleotide diversity (π) per population — Chr31",
    subtitle = "50kb sliding windows, 5kb step; MSL locus excluded; log10 y-axis",
    x        = "Population group",
    y        = "π (nucleotide diversity)"
  ) +
  theme_unified()

ggsave(file.path(FIG_DIR, "pi_per_group_chr31.jpeg"),
       p_pi_boxplot, width = 14, height = 6, dpi = 300)
cat("Saved: pi_per_group_chr31.jpeg\n")


# A2. π along chr31 — windowed line plot ---------------------------------------
# Mean π per MSL type per window — two lines only (MSL0 red, MSL4 blue).
# Matches Dan's L3 figure panel B approach.
# Individual group lines removed — per-population breakdown shown in boxplot.

pi_windowed_mean <- within_chr31 |>
  group_by(msl_type, window_mid, window_mid_mb) |>
  summarise(pi_mean = mean(pi, na.rm = TRUE), .groups = "drop")

p_pi_windowed <- ggplot(pi_windowed_mean,
                        aes(x = window_mid_mb, y = pi_mean,
                            colour = msl_type, group = msl_type)) +
  geom_vline(xintercept = MSL_MID,
             linetype = "solid", colour = "grey60", linewidth = 0.6) +
  annotate("text", x = MSL_MID + 0.02, y = Inf,
           label = "MSL", vjust = 1.5, hjust = 0,
           size = 3, colour = "grey40") +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  scale_colour_manual(
    values = c("MSL0" = COL_MSL0, "MSL4" = COL_MSL4),
    name   = "MSL type"
  ) +
  scale_y_log10(labels = scales::label_scientific()) +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  labs(
    title    = "Windowed π along Chr31",
    subtitle = "Mean per MSL type; 50kb windows, 5kb step; MSL locus (grey line) excluded",
    x        = "Position on Chr31",
    y        = "π (nucleotide diversity)"
  ) +
  theme_unified() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(FIG_DIR, "pi_windowed_chr31.jpeg"),
       p_pi_windowed, width = 14, height = 6, dpi = 300)
cat("Saved: pi_windowed_chr31.jpeg\n")


# A3. Combined: boxplot (left) + windowed line plot (right) --------------------

p_pi_combined <- p_pi_boxplot + p_pi_windowed +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title      = "Chromosome 31 Nucleotide Diversity (π) — Americas",
    subtitle   = "50kb sliding windows; MSL locus excluded",
    tag_levels = "A",
    theme      = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

ggsave(file.path(FIG_DIR, "pi_combined_chr31.jpeg"),
       p_pi_combined, width = 18, height = 6, dpi = 300)
cat("Saved: pi_combined_chr31.jpeg\n")


cat("\nAll π figures saved to:", FIG_DIR, "\n")


# =============================================================================
# SECTION B — TAJIMA'S D FIGURES
# =============================================================================

# B1. Tajima's D per group — boxplot ------------------------------------------

wilcox_td <- wilcox.test(tajima_d ~ msl_type, data = within_chr31, exact = FALSE)
cat(sprintf("\nWilcoxon TajD  MSL0 vs MSL4: W = %.0f, p = %.3e\n",
            wilcox_td$statistic, wilcox_td$p.value))

p_td_boxplot <- ggplot(within_chr31, aes(x = group, y = tajima_d, fill = msl_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4, alpha = 0.8, width = 0.6) +
  scale_fill_manual(values = c("MSL0" = COL_MSL0, "MSL4" = COL_MSL4),
                    name   = "MSL type") +
  annotate("text",
           x     = length(group_order) / 2,
           y     = max(within_chr31$tajima_d, na.rm = TRUE) * 1.05,
           label = fmt_p(wilcox_td$p.value),
           size  = 3.5, colour = "grey30") +
  labs(
    title    = "Tajima's D per population — Chr31",
    subtitle = "50kb sliding windows, 5kb step; MSL locus excluded; dashed line = 0",
    x        = "Population group",
    y        = "Tajima's D"
  ) +
  theme_unified()

ggsave(file.path(FIG_DIR, "tajima_d_per_group_chr31.jpeg"),
       p_td_boxplot, width = 14, height = 6, dpi = 300)
cat("Saved: tajima_d_per_group_chr31.jpeg\n")


# B2. Tajima's D along chr31 — windowed line plot (mean per MSL type) ----------

td_windowed_mean <- within_chr31 |>
  group_by(msl_type, window_mid, window_mid_mb) |>
  summarise(td_mean = mean(tajima_d, na.rm = TRUE), .groups = "drop")

p_td_windowed <- ggplot(td_windowed_mean,
                        aes(x = window_mid_mb, y = td_mean,
                            colour = msl_type, group = msl_type)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = MSL_MID,
             linetype = "solid", colour = "grey60", linewidth = 0.6) +
  annotate("text", x = MSL_MID + 0.02, y = Inf,
           label = "MSL", vjust = 1.5, hjust = 0,
           size = 3, colour = "grey40") +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  scale_colour_manual(values = c("MSL0" = COL_MSL0, "MSL4" = COL_MSL4),
                      name   = "MSL type") +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  labs(
    title    = "Windowed Tajima's D along Chr31",
    subtitle = "Mean per MSL type; 50kb windows, 5kb step; MSL locus (grey line) excluded",
    x        = "Position on Chr31",
    y        = "Tajima's D"
  ) +
  theme_unified() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(FIG_DIR, "tajima_d_windowed_chr31.jpeg"),
       p_td_windowed, width = 14, height = 6, dpi = 300)
cat("Saved: tajima_d_windowed_chr31.jpeg\n")


# B3. Combined: Tajima's D boxplot + windowed line plot ------------------------

p_td_combined <- p_td_boxplot + p_td_windowed +
  plot_layout(widths = c(1, 2)) +
  plot_annotation(
    title      = "Chromosome 31 Tajima's D — Americas",
    subtitle   = "50kb sliding windows; MSL locus excluded",
    tag_levels = "A",
    theme      = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

ggsave(file.path(FIG_DIR, "tajima_d_combined_chr31.jpeg"),
       p_td_combined, width = 18, height = 6, dpi = 300)
cat("Saved: tajima_d_combined_chr31.jpeg\n")

cat("\nAll Tajima's D figures saved to:", FIG_DIR, "\n")


# =============================================================================
# SECTION C — DXY FIGURES
# =============================================================================

between_chr31 <- read_tsv(
  file.path(DATA_DIR, "between_pop_chr31.tsv"),
  show_col_types = FALSE
) |>
  mutate(
    msl_type      = factor(msl_type, levels = c("MSL0-MSL0", "MSL4-MSL4")),
    window_mid_mb = window_mid / 1e6
  )

cat(sprintf("Between-pop chr31: %d rows, %d pairs\n",
            nrow(between_chr31),
            n_distinct(paste(between_chr31$pop1, between_chr31$pop2))))

msl_type_cols <- c("MSL0-MSL0" = COL_MSL0, "MSL4-MSL4" = COL_MSL4)

wilcox_dxy <- wilcox.test(dxy ~ msl_type, data = between_chr31, exact = FALSE)
cat(sprintf("\nWilcoxon Dxy  MSL0-MSL0 vs MSL4-MSL4: W = %.0f, p = %.3e\n",
            wilcox_dxy$statistic, wilcox_dxy$p.value))


# C1. Dxy — MSL0-MSL0 vs MSL4-MSL4 boxplot ------------------------------------

p_dxy_boxplot <- ggplot(between_chr31, aes(x = msl_type, y = dxy, fill = msl_type)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4, alpha = 0.8, width = 0.5) +
  scale_fill_manual(values = msl_type_cols, guide = "none") +
  annotate("text", x = 1.5,
           y = quantile(between_chr31$dxy, 0.99, na.rm = TRUE) * 1.05,
           label = fmt_p(wilcox_dxy$p.value),
           size = 3.5, colour = "grey30") +
  labs(
    title    = "Dxy by comparison type — Chr31",
    subtitle = "MSL0-MSL0 vs MSL4-MSL4 pairs; 50kb windows, 5kb step",
    x        = "Comparison type",
    y        = "Dxy (mean absolute divergence)"
  ) +
  theme_unified() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(FIG_DIR, "dxy_boxplot_chr31.jpeg"),
       p_dxy_boxplot, width = 6, height = 6, dpi = 300)
cat("Saved: dxy_boxplot_chr31.jpeg\n")


# C2. Dxy along chr31 — windowed line plot (mean per MSL type) -----------------

dxy_windowed_mean <- between_chr31 |>
  group_by(msl_type, window_mid, window_mid_mb) |>
  summarise(dxy_mean = mean(dxy, na.rm = TRUE), .groups = "drop")

p_dxy_windowed <- ggplot(dxy_windowed_mean,
                         aes(x = window_mid_mb, y = dxy_mean,
                             colour = msl_type, group = msl_type)) +
  geom_vline(xintercept = MSL_MID,
             linetype = "solid", colour = "grey60", linewidth = 0.6) +
  annotate("text", x = MSL_MID + 0.02, y = Inf,
           label = "MSL", vjust = 1.5, hjust = 0,
           size = 3, colour = "grey40") +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  scale_colour_manual(values = msl_type_cols, name = "Comparison") +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  labs(
    title    = "Windowed Dxy along Chr31",
    subtitle = "Mean per comparison type; 50kb windows, 5kb step; MSL locus (grey line) excluded",
    x        = "Position on Chr31",
    y        = "Dxy (mean absolute divergence)"
  ) +
  theme_unified() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(FIG_DIR, "dxy_windowed_chr31.jpeg"),
       p_dxy_windowed, width = 14, height = 6, dpi = 300)
cat("Saved: dxy_windowed_chr31.jpeg\n")


# C3. Combined: Dxy boxplot + windowed line plot --------------------------------

p_dxy_combined <- p_dxy_boxplot + p_dxy_windowed +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(
    title      = "Chromosome 31 Absolute Divergence (Dxy) — Americas",
    subtitle   = "50kb sliding windows; MSL locus excluded",
    tag_levels = "A",
    theme      = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

ggsave(file.path(FIG_DIR, "dxy_combined_chr31.jpeg"),
       p_dxy_combined, width = 18, height = 6, dpi = 300)
cat("Saved: dxy_combined_chr31.jpeg\n")

cat("\nAll Dxy figures saved to:", FIG_DIR, "\n")


# =============================================================================
# SECTION D — FST FIGURES
# =============================================================================

wilcox_fst <- wilcox.test(fst ~ msl_type, data = between_chr31, exact = FALSE)
cat(sprintf("\nWilcoxon Fst  MSL0-MSL0 vs MSL4-MSL4: W = %.0f, p = %.3e\n",
            wilcox_fst$statistic, wilcox_fst$p.value))


# D1. Fst — MSL0-MSL0 vs MSL4-MSL4 boxplot ------------------------------------

p_fst_boxplot <- ggplot(between_chr31, aes(x = msl_type, y = fst, fill = msl_type)) +
  geom_boxplot(outlier.size = 0.8, outlier.alpha = 0.4, alpha = 0.8, width = 0.5) +
  scale_fill_manual(values = msl_type_cols, guide = "none") +
  annotate("text", x = 1.5,
           y = quantile(between_chr31$fst, 0.99, na.rm = TRUE) * 1.05,
           label = fmt_p(wilcox_fst$p.value),
           size = 3.5, colour = "grey30") +
  labs(
    title    = "Fst by comparison type — Chr31",
    subtitle = "MSL0-MSL0 vs MSL4-MSL4 pairs; 50kb windows, 5kb step",
    x        = "Comparison type",
    y        = "Fst (Hudson's estimator)"
  ) +
  theme_unified() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(FIG_DIR, "fst_boxplot_chr31.jpeg"),
       p_fst_boxplot, width = 6, height = 6, dpi = 300)
cat("Saved: fst_boxplot_chr31.jpeg\n")


# D2. Fst along chr31 — windowed line plot (mean per MSL type) -----------------

fst_windowed_mean <- between_chr31 |>
  group_by(msl_type, window_mid, window_mid_mb) |>
  summarise(fst_mean = mean(fst, na.rm = TRUE), .groups = "drop")

p_fst_windowed <- ggplot(fst_windowed_mean,
                         aes(x = window_mid_mb, y = fst_mean,
                             colour = msl_type, group = msl_type)) +
  geom_vline(xintercept = MSL_MID,
             linetype = "solid", colour = "grey60", linewidth = 0.6) +
  annotate("text", x = MSL_MID + 0.02, y = Inf,
           label = "MSL", vjust = 1.5, hjust = 0,
           size = 3, colour = "grey40") +
  geom_line(linewidth = 0.7, alpha = 0.9) +
  scale_colour_manual(values = msl_type_cols, name = "Comparison") +
  scale_x_continuous(labels = scales::label_number(suffix = " Mb")) +
  labs(
    title    = "Windowed Fst along Chr31",
    subtitle = "Mean per comparison type; 50kb windows, 5kb step; MSL locus (grey line) excluded",
    x        = "Position on Chr31",
    y        = "Fst (Hudson's estimator)"
  ) +
  theme_unified() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

ggsave(file.path(FIG_DIR, "fst_windowed_chr31.jpeg"),
       p_fst_windowed, width = 14, height = 6, dpi = 300)
cat("Saved: fst_windowed_chr31.jpeg\n")


# D3. Combined: Fst boxplot + windowed line plot --------------------------------

p_fst_combined <- p_fst_boxplot + p_fst_windowed +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(
    title      = "Chromosome 31 Fixation Index (Fst) — Americas",
    subtitle   = "50kb sliding windows; MSL locus excluded; min 5 samples per group per SNP",
    tag_levels = "A",
    theme      = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )

ggsave(file.path(FIG_DIR, "fst_combined_chr31.jpeg"),
       p_fst_combined, width = 18, height = 6, dpi = 300)
cat("Saved: fst_combined_chr31.jpeg\n")

cat("\nAll Fst figures saved to:", FIG_DIR, "\n")


# =============================================================================
# SECTION E — PLACEHOLDER: Chr31 vs genome-wide comparison
# =============================================================================

# Add once NonChr31 TSVs are correct (π fix applied and rerun complete).