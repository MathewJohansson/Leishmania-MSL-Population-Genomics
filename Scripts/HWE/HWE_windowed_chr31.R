# HWE WINDOWED ANALYSIS — CHROMOSOME 31 (TETRAPLOID)
# Sliding window HWE across chr31, for AM and LI subsets separately.
# Per-SNP HWE p-values and odds ratios, summarised per window.
# Also runs pooled (all-samples) analysis alongside per-population.
#
# Dan's pseudocode implemented:
#   - Subset SNPs within each window
#   - Per SNP: HWE p-value + odds ratio
#   - Per window: median p, lowest p, median OR, lowest OR,
#     proportion of 1111 and 0000 genotypes relative to HWE expectation
#
# Input VCF: tetraploid chr31 filtered VCF (tetraploid calls).
# Outliers removed: ERR205787, SRR8608748.
# Population assignments from master membership TSVs (>=50% threshold, consistent with HWE_FINAL.R).
#
# DEV NOTE: Set DEV_MODE <- TRUE to test on a small chr31 region first (Dan's suggestion).
#           Set to FALSE for the full chromosome run.




# 1. LIBRARIES -----------------------------------------------------------------

library(vcfR)
library(tidyverse)
library(patchwork)




# 2. CONFIG --------------------------------------------------------------------

BASE_DIR   <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
METADATA_DIR <- file.path(BASE_DIR, "Metadata")
DATA_DIR     <- file.path(BASE_DIR, "Data")
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/HWE/Windowed Chr31")
TABLES_DIR   <- file.path(BASE_DIR, "Tables/HWE/Windowed Chr31")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)

# VCF — tetraploid chr31 calls
VCF_PATH <- file.path(DATA_DIR,
                      "tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.vcf.gz")

# Membership TSVs
LI_MEMBERSHIP <- file.path(DATA_DIR, "population_structure",
                           "population_membership_table_infantum_only_k13_LI_notation.tsv")
AM_MEMBERSHIP <- file.path(DATA_DIR, "population_structure",
                           "population_membership_table_Americas_only_k13_AM_notation.tsv")

# Outliers to remove from LI and AM analyses
OUTLIERS <- c("ERR205787", "SRR8608748")

# Window parameters
WINDOW_SIZE <- 10000   # 10kb
STEP_SIZE   <- 5000    # 5kb step (50% overlap)

# Population assignment threshold (consistent with HWE_FINAL.R)
POP_THRESHOLD <- 0.50

# DEV MODE: restrict to a small region for testing (set FALSE for full run)
DEV_MODE       <- FALSE
DEV_REGION_START <- 1100000
DEV_REGION_END   <- 1300000

cat("=== HWE WINDOWED CHR31 ANALYSIS ===\n")
cat("Window size:", WINDOW_SIZE, "bp | Step:", STEP_SIZE, "bp\n")
cat("Dev mode:", DEV_MODE, "\n\n")




# 3. LOAD VCF ------------------------------------------------------------------

cat("Loading VCF...\n")
vcf <- read.vcfR(VCF_PATH, verbose = FALSE)
cat("VCF loaded:", nrow(vcf@fix), "variants,", ncol(vcf@gt) - 1, "samples\n")

# Extract positions
chrom <- getCHROM(vcf)
pos   <- as.numeric(getPOS(vcf))

# Restrict to chr31 (should already be chr31-only, but be explicit)
chr31_mask <- chrom == "LinJ.31"
vcf <- vcf[chr31_mask, ]
pos <- pos[chr31_mask]
cat("After chr31 filter:", nrow(vcf@fix), "variants\n")

# Dev mode: restrict to small region
if (DEV_MODE) {
  dev_mask <- pos >= DEV_REGION_START & pos <= DEV_REGION_END
  vcf  <- vcf[dev_mask, ]
  pos  <- pos[dev_mask]
  cat("DEV MODE: restricted to", DEV_REGION_START, "-", DEV_REGION_END,
      "->", nrow(vcf@fix), "variants\n")
}




# 4. EXTRACT GENOTYPE MATRIX ---------------------------------------------------

# Extract raw GT strings (e.g. "0/0/1/1") and convert to alt allele dosage (0-4)
gt_strings <- extract.gt(vcf, element = "GT")

count_alt_alleles <- function(x) {
  lengths(regmatches(x, gregexpr("1", x)))
}

dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings),
  ncol = ncol(gt_strings)
)
rownames(dosage_matrix) <- rownames(gt_strings)   # SNP IDs
colnames(dosage_matrix) <- colnames(gt_strings)   # Sample IDs

# Restore NAs
dosage_matrix[is.na(gt_strings) |
                gt_strings == "./." |
                gt_strings == "./././."] <- NA

# Transpose: samples as rows, SNPs as columns (consistent with other scripts)
dosage_matrix <- t(dosage_matrix)   # [samples x SNPs]

cat("Dosage matrix:", nrow(dosage_matrix), "samples x", ncol(dosage_matrix), "SNPs\n")




# 5. LOAD MEMBERSHIP AND BUILD SUBSET SAMPLE LISTS ----------------------------

li_mem <- read_tsv(LI_MEMBERSHIP, show_col_types = FALSE) %>%
  filter(population_proportion >= POP_THRESHOLD) %>%
  filter(!sample %in% OUTLIERS)

am_mem <- read_tsv(AM_MEMBERSHIP, show_col_types = FALSE) %>%
  filter(population_proportion >= POP_THRESHOLD) %>%
  filter(!sample %in% OUTLIERS)

# Samples present in the VCF
vcf_samples <- rownames(dosage_matrix)

li_samples <- li_mem$sample[li_mem$sample %in% vcf_samples]
am_samples <- am_mem$sample[am_mem$sample %in% vcf_samples]

cat("LI samples in VCF (post-filter):", length(li_samples), "\n")
cat("AM samples in VCF (post-filter):", length(am_samples), "\n\n")




# 6. CORE FUNCTIONS ------------------------------------------------------------

# -- 6a. Per-SNP HWE test (tetraploid) ----------------------------------------
# Takes a numeric vector of dosages (0-4) for one SNP across samples.
# Returns: p_value, chi_square, df, odds_ratio (1111:0000 observed/expected),
#          prop_1111_obs, prop_0000_obs, prop_1111_exp, prop_0000_exp, n_samples

snp_hwe_tetraploid <- function(dosages) {
  dosages <- dosages[!is.na(dosages)]
  n <- length(dosages)
  
  if (n < 5) {
    return(list(p_value = NA, log_pval = NA, chi_square = NA, df = NA,
                odds_ratio = NA,
                prop_1111_obs = NA, prop_0000_obs = NA,
                prop_1111_exp = NA, prop_0000_exp = NA,
                n_samples = n))
  }
  
  # Genotype counts (dosage 0-4 maps to 0000,1000,1100,1110,1111)
  counts <- tabulate(dosages + 1, nbins = 5)  # index 1=dosage0 ... 5=dosage4
  names(counts) <- c("0000", "1000", "1100", "1110", "1111")
  
  # Allele frequency of alt (deletion) allele
  total_alleles <- n * 4
  alt_alleles   <- sum(dosages)
  q <- alt_alleles / total_alleles   # freq of alt
  p <- 1 - q                         # freq of ref
  
  # Clamp to avoid numerical issues
  p <- max(0.001, min(0.999, p))
  q <- 1 - p
  
  # Expected frequencies under tetraploid HWE (polysomic)
  exp_freq <- c(
    "0000" = q^4,
    "1000" = 4 * p * q^3,
    "1100" = 6 * p^2 * q^2,
    "1110" = 4 * p^3 * q,
    "1111" = p^4
  )
  exp_counts <- exp_freq * n
  
  # Only test classes with expected count >= 1
  valid <- exp_counts >= 1.0
  n_valid <- sum(valid)
  
  if (n_valid < 2) {
    return(list(p_value = NA, log_pval = NA, chi_square = NA, df = NA,
                odds_ratio = NA,
                prop_1111_obs = counts["1111"] / n,
                prop_0000_obs = counts["0000"] / n,
                prop_1111_exp = exp_freq["1111"],
                prop_0000_exp = exp_freq["0000"],
                n_samples = n))
  }
  
  obs_v <- counts[valid]
  exp_v <- exp_counts[valid]
  chi_sq <- sum((obs_v - exp_v)^2 / exp_v)
  df     <- max(1, n_valid - 2)
  pval     <- pchisq(chi_sq, df, lower.tail = FALSE)
  log_pval <- pchisq(chi_sq, df, lower.tail = FALSE, log.p = TRUE)  # true log p-value, avoids 1e-16 floor
  
  # Odds ratio: (obs 1111 / exp 1111) / (obs 0000 / exp 0000)
  # NA if fewer than 5 samples have non-reference dosage (unstable ratio)
  # Guards against division by zero
  obs_1111 <- counts["1111"];  exp_1111 <- exp_counts["1111"]
  obs_0000 <- counts["0000"];  exp_0000 <- exp_counts["0000"]
  n_alt <- sum(dosages > 0, na.rm = TRUE)
  
  if (n_alt >= 5 & exp_1111 > 0 & exp_0000 > 0 & obs_0000 > 0) {
    or <- (obs_1111 / exp_1111) / (obs_0000 / exp_0000)
  } else {
    or <- NA
  }
  
  list(
    p_value      = pval,
    log_pval     = log_pval,
    chi_square   = chi_sq,
    df           = df,
    odds_ratio   = or,
    prop_1111_obs = counts["1111"] / n,
    prop_0000_obs = counts["0000"] / n,
    prop_1111_exp = exp_freq["1111"],
    prop_0000_exp = exp_freq["0000"],
    n_samples    = n
  )
}


# -- 6b. Windowed summary across SNPs -----------------------------------------
# dm:   dosage matrix [samples x SNPs] (subset to relevant samples already)
# pos:  numeric vector of SNP positions (length = ncol(dm))
# Returns a data.frame with one row per window.

run_windowed_hwe <- function(dm, positions, window_size, step_size, label = "") {
  
  chr_min <- min(positions)
  chr_max <- max(positions)
  
  window_starts <- seq(chr_min, chr_max - window_size + 1, by = step_size)
  n_windows <- length(window_starts)
  
  cat(label, "- running", n_windows, "windows...\n")
  
  results <- vector("list", n_windows)
  
  for (i in seq_along(window_starts)) {
    w_start <- window_starts[i]
    w_end   <- w_start + window_size - 1
    w_mid   <- (w_start + w_end) / 2
    
    snp_idx <- which(positions >= w_start & positions <= w_end)
    n_snps  <- length(snp_idx)
    
    if (n_snps == 0) {
      results[[i]] <- data.frame(
        window_start = w_start, window_end = w_end, window_mid = w_mid,
        n_snps = 0,
        median_chisq = NA, median_OR = NA, median_neg_log10p = NA,
        median_prop_1111_obs = NA, median_prop_0000_obs = NA,
        median_prop_1111_exp = NA, median_prop_0000_exp = NA,
        mean_n_samples = NA
      )
      next
    }
    
    # Run per-SNP HWE across this window
    snp_results <- apply(dm[, snp_idx, drop = FALSE], 2, snp_hwe_tetraploid)
    
    pvals    <- sapply(snp_results, `[[`, "p_value")
    log_pvals <- sapply(snp_results, `[[`, "log_pval")
    chisqs  <- sapply(snp_results, `[[`, "chi_square")
    ors     <- sapply(snp_results, `[[`, "odds_ratio")
    p1111o  <- sapply(snp_results, `[[`, "prop_1111_obs")
    p0000o  <- sapply(snp_results, `[[`, "prop_0000_obs")
    p1111e  <- sapply(snp_results, `[[`, "prop_1111_exp")
    p0000e  <- sapply(snp_results, `[[`, "prop_0000_exp")
    ns      <- sapply(snp_results, `[[`, "n_samples")
    
    results[[i]] <- data.frame(
      window_start         = w_start,
      window_end           = w_end,
      window_mid           = w_mid,
      n_snps               = n_snps,
      median_chisq         = median(chisqs,   na.rm = TRUE),
      median_OR            = median(ors,      na.rm = TRUE),
      median_neg_log10p    = median(-log_pvals / log(10), na.rm = TRUE),  # true -log10(p) via log.p=TRUE
      median_prop_1111_obs = median(p1111o,   na.rm = TRUE),
      median_prop_0000_obs = median(p0000o,   na.rm = TRUE),
      median_prop_1111_exp = median(p1111e,   na.rm = TRUE),
      median_prop_0000_exp = median(p0000e,   na.rm = TRUE),
      mean_n_samples       = mean(ns, na.rm = TRUE)
    )
  }
  
  bind_rows(results)
}




# 7. RUN WINDOWED ANALYSIS — POOLED AND PER-POPULATION ------------------------

run_all <- function(subset_label, sample_ids, membership_df, positions) {
  
  cat("\n===", subset_label, "===\n")
  
  dm_sub <- dosage_matrix[rownames(dosage_matrix) %in% sample_ids, , drop = FALSE]
  cat("Samples in matrix:", nrow(dm_sub), "\n")
  
  # -- 7a. POOLED (all samples together) --
  cat("Running POOLED analysis...\n")
  pooled_result <- run_windowed_hwe(
    dm_sub, positions, WINDOW_SIZE, STEP_SIZE,
    label = paste(subset_label, "pooled")
  )
  pooled_result$population <- "POOLED"
  pooled_result$subset      <- subset_label
  
  # -- 7b. PER-POPULATION --
  populations <- unique(membership_df$population)
  pop_results <- list()
  
  for (pop in sort(populations)) {
    pop_samples <- membership_df$sample[membership_df$population == pop]
    pop_samples_in_vcf <- intersect(pop_samples, rownames(dm_sub))
    
    if (length(pop_samples_in_vcf) < 3) {
      cat("Skipping", pop, "(< 3 samples in VCF)\n")
      next
    }
    
    dm_pop <- dm_sub[pop_samples_in_vcf, , drop = FALSE]
    pop_res <- run_windowed_hwe(
      dm_pop, positions, WINDOW_SIZE, STEP_SIZE,
      label = paste(subset_label, pop)
    )
    pop_res$population <- pop
    pop_res$subset      <- subset_label
    pop_results[[pop]]  <- pop_res
  }
  
  bind_rows(c(list(pooled_result), pop_results))
}

# SNP positions aligned to current dosage_matrix columns
snp_positions <- pos   # already subset to DEV region if DEV_MODE

li_results <- run_all(
  "LI",
  li_samples,
  li_mem %>% filter(sample %in% li_samples),
  snp_positions
)

am_results <- run_all(
  "AM",
  am_samples,
  am_mem %>% filter(sample %in% am_samples),
  snp_positions
)

all_results <- bind_rows(li_results, am_results)




# 8. SAVE TABLES ---------------------------------------------------------------

suffix <- if (DEV_MODE) "_DEV" else ""

write_tsv(li_results, file.path(TABLES_DIR, paste0("windowed_HWE_LI", suffix, ".tsv")))
write_tsv(am_results, file.path(TABLES_DIR, paste0("windowed_HWE_AM", suffix, ".tsv")))
write_tsv(all_results, file.path(TABLES_DIR, paste0("windowed_HWE_all", suffix, ".tsv")))

cat("\nTables saved to", TABLES_DIR, "\n")




# 9. PLOTTING ------------------------------------------------------------------

# MSL region boundaries (excluded from NJ trees but shown as shading here)
MSL_START <- 1181282
MSL_END   <- 1192406

# Paul Tol palette — MSL0 = RED, MSL4 = BLUE (not needed for these plots,
# but population colours follow ADMIXTURE palette from memory)
# Pooled shown in dark grey; per-population in colour where relevant.

plot_windowed_hwe <- function(results_df, subset_label, metric, ylab,
                              show_msl_shade = TRUE) {
  
  pooled <- results_df %>% filter(population == "POOLED")
  pops   <- results_df %>% filter(population != "POOLED")
  
  # Build base plot from pooled
  p <- ggplot(pooled, aes(x = window_mid / 1e6)) +
    geom_line(aes(y = .data[[metric]]), colour = "#444444", linewidth = 0.7) +
    labs(
      title = paste(subset_label, "—", ylab, "across chr31"),
      subtitle = paste0("Window: ", WINDOW_SIZE / 1000, "kb, Step: ",
                        STEP_SIZE / 1000, "kb"),
      x = "Chr31 position (Mb)",
      y = ylab
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.minor = element_blank())
  
  # MSL region shading (only meaningful on full run — shown regardless)
  if (show_msl_shade) {
    p <- p +
      annotate("rect",
               xmin = MSL_START / 1e6, xmax = MSL_END / 1e6,
               ymin = -Inf, ymax = Inf,
               fill = "#BB5566", alpha = 0.15) +
      annotate("text",
               x = (MSL_START + MSL_END) / 2 / 1e6,
               y = Inf, vjust = 1.5, size = 3,
               label = "MSL", colour = "#BB5566")
  }
  
  p
}


# -- Generate plots for pooled results ----------------------------------------

plot_list <- list()

for (ss in c("LI", "AM")) {
  res <- all_results %>% filter(subset == ss)
  
  # Median chi-square: HWE deviation per window. Trough at MSL = low variation inside deletion.
  plot_list[[paste0(ss, "_chisq")]] <- plot_windowed_hwe(
    res, ss, "median_chisq", "Median chi-square per window"
  )
  # True median -log10(p): extracted via log.p=TRUE to avoid R's 1e-16 floor.
  # Trough at MSL = lower HWE deviation than chr31 background.
  plot_list[[paste0(ss, "_neg_log10p")]] <- plot_windowed_hwe(
    res, ss, "median_neg_log10p", "Median -log10(p) per window"
  )
  # Median odds ratio: (obs_1111/exp_1111) / (obs_0000/exp_0000).
  # Peak at MSL in LI = MSL0/MSL4 coexistence inflating genotype extremes.
  plot_list[[paste0(ss, "_OR")]] <- plot_windowed_hwe(
    res, ss, "median_OR", "Median odds ratio per window"
  )
}

# -- Save plots ----------------------------------------------------------------

for (nm in names(plot_list)) {
  fname <- file.path(FIGURES_DIR, paste0(nm, suffix, ".png"))
  ggsave(fname, plot_list[[nm]], width = 12, height = 4, dpi = 200)
  cat("Saved:", basename(fname), "\n")
}

# Combined 3-panel (chi-square + -log10p + OR) for quick visual check
for (ss in c("LI", "AM")) {
  combined <- plot_list[[paste0(ss, "_chisq")]] /
    plot_list[[paste0(ss, "_neg_log10p")]] /
    plot_list[[paste0(ss, "_OR")]]
  fname <- file.path(FIGURES_DIR, paste0(ss, "_combined_check", suffix, ".png"))
  ggsave(fname, combined, width = 12, height = 12, dpi = 200)
  cat("Saved:", basename(fname), "\n")
}

# -- Per-population chi-square plots -------------------------------------------
# Plots median chi-square per window for each population separately.
# Uses Paul Tol ADMIXTURE palette — consistent with all other manuscript figures.
# One faceted figure per subset (LI and AM), all populations overlaid as separate
# lines, plus pooled shown in dark grey for reference.

paul_tol_palette <- c(
  "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB",
  "#EE7733", "#009988", "#994F00", "#332288", "#9955BB", "#FFAABB"
)

plot_per_population <- function(results_df, subset_label, mem_df) {
  
  # Sort populations numerically (e.g. AM1, AM2, ... AM13) not lexicographically
  populations <- unique(mem_df$population)
  populations <- populations[order(as.integer(gsub("[^0-9]", "", populations)))]
  n_pops      <- length(populations)
  
  # Assign Paul Tol colours by population order
  pop_colours <- setNames(paul_tol_palette[seq_len(n_pops)], populations)
  
  pooled_data <- results_df %>%
    filter(population == "POOLED") %>%
    mutate(pop_label = "POOLED")
  
  pop_data <- results_df %>%
    filter(population %in% populations) %>%
    mutate(
      pop_label  = population,
      population = factor(population, levels = populations)  # numerical order for facets
    )
  
  # -- Panel 1: all populations overlaid on one plot --
  p_overlay <- ggplot() +
    # MSL shading
    annotate("rect",
             xmin = MSL_START / 1e6, xmax = MSL_END / 1e6,
             ymin = -Inf, ymax = Inf,
             fill = "#BB5566", alpha = 0.15) +
    annotate("text",
             x = (MSL_START + MSL_END) / 2 / 1e6,
             y = Inf, vjust = 1.5, size = 3,
             label = "MSL", colour = "#BB5566") +
    # Pooled in grey
    geom_line(data = pooled_data,
              aes(x = window_mid / 1e6, y = median_chisq),
              colour = "#444444", linewidth = 0.5, linetype = "dashed") +
    # Per-population lines
    geom_line(data = pop_data,
              aes(x = window_mid / 1e6, y = median_chisq,
                  colour = population),
              linewidth = 0.5, alpha = 0.8) +
    scale_colour_manual(values = pop_colours, name = "Population") +
    labs(
      title    = paste(subset_label, "— Median chi-square per window, by population"),
      subtitle = paste0("Window: ", WINDOW_SIZE / 1000, "kb, Step: ",
                        STEP_SIZE / 1000, "kb | Dashed grey = pooled"),
      x = "Chr31 position (Mb)",
      y = "Median chi-square"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "right",
          panel.grid.minor = element_blank())
  
  # -- Panel 2: faceted by population --
  p_facet <- ggplot(pop_data,
                    aes(x = window_mid / 1e6, y = median_chisq,
                        colour = population)) +
    annotate("rect",
             xmin = MSL_START / 1e6, xmax = MSL_END / 1e6,
             ymin = -Inf, ymax = Inf,
             fill = "#BB5566", alpha = 0.15) +
    geom_line(linewidth = 0.5) +
    scale_colour_manual(values = pop_colours, guide = "none") +
    facet_wrap(~ population, ncol = 3, scales = "free_y") +
    labs(
      title    = paste(subset_label, "— Per-population chi-square (faceted)"),
      subtitle = paste0("Window: ", WINDOW_SIZE / 1000, "kb, Step: ",
                        STEP_SIZE / 1000, "kb"),
      x = "Chr31 position (Mb)",
      y = "Median chi-square"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"))
  
  list(overlay = p_overlay, facet = p_facet)
}

cat("\nGenerating per-population plots...\n")
# Note: OR not plotted per-population — the median OR is dominated by single-SNP
# artefacts in small populations (extreme obs/exp ratio where expected counts near zero).
# Chi-square is the informative metric at population level.

for (ss in c("LI", "AM")) {
  mem <- if (ss == "LI") li_mem else am_mem
  res <- all_results %>% filter(subset == ss)
  
  plots <- plot_per_population(res, ss, mem)
  
  # Save overlay
  fname_overlay <- file.path(FIGURES_DIR, paste0(ss, "_perpop_overlay", suffix, ".png"))
  ggsave(fname_overlay, plots$overlay, width = 14, height = 5, dpi = 300)
  cat("Saved:", basename(fname_overlay), "\n")
  
  # Save faceted — height scales with number of populations
  n_pops     <- length(unique(mem$population))
  fig_height <- ceiling(n_pops / 3) * 3.5
  fname_facet <- file.path(FIGURES_DIR, paste0(ss, "_perpop_facet", suffix, ".png"))
  ggsave(fname_facet, plots$facet, width = 14, height = fig_height, dpi = 300)
  cat("Saved:", basename(fname_facet), "\n")
}

cat("DEV_MODE was:", DEV_MODE, "\n")
cat("If output looks sensible, set DEV_MODE <- FALSE and rerun for full chr31.\n")