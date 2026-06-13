# TETRAPLOID DXY / PI / FST — WINDOWED CHR31

# Input:  dosage_matrix (samples x SNPs, values 0-4, NA for missing)
#         pos           (SNP positions along chr31)
#         membership TSV (population assignments)

# Output: per-window TSV of pi, Fst, Dxy for a given population pair

# Test run: single window first, then 10 windows, then full chr31.




# 1. LIBRARIES -----------------------------------------------------------------

library(tidyverse)
library(vcfR)




# 2. PATHS ---------------------------------------------------------------------

BASE_DIR  <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
VCF_PATH  <- file.path(BASE_DIR, "Data/tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.vcf.gz")
META_PATH <- file.path(BASE_DIR, "Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
OUT_DIR   <- file.path(BASE_DIR, "Data/Pegas_Dxy")
FIG_DIR   <- file.path(BASE_DIR, "Figures/Population Diversity and Differentiation/Pegas_Dxy")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# MSL locus to exclude from all calculations
MSL_START <- 1181282
MSL_END   <- 1192406

# Window parameters (matching HWE windowed analysis)
WIN_SIZE  <- 50000L   # 50kb
WIN_STEP  <- 25000L   # 25kb step




# 3. LOAD VCF AND BUILD DOSAGE MATRIX ------------------------------------------

cat("Loading VCF...\n")
vcf <- read.vcfR(VCF_PATH)

chrom <- getCHROM(vcf)
pos   <- as.numeric(getPOS(vcf))

# Exclude MSL locus
keep_snps  <- !(chrom == "LinJ.31" & pos >= MSL_START & pos <= MSL_END)
vcf_no_msl <- vcf[keep_snps, ]
pos_no_msl <- pos[keep_snps]

# Extract dosage (0-4 alt allele count)
count_alt_alleles <- function(x) lengths(regmatches(x, gregexpr("1", x)))

gt_strings    <- extract.gt(vcf_no_msl, element = "GT")
dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings), ncol = ncol(gt_strings)
)
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)
dosage_matrix[is.na(gt_strings) | gt_strings %in% c("./.", "./././.")] <- NA
dosage_matrix <- t(dosage_matrix)   # now: samples x SNPs

cat(sprintf("Dosage matrix: %d samples x %d SNPs\n", nrow(dosage_matrix), ncol(dosage_matrix)))




# 4. TETRAPLOID POPULATION GENETIC FUNCTIONS -----------------------------------

#' Tetraploid nucleotide diversity (pi)
#' For each SNP, pi = (n / (n-1)) * 2 * p * (1-p) where p = mean dosage / 4
#' Mean taken across SNPs in window.
#' @param dos_window  matrix: samples x SNPs (values 0-4, NA allowed)
#' @return scalar pi estimate
compute_pi_tet <- function(dos_window) {
  if (is.null(dos_window) || ncol(dos_window) == 0) return(NA_real_)
  per_snp <- apply(dos_window, 2, function(x) {
    x   <- x[!is.na(x)]
    n   <- length(x)
    if (n < 2) return(NA_real_)
    p   <- mean(x) / 4          # alt allele frequency
    (n / (n - 1)) * 2 * p * (1 - p)
  })
  mean(per_snp, na.rm = TRUE)
}

#' Tetraploid Dxy — mean absolute pairwise divergence between two populations
#' For each SNP: mean of all cross-population pairwise |dosage_i - dosage_j| / 4
#' Mean taken across SNPs in window.
#' @param dos_window  matrix: samples x SNPs (values 0-4, NA allowed)
#' @param ids1        sample IDs for population 1
#' @param ids2        sample IDs for population 2
#' @return scalar Dxy estimate
compute_dxy_tet <- function(dos_window, ids1, ids2) {
  p1 <- intersect(ids1, rownames(dos_window))
  p2 <- intersect(ids2, rownames(dos_window))
  if (length(p1) < 1 || length(p2) < 1) return(NA_real_)
  
  mat1 <- dos_window[p1, , drop = FALSE]
  mat2 <- dos_window[p2, , drop = FALSE]
  
  per_snp <- sapply(seq_len(ncol(dos_window)), function(j) {
    a <- mat1[, j]; b <- mat2[, j]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    # All pairwise |a_i - b_j| / 4, averaged
    pairs <- outer(a, b, function(x, y) abs(x - y) / 4)
    mean(pairs, na.rm = TRUE)
  })
  mean(per_snp, na.rm = TRUE)
}

#' Tetraploid Fst (Weir & Cockerham approximation via allele frequencies)
#' Uses Hudson's estimator: Fst = (pi_between - pi_within) / pi_between
#' @param dos_window  matrix: samples x SNPs (values 0-4, NA allowed)
#' @param ids1        sample IDs for population 1
#' @param ids2        sample IDs for population 2
#' @return scalar Fst estimate (clamped to [0, 1])
compute_fst_tet <- function(dos_window, ids1, ids2) {
  p1 <- intersect(ids1, rownames(dos_window))
  p2 <- intersect(ids2, rownames(dos_window))
  if (length(p1) < 2 || length(p2) < 2) return(NA_real_)
  
  mat1 <- dos_window[p1, , drop = FALSE]
  mat2 <- dos_window[p2, , drop = FALSE]
  
  per_snp <- sapply(seq_len(ncol(dos_window)), function(j) {
    a <- mat1[, j] / 4; b <- mat2[, j] / 4   # convert to [0,1] freq
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) < 2 || length(b) < 2) return(NA_real_)
    
    p1_freq <- mean(a); p2_freq <- mean(b)
    n1 <- length(a);    n2 <- length(b)
    
    # Hudson Fst estimator
    pi1     <- 2 * p1_freq * (1 - p1_freq) * n1 / (n1 - 1)
    pi2     <- 2 * p2_freq * (1 - p2_freq) * n2 / (n2 - 1)
    pi_w    <- (pi1 + pi2) / 2
    pi_b    <- p1_freq * (1 - p2_freq) + p2_freq * (1 - p1_freq)
    if (pi_b == 0) return(NA_real_)
    (pi_b - pi_w) / pi_b
  })
  
  fst <- mean(per_snp, na.rm = TRUE)
  max(0, min(1, fst))   # clamp to [0,1]
}




# 5. WINDOWED ANALYSIS FUNCTION ------------------------------------------------

#' Run windowed pi / Dxy / Fst along chr31 for a given population pair.
#' @param ids1       sample IDs for population 1
#' @param ids2       sample IDs for population 2
#' @param label1     label for population 1 (for output)
#' @param label2     label for population 2 (for output)
#' @param test_n     if not NULL, only run on first test_n windows (for testing)
#' @return data.frame with one row per window
run_windowed_metrics <- function(ids1, ids2, label1, label2, test_n = NULL) {
  
  cat(sprintf("\nRunning: %s (n=%d) vs %s (n=%d)\n",
              label1, length(ids1), label2, length(ids2)))
  
  chr_min <- min(pos_no_msl)
  chr_max <- max(pos_no_msl)
  
  # Build window start positions
  starts <- seq(chr_min, chr_max - WIN_SIZE + 1, by = WIN_STEP)
  if (!is.null(test_n)) {
    starts <- head(starts, test_n)
    cat(sprintf("TEST MODE: running first %d windows only\n", test_n))
  }
  
  results <- vector("list", length(starts))
  
  for (i in seq_along(starts)) {
    w_start <- starts[i]
    w_end   <- w_start + WIN_SIZE - 1
    w_mid   <- (w_start + w_end) / 2
    
    # SNP indices in this window
    in_win  <- which(pos_no_msl >= w_start & pos_no_msl <= w_end)
    n_snps  <- length(in_win)
    
    if (n_snps == 0) {
      results[[i]] <- data.frame(
        pop1 = label1, pop2 = label2,
        window_start = w_start, window_end = w_end, window_mid = w_mid,
        n_snps = 0L, pi1 = NA, pi2 = NA, dxy = NA, fst = NA
      )
      next
    }
    
    dos_win <- dosage_matrix[, in_win, drop = FALSE]
    
    results[[i]] <- data.frame(
      pop1         = label1,
      pop2         = label2,
      window_start = w_start,
      window_end   = w_end,
      window_mid   = w_mid,
      n_snps       = n_snps,
      pi1          = compute_pi_tet(dos_win[intersect(ids1, rownames(dos_win)), , drop = FALSE]),
      pi2          = compute_pi_tet(dos_win[intersect(ids2, rownames(dos_win)), , drop = FALSE]),
      dxy          = compute_dxy_tet(dos_win, ids1, ids2),
      fst          = compute_fst_tet(dos_win, ids1, ids2)
    )
    
    if (i %% 50 == 0) cat(sprintf("  Window %d / %d\n", i, length(starts)))
  }
  
  bind_rows(results)
}




# 6. MSL0 / MSL4 CLASSIFICATION FROM METADATA TSV ------------------------------

# Uses the normalised MSL_coverage column from the sample metadata file.
# This is normalised coverage (copy-number scale, 0-4), not raw read counts.
# Thresholds: MSL0 = < 0.5; MSL4 = >= 3.5; intermediate = excluded.
# Source: vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv

cat("Loading MSL classification from metadata TSV...\n")
meta_full <- read_tsv(META_PATH, show_col_types = FALSE)

msl_class <- meta_full |>
  select(sample_id, MSL_coverage) |>
  mutate(MSL_coverage = as.numeric(MSL_coverage))

msl0_ids <- msl_class$sample_id[msl_class$MSL_coverage <  0.5]
msl4_ids <- msl_class$sample_id[msl_class$MSL_coverage >= 3.5]
msl_inter <- msl_class$sample_id[msl_class$MSL_coverage >= 0.5 & msl_class$MSL_coverage < 3.5]

cat(sprintf("MSL0 (coverage < 0.5):    n=%d\n", length(msl0_ids)))
cat(sprintf("MSL4 (coverage >= 3.5):   n=%d\n", length(msl4_ids)))
cat(sprintf("Intermediate (excluded):  n=%d\n", length(msl_inter)))




# 7. DEFINE POPULATION GROUPS --------------------------------------------------

# All groups are GEOGRAPHIC (state-level).
# WFR1 (PI/MA) is the reference region compared against all other qualifying regions.
# Qualifying regions: have >= 1 MSL0 AND >= 1 MSL4 sample.
# From metadata:
#   WFR1 (Piauí + Maranhão): MSL0=85, MSL4=75
#   WFR2 (Mato Grosso do Sul): MSL0=1, MSL4=11  
#   Minas Gerais:              MSL0=35, MSL4=2
#   Mato Grosso:               MSL0=16, MSL4=2
#   Santa Catarina:            MSL0=17, MSL4=2
#   Pará:                      MSL0=2,  MSL4=1
# META_PATH already loaded above as meta_full.

get_region_ids <- function(states) {
  meta_full$sample_id[meta_full$source_state_or_region %in% states &
                        !is.na(meta_full$source_state_or_region)]
}

# WFR1 — reference region
ids_WFR1_all  <- get_region_ids(c("Piauí", "Maranhão"))
ids_WFR1_MSL0 <- intersect(ids_WFR1_all, msl0_ids)
ids_WFR1_MSL4 <- intersect(ids_WFR1_all, msl4_ids)

# WFR2 — Mato Grosso do Sul
ids_WFR2_all  <- get_region_ids("Mato Grosso do Sul")
ids_WFR2_MSL0 <- intersect(ids_WFR2_all, msl0_ids)  
ids_WFR2_MSL4 <- intersect(ids_WFR2_all, msl4_ids)

# Minas Gerais
ids_MG_all    <- get_region_ids("Minas Gerais")
ids_MG_MSL0   <- intersect(ids_MG_all, msl0_ids)
ids_MG_MSL4   <- intersect(ids_MG_all, msl4_ids)

# Mato Grosso
ids_MT_all    <- get_region_ids("Mato Grosso")
ids_MT_MSL0   <- intersect(ids_MT_all, msl0_ids)
ids_MT_MSL4   <- intersect(ids_MT_all, msl4_ids)

# Santa Catarina
ids_SC_all    <- get_region_ids("Santa Catarina")
ids_SC_MSL0   <- intersect(ids_SC_all, msl0_ids)
ids_SC_MSL4   <- intersect(ids_SC_all, msl4_ids)

# Pará
ids_PA_all    <- get_region_ids("Pará")
ids_PA_MSL0   <- intersect(ids_PA_all, msl0_ids)
ids_PA_MSL4   <- intersect(ids_PA_all, msl4_ids)

# nonWFR geo — all Brazil excluding PI, MA, MS (WFR regions)
# Includes MG, MT, SC, PA, and all other Brazilian states.
# Dan confirmed this pooled comparison as a key test.
wfr_excl_states   <- c("Piauí", "Maranhão", "Mato Grosso do Sul")
ids_nonWFR_all    <- meta_full$sample_id[
  meta_full$source_country == "Brazil" &
    !meta_full$source_state_or_region %in% wfr_excl_states &
    !is.na(meta_full$source_state_or_region)
]
ids_nonWFR_MSL0   <- intersect(ids_nonWFR_all, msl0_ids)
ids_nonWFR_MSL4   <- intersect(ids_nonWFR_all, msl4_ids)

cat(sprintf("\nWFR1 (PI/MA):       total=%d  |  MSL0=%d  MSL4=%d\n",
            length(ids_WFR1_all), length(ids_WFR1_MSL0), length(ids_WFR1_MSL4)))
cat(sprintf("WFR2 (MS):          total=%d  |  MSL0=%d  MSL4=%d  [MSL0 n=1: Dxy runs, pi/Fst return NA]\n",
            length(ids_WFR2_all), length(ids_WFR2_MSL0), length(ids_WFR2_MSL4)))
cat(sprintf("Minas Gerais (MG):  total=%d  |  MSL0=%d  MSL4=%d\n",
            length(ids_MG_all), length(ids_MG_MSL0), length(ids_MG_MSL4)))
cat(sprintf("Mato Grosso (MT):   total=%d  |  MSL0=%d  MSL4=%d\n",
            length(ids_MT_all), length(ids_MT_MSL0), length(ids_MT_MSL4)))
cat(sprintf("Santa Catarina (SC):total=%d  |  MSL0=%d  MSL4=%d\n",
            length(ids_SC_all), length(ids_SC_MSL0), length(ids_SC_MSL4)))
cat(sprintf("Pará (PA):          total=%d  |  MSL0=%d  MSL4=%d\n",
            length(ids_PA_all), length(ids_PA_MSL0), length(ids_PA_MSL4)))
cat(sprintf("nonWFR geo (all BR excl PI/MA/MS): MSL0=%d  MSL4=%d\n",
            length(ids_nonWFR_MSL0), length(ids_nonWFR_MSL4)))




# 8. TEST RUN — SINGLE WINDOW --------------------------------------------------

# Tests WFR1_MSL0 vs MG_MSL0 — largest non-WFR group, key comparison.
# Check: pi ~0.001-0.01, dxy >= pi, fst in [0,1].

cat("\n--- SINGLE WINDOW TEST (WFR1_MSL0 vs MG_MSL0) ---\n")
test_1 <- run_windowed_metrics(ids_WFR1_MSL0, ids_MG_MSL0,
                               "WFR1_MSL0", "MG_MSL0",
                               test_n = 1)
print(test_1)
cat("Check pi1, pi2 ~0.001-0.01; dxy >= pi; fst in [0,1].\n\n")




# 9. FULL CHR31 RUN — ALL COMPARISONS -----------------------------------------

# Design: WFR1 (PI/MA) vs each qualifying region — paired MSL0-MSL0 and MSL4-MSL4.
# Qualifying regions: WFR2/MS, Minas Gerais, Mato Grosso, Santa Catarina, Pará.

# Comparison list:
#   C1:  WFR1_MSL0 vs WFR2_MSL0    — Dxy runs (n=1 sufficient); pi2 and Fst return NA (n<2)
#   C2:  WFR1_MSL4 vs WFR2_MSL4    — MSL4-MSL4 paired with C1
#   C3:  WFR1_MSL0 vs MG_MSL0      — key test
#   C4:  WFR1_MSL4 vs MG_MSL4      — MSL4-MSL4 paired with C3
#   C5:  WFR1_MSL0 vs MT_MSL0      — key test
#   C6:  WFR1_MSL4 vs MT_MSL4      — MSL4-MSL4 paired with C5
#   C7:  WFR1_MSL0 vs SC_MSL0      — key test
#   C8:  WFR1_MSL4 vs SC_MSL4      — MSL4-MSL4 paired with C7
#   C9:  WFR1_MSL0 vs PA_MSL0      — key test
#   C10: WFR1_MSL4 vs PA_MSL4      — MSL4-MSL4 paired with C9

# WFR2 — C1 Dxy runs with n=1 MSL0; pi2 and Fst return NA by design
comp1  <- run_windowed_metrics(ids_WFR1_MSL0, ids_WFR2_MSL0, "WFR1_MSL0", "WFR2_MSL0")
comp2  <- run_windowed_metrics(ids_WFR1_MSL4, ids_WFR2_MSL4, "WFR1_MSL4", "WFR2_MSL4")

# Minas Gerais
comp3  <- run_windowed_metrics(ids_WFR1_MSL0, ids_MG_MSL0,   "WFR1_MSL0", "MG_MSL0")
comp4  <- run_windowed_metrics(ids_WFR1_MSL4, ids_MG_MSL4,   "WFR1_MSL4", "MG_MSL4")

# Mato Grosso
comp5  <- run_windowed_metrics(ids_WFR1_MSL0, ids_MT_MSL0,   "WFR1_MSL0", "MT_MSL0")
comp6  <- run_windowed_metrics(ids_WFR1_MSL4, ids_MT_MSL4,   "WFR1_MSL4", "MT_MSL4")

# Santa Catarina
comp7  <- run_windowed_metrics(ids_WFR1_MSL0, ids_SC_MSL0,   "WFR1_MSL0", "SC_MSL0")
comp8  <- run_windowed_metrics(ids_WFR1_MSL4, ids_SC_MSL4,   "WFR1_MSL4", "SC_MSL4")

# Pará
comp9  <- run_windowed_metrics(ids_WFR1_MSL0, ids_PA_MSL0,   "WFR1_MSL0", "PA_MSL0")
comp10 <- run_windowed_metrics(ids_WFR1_MSL4, ids_PA_MSL4,   "WFR1_MSL4", "PA_MSL4")

# nonWFR geo — pooled (all Brazil excl PI, MA, MS)
# Dan confirmed this as a key comparison: WFR1 vs all non-WFR regions pooled.
#   C11: WFR1_MSL0 vs nonWFR_geo_MSL0  — key pooled test
#   C12: WFR1_MSL4 vs nonWFR_geo_MSL4  — MSL4-MSL4 paired control
comp11 <- run_windowed_metrics(ids_WFR1_MSL0, ids_nonWFR_MSL0, "WFR1_MSL0", "nonWFR_geo_MSL0")
comp12 <- run_windowed_metrics(ids_WFR1_MSL4, ids_nonWFR_MSL4, "WFR1_MSL4", "nonWFR_geo_MSL4")

all_comps <- bind_rows(comp1, comp2, comp3, comp4, comp5, comp6,
                       comp7, comp8, comp9, comp10, comp11, comp12)
write_tsv(all_comps, file.path(OUT_DIR, paste0("windowed_dxy_fst_chr31_.tsv")))
cat("Full chr31 results saved.\n")




# 10. SUMMARY — MEAN DXY BY REGION FOR EACH COMPARISON -------------------------

# Produces mean Dxy for background, MSL flanking, and MSL locus per comparison.
# Regions:
#   background   = all windows outside MSL locus +/- 50kb flanking zone
#   MSL_flanking = windows within 50kb either side of MSL locus
#   MSL_locus    = windows whose midpoint falls within 1,181,282-1,192,406

# Run after Section 10. Update filename to match today's TSV.

summarise_dxy <- function(df, label) {
  msl_start <- 1181282
  msl_end   <- 1192406
  flank_pad <- 50000
  
  df <- df |>
    mutate(region = case_when(
      window_mid >= msl_start & window_mid <= msl_end ~ "MSL_locus",
      window_mid >= (msl_start - flank_pad) & window_mid <= (msl_end + flank_pad) ~ "MSL_flanking",
      TRUE ~ "background"
    ))
  
  df |>
    group_by(region) |>
    summarise(
      n_windows = n(),
      mean_dxy  = round(mean(dxy, na.rm = TRUE), 6),
      .groups   = "drop"
    ) |>
    mutate(comparison = label) |>
    select(comparison, region, n_windows, mean_dxy)
}

# Load saved TSV from Section 10
all_comps <- read_tsv(file.path(OUT_DIR, paste0("windowed_dxy_fst_chr31_.tsv")),
                      show_col_types = FALSE)

all_summaries <- bind_rows(
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL0" & pop2 == "WFR2_MSL0"), "C1  — WFR1_MSL0 vs WFR2_MSL0 (key test; pi2/Fst NA — WFR2 MSL0 n=1)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL4" & pop2 == "WFR2_MSL4"), "C2  — WFR1_MSL4 vs WFR2_MSL4 (paired control)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL0" & pop2 == "MG_MSL0"),   "C3  — WFR1_MSL0 vs MG_MSL0 (key test)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL4" & pop2 == "MG_MSL4"),   "C4  — WFR1_MSL4 vs MG_MSL4 (paired control)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL0" & pop2 == "MT_MSL0"),   "C5  — WFR1_MSL0 vs MT_MSL0 (key test)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL4" & pop2 == "MT_MSL4"),   "C6  — WFR1_MSL4 vs MT_MSL4 (paired control)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL0" & pop2 == "SC_MSL0"),   "C7  — WFR1_MSL0 vs SC_MSL0 (key test)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL4" & pop2 == "SC_MSL4"),   "C8  — WFR1_MSL4 vs SC_MSL4 (paired control)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL0" & pop2 == "PA_MSL0"),   "C9  — WFR1_MSL0 vs PA_MSL0 (key test)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL4" & pop2 == "PA_MSL4"),   "C10 — WFR1_MSL4 vs PA_MSL4 (paired control)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL0" & pop2 == "nonWFR_geo_MSL0"), "C11 — WFR1_MSL0 vs nonWFR_geo_MSL0 (pooled key test)"),
  summarise_dxy(filter(all_comps, pop1 == "WFR1_MSL4" & pop2 == "nonWFR_geo_MSL4"), "C12 — WFR1_MSL4 vs nonWFR_geo_MSL4 (pooled paired control)")
)

print(all_summaries, n = Inf)
write_tsv(all_summaries, file.path(OUT_DIR, paste0("dxy_summary_by_region_.tsv")))
cat("Summary saved.\n")




# 11. WINDOWED DXY PLOTS -------------------------------------------------------

# Each region pair produces one figure with two lines:
#   MSL0-MSL0 (RED #BB5566) and MSL4-MSL4 (BLUE #4477AA) on the same panel.
# Five figures total (one per qualifying region).
# WFR2 MSL0-MSL0 line shows Dxy only — pi2 and Fst are NA (n=1) and not plotted.

library(ggplot2)
library(patchwork)

# Reload TSV (same file as section 11)
all_comps <- read_tsv(file.path(OUT_DIR, paste0("windowed_dxy_fst_chr31_.tsv")),
                      show_col_types = FALSE)

MSL_MID   <- (MSL_START + MSL_END) / 2 / 1e6
COL_MSL0  <- "#BB5566"   # RED
COL_MSL4  <- "#4477AA"   # BLUE


# Helper: paired plot — MSL0-MSL0 and MSL4-MSL4 on same panel
plot_paired_dxy <- function(data, msl0_pop1, msl0_pop2, msl4_pop1, msl4_pop2,
                            region_label, include_msl0 = TRUE) {
  
  d4 <- data |> filter(pop1 == msl4_pop1, pop2 == msl4_pop2) |>
    mutate(MSL_type = "MSL4-MSL4")
  
  y_max <- max(d4$dxy, na.rm = TRUE) * 1.1
  
  if (include_msl0) {
    d0 <- data |> filter(pop1 == msl0_pop1, pop2 == msl0_pop2) |>
      mutate(MSL_type = "MSL0-MSL0")
    y_max <- max(c(d0$dxy, d4$dxy), na.rm = TRUE) * 1.1
    plot_data <- bind_rows(d0, d4)
  } else {
    plot_data <- d4
  }
  
  ggplot(plot_data, aes(x = window_mid / 1e6, y = dxy,
                        colour = MSL_type, group = MSL_type)) +
    geom_vline(xintercept = MSL_MID, colour = "grey40",
               linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = MSL_MID + 0.02, y = y_max * 0.95,
             label = "MSL", colour = "grey40", size = 3,
             fontface = "bold", hjust = 0) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 1.5, alpha = 0.7) +
    scale_colour_manual(
      values = c("MSL0-MSL0" = COL_MSL0, "MSL4-MSL4" = COL_MSL4),
      name   = NULL
    ) +
    scale_x_continuous(name = "Chr31 position (Mb)") +
    scale_y_continuous(name = "Dxy", limits = c(0, y_max)) +
    labs(title = region_label) +
    theme_classic(base_size = 11) +
    theme(
      plot.title   = element_text(size = 10, face = "bold"),
      legend.position = "bottom"
    )
}


# WFR1 vs Minas Gerais 
p_MG <- plot_paired_dxy(all_comps,
                        "WFR1_MSL0", "MG_MSL0",
                        "WFR1_MSL4", "MG_MSL4",
                        "WFR1 vs Minas Gerais")
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_WFR1_vs_MG_.jpeg")),
       p_MG, width = 12, height = 5, dpi = 300)
cat("WFR1 vs MG plot saved.\n")


# WFR1 vs Mato Grosso 
p_MT <- plot_paired_dxy(all_comps,
                        "WFR1_MSL0", "MT_MSL0",
                        "WFR1_MSL4", "MT_MSL4",
                        "WFR1 vs Mato Grosso")
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_WFR1_vs_MT_.jpeg")),
       p_MT, width = 12, height = 5, dpi = 300)
cat("WFR1 vs MT plot saved.\n")


# WFR1 vs Santa Catarina
p_SC <- plot_paired_dxy(all_comps,
                        "WFR1_MSL0", "SC_MSL0",
                        "WFR1_MSL4", "SC_MSL4",
                        "WFR1 vs Santa Catarina")
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_WFR1_vs_SC_.jpeg")),
       p_SC, width = 12, height = 5, dpi = 300)
cat("WFR1 vs SC plot saved.\n")


# WFR1 vs Pará
p_PA <- plot_paired_dxy(all_comps,
                        "WFR1_MSL0", "PA_MSL0",
                        "WFR1_MSL4", "PA_MSL4",
                        "WFR1 vs Pará")
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_WFR1_vs_PA_.jpeg")),
       p_PA, width = 12, height = 5, dpi = 300)
cat("WFR1 vs PA plot saved.\n")


# WFR1 vs WFR2 — MSL0-MSL0 Dxy runs (n=1); pi2/Fst NA by design
p_WFR2 <- plot_paired_dxy(all_comps,
                          "WFR1_MSL0", "WFR2_MSL0",
                          "WFR1_MSL4", "WFR2_MSL4",
                          "WFR1 vs WFR2/MS (MSL0-MSL0 Dxy only — WFR2 MSL0 n=1)")
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_WFR1_vs_WFR2_.jpeg")),
       p_WFR2, width = 12, height = 5, dpi = 300)
cat("WFR1 vs WFR2 plot saved.\n")

# --- WFR1 vs nonWFR geo (pooled) ---
p_nonWFR <- plot_paired_dxy(all_comps,
                            "WFR1_MSL0", "nonWFR_geo_MSL0",
                            "WFR1_MSL4", "nonWFR_geo_MSL4",
                            "WFR1 vs nonWFR (all Brazil excl PI/MA/MS)")
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_WFR1_vs_nonWFR_geo_.jpeg")),
       p_nonWFR, width = 12, height = 5, dpi = 300)
cat("WFR1 vs nonWFR geo plot saved.\n")


# Combined panel: MG, MT, SC, PA (per-region)
combined_regions <- (p_MG | p_MT) / (p_SC | p_PA) +
  plot_annotation(
    title    = "Dxy along Chr31 — WFR1 vs qualifying regions",
    subtitle = "RED = MSL0-MSL0; BLUE = MSL4-MSL4. MSL locus position shown (dashed).",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )
ggsave(file.path(FIG_DIR, paste0("dxy_chr31_combined_.jpeg")),
       combined_regions, width = 16, height = 10, dpi = 300)
cat("Combined regions panel saved.\n")

cat("All Dxy plots saved to", FIG_DIR, "\n")







# 12. DXY RATIO (Dxy00 / Dxy44) ------------------------------------------------

# For each region pair, compute per-window ratio: Dxy_MSL0-MSL0 / Dxy_MSL4-MSL4.
# Ratio < 1: MSL0 strains more similar than MSL4 at that window.
# Ratio at MSL locus relative to background ratio = sweep signal.
# All 5 region pairs included — PA MSL4 and WFR2 MSL0 now have valid Dxy (n=1 ok).
# WFR2 ratio interpreted cautiously (WFR2 MSL0 n=1 — single sample).

# Reload TSV if not already in environment
all_comps <- read_tsv(file.path(OUT_DIR, paste0("windowed_dxy_fst_chr31_.tsv")),
                      show_col_types = FALSE)

# Region pairs for ratio calculation
ratio_regions <- list(
  list(label = "WFR1 vs WFR2", msl0_p1 = "WFR1_MSL0", msl0_p2 = "WFR2_MSL0",
       msl4_p1 = "WFR1_MSL4", msl4_p2 = "WFR2_MSL4", note = "WFR2 MSL0 n=1"),
  list(label = "WFR1 vs MG",   msl0_p1 = "WFR1_MSL0", msl0_p2 = "MG_MSL0",
       msl4_p1 = "WFR1_MSL4", msl4_p2 = "MG_MSL4",   note = ""),
  list(label = "WFR1 vs MT",   msl0_p1 = "WFR1_MSL0", msl0_p2 = "MT_MSL0",
       msl4_p1 = "WFR1_MSL4", msl4_p2 = "MT_MSL4",   note = ""),
  list(label = "WFR1 vs SC",   msl0_p1 = "WFR1_MSL0", msl0_p2 = "SC_MSL0",
       msl4_p1 = "WFR1_MSL4", msl4_p2 = "SC_MSL4",   note = ""),
  list(label = "WFR1 vs PA",   msl0_p1 = "WFR1_MSL0", msl0_p2 = "PA_MSL0",
       msl4_p1 = "WFR1_MSL4", msl4_p2 = "PA_MSL4",   note = "PA MSL4 n=1"),
  list(label = "WFR1 vs nonWFR geo", msl0_p1 = "WFR1_MSL0", msl0_p2 = "nonWFR_geo_MSL0",
       msl4_p1 = "WFR1_MSL4",        msl4_p2 = "nonWFR_geo_MSL4", note = "")
)

compute_ratio <- function(comps, msl0_p1, msl0_p2, msl4_p1, msl4_p2, label) {
  d0 <- comps |> filter(pop1 == msl0_p1, pop2 == msl0_p2) |>
    select(window_start, window_end, window_mid, n_snps, dxy) |>
    rename(dxy_msl0 = dxy)
  d4 <- comps |> filter(pop1 == msl4_p1, pop2 == msl4_p2) |>
    select(window_mid, dxy) |>
    rename(dxy_msl4 = dxy)
  d0 |>
    left_join(d4, by = "window_mid") |>
    mutate(
      region_pair = label,
      dxy_ratio   = dxy_msl0 / dxy_msl4   # NA if either is NA; Inf if dxy_msl4 == 0
    )
}

ratio_data <- bind_rows(lapply(ratio_regions, function(r) {
  compute_ratio(all_comps, r$msl0_p1, r$msl0_p2, r$msl4_p1, r$msl4_p2, r$label)
}))

write_tsv(ratio_data, file.path(OUT_DIR, paste0("dxy_ratio_chr31_.tsv")))
cat("Dxy ratio TSV saved.\n")


# Summary: mean ratio by region (background vs MSL flanking vs MSL locus)
MSL_START_R <- 1181282
MSL_END_R   <- 1192406
FLANK_PAD_R <- 50000

ratio_summary <- ratio_data |>
  mutate(region = case_when(
    window_mid >= MSL_START_R & window_mid <= MSL_END_R ~ "MSL_locus",
    window_mid >= (MSL_START_R - FLANK_PAD_R) &
      window_mid <= (MSL_END_R + FLANK_PAD_R)           ~ "MSL_flanking",
    TRUE                                                 ~ "background"
  )) |>
  group_by(region_pair, region) |>
  summarise(
    n_windows  = n(),
    mean_ratio = round(mean(dxy_ratio, na.rm = TRUE), 4),
    .groups    = "drop"
  )

print(ratio_summary, n = Inf)
write_tsv(ratio_summary, file.path(OUT_DIR, paste0("dxy_ratio_summary_chr31_.tsv")))
cat("Dxy ratio summary saved.\n")


# Ratio plots — one per region pair, horizontal reference line at ratio = 1
plot_ratio <- function(data, region_label) {
  d <- data |> filter(region_pair == region_label)
  y_max <- max(c(d$dxy_ratio, 2), na.rm = TRUE) * 1.05
  y_max <- min(y_max, 5)   # cap extreme values for readability
  
  ggplot(d, aes(x = window_mid / 1e6, y = dxy_ratio)) +
    geom_hline(yintercept = 1, colour = "grey60", linetype = "dashed",
               linewidth = 0.5) +
    geom_vline(xintercept = MSL_MID, colour = "grey40",
               linetype = "dashed", linewidth = 0.6) +
    annotate("text", x = MSL_MID + 0.02, y = y_max * 0.95,
             label = "MSL", colour = "grey40", size = 3,
             fontface = "bold", hjust = 0) +
    geom_line(colour = "#BB5566", linewidth = 0.6) +
    geom_point(colour = "#BB5566", size = 1.5, alpha = 0.7) +
    geom_smooth(method = "loess", span = 0.3, se = FALSE,
                colour = "#4477AA", linewidth = 0.8, linetype = "solid") +
    scale_x_continuous(name = "Chr31 position (Mb)") +
    scale_y_continuous(name = "Dxy ratio (MSL0-MSL0 / MSL4-MSL4)",
                       limits = c(0, y_max)) +
    labs(title = region_label,
         subtitle = "Ratio < 1: MSL0 more similar than MSL4. Dashed line = ratio of 1 (no difference).") +
    theme_classic(base_size = 11) +
    theme(
      plot.title    = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, colour = "grey40")
    )
}

p_ratio_WFR2 <- plot_ratio(ratio_data, "WFR1 vs WFR2")
p_ratio_MG   <- plot_ratio(ratio_data, "WFR1 vs MG")
p_ratio_MT   <- plot_ratio(ratio_data, "WFR1 vs MT")
p_ratio_SC   <- plot_ratio(ratio_data, "WFR1 vs SC")
p_ratio_PA   <- plot_ratio(ratio_data, "WFR1 vs PA")

# Individual saves
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_WFR1_vs_WFR2_.jpeg")),
       p_ratio_WFR2, width = 12, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_WFR1_vs_MG_.jpeg")),
       p_ratio_MG, width = 12, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_WFR1_vs_MT_.jpeg")),
       p_ratio_MT, width = 12, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_WFR1_vs_SC_.jpeg")),
       p_ratio_SC, width = 12, height = 5, dpi = 300)
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_WFR1_vs_PA_.jpeg")),
       p_ratio_PA, width = 12, height = 5, dpi = 300)

p_ratio_nonWFR <- plot_ratio(ratio_data, "WFR1 vs nonWFR geo")
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_WFR1_vs_nonWFR_geo_.jpeg")),
       p_ratio_nonWFR, width = 12, height = 5, dpi = 300)
cat("Individual ratio plots saved.\n")

# Combined panel — 2x2: MG, MT, SC, nonWFR geo (key comparisons)
ratio_combined <- (p_ratio_MG | p_ratio_MT) / (p_ratio_SC | p_ratio_nonWFR) +
  plot_annotation(
    title    = "Dxy ratio (MSL0-MSL0 / MSL4-MSL4) along Chr31 — WFR1 vs qualifying regions",
    subtitle = "Ratio < 1 = MSL0 strains more similar than MSL4 at that window. Loess smooth in blue.",
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  )
ggsave(file.path(FIG_DIR, paste0("dxy_ratio_combined_.jpeg")),
       ratio_combined, width = 16, height = 10, dpi = 300)
cat("Combined ratio panel saved.\n")

cat("All ratio outputs saved to", OUT_DIR, "and", FIG_DIR, "\n")





# 13. PERCENTILE WINDOW STATISTICAL TEST (5% and 10%) --------------------------

# Per MSL0-MSL0 comparison: do MSL flanking windows fall within the lowest
# 5% or 10% of Dxy values along chr31?
# Run on MSL0-MSL0 only — sweep signal expected there, not MSL4-MSL4.

all_comps <- read_tsv(file.path(OUT_DIR, "windowed_dxy_fst_chr31_.tsv"),
                      show_col_types = FALSE)

MSL_START_P <- 1181282
MSL_END_P   <- 1192406
FLANK_PAD_P <- 50000

msl0_comps <- list(
  list(p1 = "WFR1_MSL0", p2 = "WFR2_MSL0", label = "WFR1 vs WFR2"),
  list(p1 = "WFR1_MSL0", p2 = "MG_MSL0",   label = "WFR1 vs MG"),
  list(p1 = "WFR1_MSL0", p2 = "MT_MSL0",   label = "WFR1 vs MT"),
  list(p1 = "WFR1_MSL0", p2 = "SC_MSL0",   label = "WFR1 vs SC"),
  list(p1 = "WFR1_MSL0", p2 = "PA_MSL0",           label = "WFR1 vs PA"),
  list(p1 = "WFR1_MSL0", p2 = "nonWFR_geo_MSL0",   label = "WFR1 vs nonWFR geo")
)

run_pct_test <- function(comp, threshold) {
  d <- all_comps |>
    filter(pop1 == comp$p1, pop2 == comp$p2) |>
    mutate(zone = case_when(
      window_mid >= MSL_START_P & window_mid <= MSL_END_P ~ "MSL_locus",
      window_mid >= (MSL_START_P - FLANK_PAD_P) &
        window_mid <= (MSL_END_P + FLANK_PAD_P) ~ "MSL_flanking",
      TRUE ~ "background"
    ))
  
  pct_val <- quantile(d$dxy, threshold, na.rm = TRUE)
  d       <- d |> mutate(below_pct = !is.na(dxy) & dxy <= pct_val)
  
  msl_win <- d |> filter(zone == "MSL_flanking")
  n_msl   <- nrow(msl_win)
  n_below <- sum(msl_win$below_pct, na.rm = TRUE)
  
  data.frame(
    comparison    = comp$label,
    threshold_pct = paste0(threshold * 100, "%"),
    pct_dxy_value = round(pct_val, 6),
    n_msl_windows = n_msl,
    n_below       = n_below,
    pct_below     = round(100 * n_below / n_msl, 1),
    result        = paste0(n_below, "/", n_msl,
                           " MSL flanking windows in lowest ",
                           threshold * 100, "%"),
    stringsAsFactors = FALSE
  )
}

pct_results <- do.call(rbind, lapply(msl0_comps, function(comp) {
  rbind(
    run_pct_test(comp, 0.05),
    run_pct_test(comp, 0.10)
  )
}))

print(pct_results)
write.table(pct_results, file.path(OUT_DIR, "dxy_pct_test_chr31_.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
cat("Percentile window test results saved.\n")



