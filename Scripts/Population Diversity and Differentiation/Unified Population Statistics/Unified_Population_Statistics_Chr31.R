

# UNIFIED POPULATION STATISTICS — CHR31 (PSEUDODIPLOID / TETRAPLOID)
# Overhauled to fix mathematical biases, window-averaging errors, and ploidy scaling.



# 1. LIBRARIES -----------------------------------------------------------------
suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(vcfR)
})



# 2. COMMAND-LINE ARGUMENTS ----------------------------------------------------
parser <- ArgumentParser(description = "Windowed population stats along chr31.")
parser$add_argument("--vcf_path", type = "character", default = "Data/tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.vcf.gz")
parser$add_argument("--meta_path", type = "character", default = "Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
parser$add_argument("--mem_path", type = "character", default = "Data/population_structure/population_membership_table_Americas_only_k13_AM_notation.tsv")
parser$add_argument("--out_dir", type = "character", default = "Data/Unified_Popgen")
args <- parser$parse_args()



# 3. CONSTANTS -----------------------------------------------------------------
MSL_START  <- 1181282L
MSL_END    <- 1192406L
WIN_SIZE   <- 50000L
WIN_STEP   <-  5000L

# CRITICAL PLOIDY SETTING: Change to 4L ONLY if VCF genotypes look like "1/1/0/1".
# If your data was parsed/called as pseudodiploid ("0/1", "1/1"), this MUST remain 2L.
PLOIDY     <- 2L         

MIN_N_SNPS   <- 3L         # Minimum SNPs to evaluate a window
MIN_FST_SNPS <- 10L        # Minimum valid SNPs required for robust window Fst
BATCH_SIZE   <- 50L        
MIN_GROUP_N  <- 5L         # Enforced minimum samples per group



# 4. OUTPUT DIRECTORY AND LOGGING ----------------------------------------------
dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)
timestamp <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
log_path  <- file.path(args$out_dir, paste0("unified_popgen_chr31_", timestamp, ".log"))
log_con   <- file(log_path, open = "wt")
sink(log_con, split = TRUE)

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"), sep = "")
}

log_msg("====================================================")
log_msg(paste0("Ploidy Model Enforced: ", PLOIDY))
log_msg("====================================================")



# 5. LOAD VCF AND BUILD DOSAGE MATRIX ------------------------------------------
log_msg("Loading VCF...")
vcf <- read.vcfR(args$vcf_path, verbose = FALSE)
chrom <- getCHROM(vcf)
pos   <- as.numeric(getPOS(vcf))

keep_snps  <- !(chrom == "LinJ.31" & pos >= MSL_START & pos <= MSL_END)
vcf_no_msl <- vcf[keep_snps, ]
pos_no_msl <- pos[keep_snps]

count_alt_alleles <- function(x) lengths(regmatches(x, gregexpr("1", x)))
gt_strings    <- extract.gt(vcf_no_msl, element = "GT")
dosage_matrix <- matrix(count_alt_alleles(gt_strings), nrow = nrow(gt_strings), ncol = ncol(gt_strings))
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)

dosage_matrix[is.na(gt_strings) | gt_strings %in% c("./.", "./././.")] <- NA
dosage_matrix <- t(dosage_matrix)



# 6. LOAD METADATA AND DEFINE ELIGIBLE GROUPS ----------------------------------
meta <- read_tsv(args$meta_path, show_col_types = FALSE) |> select(sample_id, MSL_coverage_rounded)
mem  <- read_tsv(args$mem_path, show_col_types = FALSE) |> select(sample, population)

pop_msl <- mem |>
  left_join(meta, by = c("sample" = "sample_id")) |>
  filter(MSL_coverage_rounded %in% c(0L, 4L)) |>
  mutate(group = paste0(population, "_", MSL_coverage_rounded))

group_counts <- pop_msl |> count(group) |> filter(n >= MIN_GROUP_N) |> arrange(group)
eligible_groups <- group_counts$group

group_samples <- setNames(lapply(eligible_groups, function(g) {
  intersect(pop_msl$sample[pop_msl$group == g], rownames(dosage_matrix))
}), eligible_groups)



# 7. GENERATE POPULATION PAIRS -------------------------------------------------
msl0_groups <- eligible_groups[grepl("_0$", eligible_groups)]
msl4_groups <- eligible_groups[grepl("_4$", eligible_groups)]
all_pairs   <- c(combn(msl0_groups, 2, simplify = FALSE), combn(msl4_groups, 2, simplify = FALSE))



# 8. STAT FUNCTIONS (CORRECTED MATHEMATICS) ------------------------------------

# 8a. Nucleotide diversity (pi)
compute_pi <- function(dos_win, ploidy = PLOIDY) {
  if (is.null(dos_win) || ncol(dos_win) == 0) return(NA_real_)
  
  # Total alleles sampled per individual = sample count * ploidy
  per_snp <- apply(dos_win, 2, function(x) {
    x <- x[!is.na(x)]
    n_alleles <- length(x) * ploidy
    if (n_alleles < 2) return(0)
    p <- sum(x) / n_alleles
    if (p <= 0 || p >= 1) return(0)
    (n_alleles / (n_alleles - 1)) * 2 * p * (1 - p)
  })
  sum(per_snp) / WIN_SIZE
}

# 8b. Tajima's D
compute_tajima <- function(dos_win, ploidy = PLOIDY) {
  if (is.null(dos_win) || ncol(dos_win) == 0) return(NA_real_)
  
  # Drop columns with excessive missing data to preserve n stability across the window
  valid_counts <- colSums(!is.na(dos_win))
  max_samples  <- nrow(dos_win)
  dos_win      <- dos_win[, valid_counts >= (max_samples * 0.8), drop = FALSE]
  
  if (ncol(dos_win) == 0) return(NA_real_)
  
  # Sample size tracking in units of independent chromosomes/alleles
  n_alleles <- max_samples * ploidy
  if (n_alleles < 4) return(NA_real_)
  
  per_snp <- apply(dos_win, 2, function(x) {
    x <- x[!is.na(x)]
    n_local <- length(x) * ploidy
    if (n_local < 2) return(c(seg = 0L, k_snp = 0))
    p <- sum(x) / n_local
    if (p <= 0 || p >= 1) return(c(seg = 0L, k_snp = 0))
    k_snp <- (n_local / (n_local - 1)) * 2 * p * (1 - p)
    c(seg = 1L, k_snp = k_snp)
  })
  
  S <- sum(per_snp["seg", ])
  if (S < 2) return(NA_real_)
  k <- sum(per_snp["k_snp", ])
  
  # Math Correction: Variance factors scaled strictly to independent chromosome count
  a1  <- sum(1 / seq_len(n_alleles - 1))
  a2  <- sum(1 / seq_len(n_alleles - 1)^2)
  b1  <- (n_alleles + 1) / (3 * (n_alleles - 1))
  b2  <- 2 * (n_alleles^2 + n_alleles + 3) / (9 * n_alleles * (n_alleles - 1))
  c1  <- b1 - 1 / a1
  c2  <- b2 - (n_alleles + 2) / (a1 * n_alleles) + a2 / a1^2
  e1  <- c1 / a1
  e2  <- c2 / (a1^2 + a2)
  
  var_d <- e1 * S + e2 * S * (S - 1)
  if (var_d <= 0) return(NA_real_)
  
  (k - S / a1) / sqrt(var_d)
}

# 8c. Dxy (Ratio of Averages implementation)
compute_window_dxy <- function(dos_win, ids1, ids2, ploidy = PLOIDY) {
  p1 <- intersect(ids1, rownames(dos_win))
  p2 <- intersect(ids2, rownames(dos_win))
  if (length(p1) < 1 || length(p2) < 1) return(NA_real_)
  
  mat1 <- dos_win[p1, , drop = FALSE]
  mat2 <- dos_win[p2, , drop = FALSE]
  
  snp_dxy <- sapply(seq_len(ncol(dos_win)), function(j) {
    a <- mat1[, j]; b <- mat2[, j]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    mean(outer(a, b, function(x, y) abs(x - y) / ploidy))
  })
  
  # Average over valid SNPs present inside the window
  mean(snp_dxy, na.rm = TRUE)
}

# 8d. Fst (Hudson Ratio of Averages implementation with unequal sample weights)
compute_window_fst <- function(dos_win, ids1, ids2, ploidy = PLOIDY) {
  p1 <- intersect(ids1, rownames(dos_win))
  p2 <- intersect(ids2, rownames(dos_win))
  if (length(p1) < MIN_GROUP_N || length(p2) < MIN_GROUP_N) return(NA_real_)
  
  mat1 <- dos_win[p1, , drop = FALSE]
  mat2 <- dos_win[p2, , drop = FALSE]
  
  num_sum <- 0
  den_sum <- 0
  count_valid <- 0
  
  for (j in seq_len(ncol(dos_win))) {
    a <- mat1[, j]; b <- mat2[, j]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    
    if (length(a) < MIN_GROUP_N || length(b) < MIN_GROUP_N) next
    
    n1 <- length(a) * ploidy
    n2 <- length(b) * ploidy
    f1 <- sum(a) / n1
    f2 <- sum(b) / n2
    
    # Corrected Unbiased Hudson Components
    num <- (f1 - f2)^2 - (f1 * (1 - f1) / (n1 - 1)) - (f2 * (1 - f2) / (n2 - 1))
    den <- f1 * (1 - f2) + f2 * (1 - f1)
    
    if (den > 0) {
      num_sum <- num_sum + num
      den_sum <- den_sum + den
      count_valid <- count_valid + 1
    }
  }
  
  if (count_valid < MIN_FST_SNPS || den_sum == 0) return(NA_real_)
  
  fst <- num_sum / den_sum
  return(max(0, min(1, fst))) # Clamp to realistic 0-1 bounds
}



# 9. WINDOW COORDINATES --------------------------------------------------------
chr_min    <- min(pos_no_msl)
chr_max    <- max(pos_no_msl)
win_starts <- seq(chr_min, chr_max - WIN_SIZE + 1L, by = WIN_STEP)



# 10. INITIALISE OUTPUT FILES --------------------------------------------------
within_path  <- file.path(args$out_dir, "within_pop_chr31.tsv")
between_path <- file.path(args$out_dir, "between_pop_chr31.tsv")

write_tsv(tibble(group = character(), window_start = integer(), window_end = integer(), window_mid = numeric(), n_snps = integer(), pi = numeric(), tajima_d = numeric()), within_path)
write_tsv(tibble(pop1 = character(), pop2 = character(), msl_type = character(), window_start = integer(), window_end = integer(), window_mid = numeric(), n_snps = integer(), dxy = numeric(), fst = numeric()), between_path)



# 11. MAIN LOOP ----------------------------------------------------------------
within_buffer  <- list()
between_buffer <- list()
n_done         <- 0L

for (i in seq_along(win_starts)) {
  w_start <- win_starts[i]
  w_end   <- w_start + WIN_SIZE - 1L
  w_mid   <- (w_start + w_end) / 2
  
  in_win <- which(pos_no_msl >= w_start & pos_no_msl <= w_end)
  n_snps <- length(in_win)
  dos_win <- if (n_snps >= MIN_N_SNPS) dosage_matrix[, in_win, drop = FALSE] else NULL
  
  # (a) Within-population
  for (grp in eligible_groups) {
    ids     <- group_samples[[grp]]
    dos_grp <- if (!is.null(dos_win)) dos_win[intersect(ids, rownames(dos_win)), , drop = FALSE] else NULL
    
    pi_val <- if (!is.null(dos_grp) && nrow(dos_grp) >= MIN_GROUP_N) compute_pi(dos_grp) else NA_real_
    td_val <- if (!is.null(dos_grp) && nrow(dos_grp) >= MIN_GROUP_N) compute_tajima(dos_grp) else NA_real_
    
    within_buffer[[length(within_buffer) + 1]] <- tibble(group = grp, window_start = w_start, window_end = w_end, window_mid = w_mid, n_snps = n_snps, pi = pi_val, tajima_d = td_val)
  }
  
  # (b) Between-population
  for (pair in all_pairs) {
    grp1     <- pair[[1]]
    grp2     <- pair[[2]]
    ids1     <- group_samples[[grp1]]
    ids2     <- group_samples[[grp2]]
    msl_type <- if (grepl("_0$", grp1)) "MSL0-MSL0" else "MSL4-MSL4"
    
    dxy_val <- if (!is.null(dos_win)) compute_window_dxy(dos_win, ids1, ids2) else NA_real_
    fst_val <- if (!is.null(dos_win)) compute_window_fst(dos_win, ids1, ids2) else NA_real_
    
    between_buffer[[length(between_buffer) + 1]] <- tibble(pop1 = grp1, pop2 = grp2, msl_type = msl_type, window_start = w_start, window_end = w_end, window_mid = w_mid, n_snps = n_snps, dxy = dxy_val, fst = fst_val)
  }
  
  n_done <- n_done + 1L
  if (n_done %% BATCH_SIZE == 0L || i == length(win_starts)) {
    write_tsv(bind_rows(within_buffer),  within_path,  append = TRUE)
    write_tsv(bind_rows(between_buffer), between_path, append = TRUE)
    within_buffer  <- list()
    between_buffer <- list()
    log_msg(sprintf("Windows processed: %d / %d", n_done, length(win_starts)))
  }
}

log_msg("Run Complete.")
sink(); close(log_con)



