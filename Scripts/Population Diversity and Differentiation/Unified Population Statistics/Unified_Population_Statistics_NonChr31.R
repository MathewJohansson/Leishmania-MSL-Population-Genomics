# UNIFIED PER-POPULATION STATISTICS — NON-CHR31 (DIPLOID, ALL CHROMOSOMES)

# Companion script to Unified_Population_Statistics_Chr31.R. Computes the same four
# statistics in 50kb sliding windows (5kb step) across all non-chr31
# chromosomes, using the diploid whole-genome VCF.

# The primary purpose of this script is to provide a genome-wide baseline
# for comparison against the chr31 signal — specifically, to test whether
# the MSL0-MSL0 Dxy depression seen on chr31 is chr31-specific or reflects
# a genome-wide demographic effect.

# KEY DIFFERENCE FROM CHR31 SCRIPT: PLOIDY = 2
#   Dosage values are 0, 1, or 2 (not 0–4).
#   Allele frequency = dosage / 2.
#   All stat functions are otherwise identical — they use the ploidy argument.

# OTHER DIFFERENCES FROM CHR31 SCRIPT:
#   - VCF is the diploid non-chr31 file (chr31 already absent; no MSL exclusion)
#   - Window loop iterates per chromosome (positions restart at 1 per chromosome)
#   - Output TSVs include a `chrom` column
#   - Output filenames use "_nonchr31" suffix

# STATISTICS, GROUPS, AND PAIRS: identical to chr31 script.
#   Within-population (per group per window): pi, tajima_d
#   Between-population (per same-MSL-type pair per window): dxy, fst
#   Eligible groups (n=13): AM1_0 AM2_0 AM2_4 AM3_0 AM4_4 AM5_0 AM6_0
#                            AM7_0 AM8_4 AM9_0 AM10_0 AM11_0 AM12_0
#   Pairs: 45 MSL0-MSL0 + 3 MSL4-MSL4 = 48 total

# INPUTS (command-line args; defaults = Viking paths since this runs on Viking)
#   --vcf_path  : diploid non-chr31 VCF (.vcf.gz)
#   --meta_path : metadata TSV with MSL_coverage_rounded per sample
#   --mem_path  : AM membership TSV (sample, population, proportion)
#   --out_dir   : directory to write TSV outputs and log

# OUTPUTS (in --out_dir)
#   within_pop_nonchr31.tsv   — pi and tajima_d per group per chromosome per window
#   between_pop_nonchr31.tsv  — dxy and fst per pair per chromosome per window
#   unified_popgen_nonchr31_TIMESTAMP.log

# Run on Viking:
#   Rscript Unified_Population_Statistics_Other_Chrs.R \
#     --vcf_path  /mnt/scratch/.../vcf_files/Linfantum_samples_no_outgroup_biallelic.2025-06-10.noLinJ.31.SNPsonly.vcf.gz \
#     --meta_path /mnt/scratch/.../scripting/vcf_samples29_...tsv \
#     --mem_path  /mnt/scratch/.../scripting/population_membership_table_Americas_only_k13_AM_notation.tsv \
#     --out_dir   /mnt/scratch/.../popgen_unified/

# NOTE: The non-chr31 VCF is large (whole genome). Loading it with read.vcfR
# may take several minutes and ~8–16GB RAM. Run via SLURM on Viking.




# 1. LIBRARIES -----------------------------------------------------------------

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(vcfR)
})




# 2. COMMAND-LINE ARGUMENTS ----------------------------------------------------

parser <- ArgumentParser(description = "Windowed population stats across non-chr31 chromosomes (diploid).")

parser$add_argument(
  "--vcf_path", type = "character",
  default = "/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/vcf_files/Linfantum_samples_no_outgroup_biallelic.2025-06-10.noLinJ.31.SNPsonly.vcf.gz",
  help    = "Path to diploid non-chr31 VCF (.vcf.gz)"
)
parser$add_argument(
  "--meta_path", type = "character",
  default = "/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/scripting/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv",
  help    = "Path to metadata TSV (must contain sample_id and MSL_coverage_rounded)"
)
parser$add_argument(
  "--mem_path", type = "character",
  default = "/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/scripting/population_membership_table_Americas_only_k13_AM_notation.tsv",
  help    = "Path to AM membership TSV (must contain sample and population columns)"
)
parser$add_argument(
  "--out_dir", type = "character",
  default = "/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/popgen_unified",
  help    = "Output directory for TSVs and log file"
)

args <- parser$parse_args()




# 3. CONSTANTS -----------------------------------------------------------------

# PLOIDY = 2: this is the only fundamental difference from the chr31 script.
# All stat functions divide dosage by ploidy, so they work correctly for both.
PLOIDY     <- 2L         # diploid: allele frequency = dosage / 2

WIN_SIZE   <- 50000L     # 50kb window size  — matches chr31 script
WIN_STEP   <-  5000L     # 5kb step          — matches chr31 script

MIN_N_SNPS <- 3L         # minimum SNPs in window to compute any statistic
BATCH_SIZE <- 50L        # write to TSV every this many windows (across all chromosomes)
MIN_GROUP_N <- 5L        # minimum samples per group to be eligible

# No MSL_START / MSL_END constants here — chr31 is already absent from this VCF.




# 4. OUTPUT DIRECTORY AND LOGGING ----------------------------------------------

dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y-%m-%dT%H%M%S")
log_path  <- file.path(args$out_dir, paste0("Unified_Population_Statistics_Other_Chrs_", timestamp, ".log"))
log_con   <- file(log_path, open = "wt")
sink(log_con, split = TRUE)

log_msg <- function(...) {
  cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"), sep = "")
}

log_msg("====================================================")
log_msg("Unified_Population_Statistics_Other_Chrs.R")
log_msg(paste0("Started:    ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
log_msg(paste0("VCF:        ", args$vcf_path))
log_msg(paste0("Metadata:   ", args$meta_path))
log_msg(paste0("Membership: ", args$mem_path))
log_msg(paste0("Out dir:    ", args$out_dir))
log_msg(paste0("Window:     ", WIN_SIZE / 1000, "kb size, ", WIN_STEP / 1000, "kb step"))
log_msg(paste0("Ploidy:     ", PLOIDY, " (diploid — dosage 0-2, freq = dosage / 2)"))
log_msg("====================================================")




# 5. LOAD VCF AND BUILD DOSAGE MATRIX ------------------------------------------

# Identical logic to the chr31 script, but:
#   - No MSL exclusion step (chr31 is already absent from this VCF)
#   - Missing GT string is "./." (diploid) not "./././." (tetraploid)
#   - Dosage values will be 0, 1, or 2 (not 0–4)

# We retain both chrom_vec and pos_vec after loading — needed because the
# window loop iterates per chromosome (positions restart at 1 per chromosome).

log_msg("Loading VCF... (this may take several minutes for the whole-genome file)")
vcf <- read.vcfR(args$vcf_path, verbose = FALSE)

chrom_vec <- getCHROM(vcf)
pos_vec   <- as.numeric(getPOS(vcf))

log_msg(sprintf("VCF loaded: %d SNPs across %d chromosomes",
                length(pos_vec), length(unique(chrom_vec))))

# Count alt alleles per genotype string (same function as chr31 script)
# "0/0" → 0, "0/1" → 1, "1/1" → 2
count_alt_alleles <- function(x) lengths(regmatches(x, gregexpr("1", x)))

gt_strings    <- extract.gt(vcf, element = "GT")
dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings),
  ncol = ncol(gt_strings)
)
rownames(dosage_matrix) <- rownames(gt_strings)   # SNP identifiers
colnames(dosage_matrix) <- colnames(gt_strings)   # sample IDs

# Set missing diploid genotypes to NA ("./." only — no tetraploid pattern)
dosage_matrix[is.na(gt_strings) | gt_strings == "./."] <- NA

# Transpose: rows = samples, columns = SNPs
dosage_matrix <- t(dosage_matrix)

log_msg(sprintf("Dosage matrix built: %d samples x %d SNPs",
                nrow(dosage_matrix), ncol(dosage_matrix)))
log_msg(sprintf("Dosage range check: min=%d, max=%d (expect 0–2 for diploid)",
                min(dosage_matrix, na.rm = TRUE),
                max(dosage_matrix, na.rm = TRUE)))




# 6. LOAD METADATA AND DEFINE ELIGIBLE GROUPS ----------------------------------

# Identical to chr31 script — same metadata, same membership TSV, same groups.
# Groups are defined from population × MSL-type (e.g. AM2_0, AM8_4).

log_msg("Loading metadata and membership TSVs...")

meta <- read_tsv(args$meta_path, show_col_types = FALSE) |>
  select(sample_id, MSL_coverage_rounded)

mem <- read_tsv(args$mem_path, show_col_types = FALSE) |>
  select(sample, population)

pop_msl <- mem |>
  left_join(meta, by = c("sample" = "sample_id")) |>
  filter(MSL_coverage_rounded %in% c(0L, 4L)) |>
  mutate(group = paste0(population, "_", MSL_coverage_rounded))

group_counts <- pop_msl |>
  count(group) |>
  filter(n >= MIN_GROUP_N) |>
  arrange(group)

eligible_groups <- group_counts$group

log_msg(sprintf("Eligible groups (n >= %d): %d", MIN_GROUP_N, length(eligible_groups)))
log_msg(paste(" ", eligible_groups, collapse = "\n"))

group_samples <- setNames(
  lapply(eligible_groups, function(g) {
    ids <- pop_msl$sample[pop_msl$group == g]
    intersect(ids, rownames(dosage_matrix))
  }),
  eligible_groups
)

for (g in eligible_groups) {
  n_meta <- sum(pop_msl$group == g)
  n_vcf  <- length(group_samples[[g]])
  if (n_vcf < n_meta) {
    log_msg(sprintf("  WARNING: %s — %d in metadata but only %d in VCF", g, n_meta, n_vcf))
  }
  log_msg(sprintf("  %s: n = %d", g, n_vcf))
}




# 7. GENERATE POPULATION PAIRS -------------------------------------------------

# Identical to chr31 script.

msl0_groups <- eligible_groups[grepl("_0$", eligible_groups)]
msl4_groups <- eligible_groups[grepl("_4$", eligible_groups)]

msl0_pairs  <- combn(msl0_groups, 2, simplify = FALSE)
msl4_pairs  <- combn(msl4_groups, 2, simplify = FALSE)
all_pairs   <- c(msl0_pairs, msl4_pairs)

log_msg(sprintf("MSL0 groups: %d  →  %d MSL0-MSL0 pairs",
                length(msl0_groups), length(msl0_pairs)))
log_msg(sprintf("MSL4 groups: %d  →  %d MSL4-MSL4 pairs",
                length(msl4_groups), length(msl4_pairs)))
log_msg(sprintf("Total pairs: %d", length(all_pairs)))




# 8. STAT FUNCTIONS ------------------------------------------------------------

# Identical to chr31 script. PLOIDY = 2 is passed as default argument,
# so all allele frequencies are computed as dosage / 2 throughout.

# 8a. Nucleotide diversity (pi)
# π is per-base: sum of per-site heterozygosity across all SNPs in the window,
# divided by window length (WIN_SIZE). Monomorphic and missing SNPs contribute
# 0 — they are real bases in the window and must be counted in the denominator.
compute_pi <- function(dos_win, ploidy = PLOIDY) {
  if (is.null(dos_win) || ncol(dos_win) == 0) return(NA_real_)
  per_snp <- apply(dos_win, 2, function(x) {
    x <- x[!is.na(x)]
    n <- length(x)
    if (n < 2) return(0)
    p <- mean(x) / ploidy
    if (p <= 0 || p >= 1) return(0)
    (n / (n - 1)) * 2 * p * (1 - p)
  })
  sum(per_snp) / WIN_SIZE
}

# 8b. Tajima's D
compute_tajima <- function(dos_win, ploidy = PLOIDY) {
  if (is.null(dos_win) || ncol(dos_win) == 0) return(NA_real_)
  n <- nrow(dos_win)
  if (n < 4) return(NA_real_)
  per_snp <- apply(dos_win, 2, function(x) {
    x   <- x[!is.na(x)]
    n_x <- length(x)
    if (n_x < 2) return(c(seg = 0L, k_snp = 0))
    p <- mean(x) / ploidy
    if (p <= 0 || p >= 1) return(c(seg = 0L, k_snp = 0))
    c(seg = 1L, k_snp = (n_x / (n_x - 1)) * 2 * p * (1 - p))
  })
  S <- sum(per_snp["seg", ])
  if (S < 2) return(NA_real_)
  k <- sum(per_snp["k_snp", ])
  a1  <- sum(1 / seq_len(n - 1))
  a2  <- sum(1 / seq_len(n - 1)^2)
  b1  <- (n + 1) / (3 * (n - 1))
  b2  <- 2 * (n^2 + n + 3) / (9 * n * (n - 1))
  c1  <- b1 - 1 / a1
  c2  <- b2 - (n + 2) / (a1 * n) + a2 / a1^2
  e1  <- c1 / a1
  e2  <- c2 / (a1^2 + a2)
  var_d <- e1 * S + e2 * S * (S - 1)
  if (var_d <= 0) return(NA_real_)
  (k - S / a1) / sqrt(var_d)
}

# 8c. Dxy
compute_dxy <- function(dos_win, ids1, ids2, ploidy = PLOIDY) {
  p1 <- intersect(ids1, rownames(dos_win))
  p2 <- intersect(ids2, rownames(dos_win))
  if (length(p1) < 1 || length(p2) < 1) return(NA_real_)
  mat1 <- dos_win[p1, , drop = FALSE]
  mat2 <- dos_win[p2, , drop = FALSE]
  per_snp <- sapply(seq_len(ncol(dos_win)), function(j) {
    a <- mat1[, j]; b <- mat2[, j]
    a <- a[!is.na(a)]; b <- b[!is.na(b)]
    if (length(a) == 0 || length(b) == 0) return(NA_real_)
    mean(outer(a, b, function(x, y) abs(x - y) / ploidy), na.rm = TRUE)
  })
  mean(per_snp, na.rm = TRUE)
}

# 8d. Fst
compute_fst <- function(dos_win, ids1, ids2, ploidy = PLOIDY) {
  p1 <- intersect(ids1, rownames(dos_win))
  p2 <- intersect(ids2, rownames(dos_win))
  if (length(p1) < 2 || length(p2) < 2) return(NA_real_)
  mat1 <- dos_win[p1, , drop = FALSE]
  mat2 <- dos_win[p2, , drop = FALSE]
  per_snp <- sapply(seq_len(ncol(dos_win)), function(j) {
    a <- mat1[, j] / ploidy; b <- mat2[, j] / ploidy
    a <- a[!is.na(a)];       b <- b[!is.na(b)]
    if (length(a) < 2 || length(b) < 2) return(NA_real_)
    f1 <- mean(a); f2 <- mean(b)
    n1 <- length(a); n2 <- length(b)
    pi1  <- 2 * f1 * (1 - f1) * n1 / (n1 - 1)
    pi2  <- 2 * f2 * (1 - f2) * n2 / (n2 - 1)
    pi_w <- (pi1 + pi2) / 2
    pi_b <- f1 * (1 - f2) + f2 * (1 - f1)
    if (pi_b == 0) return(NA_real_)
    (pi_b - pi_w) / pi_b
  })
  fst <- mean(per_snp, na.rm = TRUE)
  max(0, min(1, fst))
}




# 9. CHROMOSOMES TO PROCESS ----------------------------------------------------

# Extract the list of chromosomes present in the VCF and sort them.
# Chr31 (LinJ.31) should not be present — the VCF filename confirms this
# (noLinJ.31), but we check and warn if it appears unexpectedly.

chroms_available <- sort(unique(chrom_vec))

if ("LinJ.31" %in% chroms_available) {
  log_msg("WARNING: LinJ.31 found in VCF — expected noLinJ.31 file. Check VCF path.")
}

log_msg(sprintf("Chromosomes to process: %d", length(chroms_available)))
log_msg(paste(" ", chroms_available, collapse = "\n"))

# Estimate total windows across all chromosomes
total_windows_est <- sum(sapply(chroms_available, function(chr) {
  p       <- pos_vec[chrom_vec == chr]
  chr_end <- max(p) - WIN_SIZE + 1L
  if (chr_end < min(p)) return(0L)
  length(seq(min(p), chr_end, by = WIN_STEP))
}))
log_msg(sprintf("Estimated total windows across all chromosomes: ~%d", total_windows_est))




# 10. INITIALISE OUTPUT FILES --------------------------------------------------

# Same structure as chr31 script but with an additional `chrom` column.
# Headers written first; results appended in batches.

within_path  <- file.path(args$out_dir, "within_pop_nonchr31.tsv")
between_path <- file.path(args$out_dir, "between_pop_nonchr31.tsv")

write_tsv(
  tibble(
    chrom        = character(),   # chromosome — additional vs chr31 script
    group        = character(),
    window_start = integer(),
    window_end   = integer(),
    window_mid   = numeric(),
    n_snps       = integer(),
    pi           = numeric(),
    tajima_d     = numeric()
  ),
  within_path
)

write_tsv(
  tibble(
    chrom        = character(),   # chromosome — additional vs chr31 script
    pop1         = character(),
    pop2         = character(),
    msl_type     = character(),
    window_start = integer(),
    window_end   = integer(),
    window_mid   = numeric(),
    n_snps       = integer(),
    dxy          = numeric(),
    fst          = numeric()
  ),
  between_path
)

log_msg("Output files initialised:")
log_msg(sprintf("  Within:  %s", within_path))
log_msg(sprintf("  Between: %s", between_path))




# 11. MAIN LOOP ----------------------------------------------------------------

# Outer loop: chromosomes.
# Middle loop: windows within each chromosome.
# Inner loops: (a) eligible groups → pi and tajima_d
#              (b) same-MSL-type pairs → dxy and fst

# Chromosomes are processed sequentially. Windows are defined per chromosome
# so that no window spans a chromosome boundary.

# The batch write counter (n_done) runs across all chromosomes — results are
# flushed every BATCH_SIZE windows regardless of chromosome boundaries.

log_msg("Starting main loop...")
run_start <- proc.time()

within_buffer  <- list()
between_buffer <- list()
n_done         <- 0L          # total windows processed across all chromosomes
total_wins     <- 0L          # running total of windows (for final summary)

for (chr in chroms_available) {
  
  # Column indices in the dosage matrix for SNPs on this chromosome
  chr_cols <- which(chrom_vec == chr)
  pos_chr  <- pos_vec[chr_cols]
  
  if (length(pos_chr) < 2) {
    log_msg(sprintf("  Skipping %s — fewer than 2 SNPs", chr))
    next
  }
  
  chr_end <- max(pos_chr) - WIN_SIZE + 1L
  if (chr_end < min(pos_chr)) {
    log_msg(sprintf("  Skipping %s — chromosome shorter than window size", chr))
    next
  }
  
  win_starts_chr <- seq(min(pos_chr), chr_end, by = WIN_STEP)
  
  if (length(win_starts_chr) == 0) {
    log_msg(sprintf("  Skipping %s — no valid windows", chr))
    next
  }
  
  log_msg(sprintf("Processing %s: %d SNPs, %d windows",
                  chr, length(pos_chr), length(win_starts_chr)))
  
  for (i in seq_along(win_starts_chr)) {
    
    w_start <- win_starts_chr[i]
    w_end   <- w_start + WIN_SIZE - 1L
    w_mid   <- (w_start + w_end) / 2
    
    # Which of this chromosome's SNPs fall in this window?
    # Then map back to dosage matrix column indices.
    in_chr_win <- which(pos_chr >= w_start & pos_chr <= w_end)
    in_win     <- chr_cols[in_chr_win]
    n_snps     <- length(in_win)
    
    dos_win <- if (n_snps >= MIN_N_SNPS) dosage_matrix[, in_win, drop = FALSE] else NULL
    
    
    # --- (a) Within-population: pi and tajima_d for each eligible group ------
    
    for (grp in eligible_groups) {
      ids     <- group_samples[[grp]]
      dos_grp <- if (!is.null(dos_win)) {
        dos_win[intersect(ids, rownames(dos_win)), , drop = FALSE]
      } else NULL
      
      pi_val <- if (!is.null(dos_grp) && nrow(dos_grp) >= 2) {
        compute_pi(dos_grp)
      } else NA_real_
      
      td_val <- if (!is.null(dos_grp) && nrow(dos_grp) >= 4) {
        compute_tajima(dos_grp)
      } else NA_real_
      
      within_buffer[[length(within_buffer) + 1]] <- tibble(
        chrom        = chr,
        group        = grp,
        window_start = w_start,
        window_end   = w_end,
        window_mid   = w_mid,
        n_snps       = n_snps,
        pi           = pi_val,
        tajima_d     = td_val
      )
    }
    
    
    # --- (b) Between-population: dxy and fst for each same-MSL-type pair ----
    
    for (pair in all_pairs) {
      grp1     <- pair[[1]]
      grp2     <- pair[[2]]
      ids1     <- group_samples[[grp1]]
      ids2     <- group_samples[[grp2]]
      msl_type <- if (grepl("_0$", grp1)) "MSL0-MSL0" else "MSL4-MSL4"
      
      dxy_val <- if (!is.null(dos_win)) compute_dxy(dos_win, ids1, ids2) else NA_real_
      fst_val <- if (!is.null(dos_win)) compute_fst(dos_win, ids1, ids2) else NA_real_
      
      between_buffer[[length(between_buffer) + 1]] <- tibble(
        chrom        = chr,
        pop1         = grp1,
        pop2         = grp2,
        msl_type     = msl_type,
        window_start = w_start,
        window_end   = w_end,
        window_mid   = w_mid,
        n_snps       = n_snps,
        dxy          = dxy_val,
        fst          = fst_val
      )
    }
    
    n_done     <- n_done + 1L
    total_wins <- total_wins + 1L
    
    # Flush buffers every BATCH_SIZE windows
    if (n_done %% BATCH_SIZE == 0L) {
      write_tsv(bind_rows(within_buffer),  within_path,  append = TRUE)
      write_tsv(bind_rows(between_buffer), between_path, append = TRUE)
      within_buffer  <- list()
      between_buffer <- list()
      n_done         <- 0L
      log_msg(sprintf("  Windows processed so far: %d (on %s)", total_wins, chr))
    }
  }
}

# Final flush for any remaining buffered rows
if (length(within_buffer) > 0 || length(between_buffer) > 0) {
  write_tsv(bind_rows(within_buffer),  within_path,  append = TRUE)
  write_tsv(bind_rows(between_buffer), between_path, append = TRUE)
}




# 12. CLEANUP AND SUMMARY ------------------------------------------------------

elapsed <- round((proc.time() - run_start)["elapsed"])

log_msg("====================================================")
log_msg("Run complete.")
log_msg(sprintf("Chromosomes processed: %d", length(chroms_available)))
log_msg(sprintf("Total windows:         %d", total_wins))
log_msg(sprintf("Groups per window:     %d (within-pop)", length(eligible_groups)))
log_msg(sprintf("Pairs per window:      %d (between-pop)", length(all_pairs)))
log_msg(sprintf("Within-pop rows:       %d", total_wins * length(eligible_groups)))
log_msg(sprintf("Between-pop rows:      %d", total_wins * length(all_pairs)))
log_msg(sprintf("Within output:         %s", within_path))
log_msg(sprintf("Between output:        %s", between_path))
log_msg(sprintf("Elapsed:               %d seconds (~%.1f minutes)", elapsed, elapsed / 60))
log_msg("====================================================")

sink()
close(log_con)