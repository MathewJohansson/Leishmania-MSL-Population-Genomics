

# PURPOSE: Verify R implementations of π, Tajima's D, and Fst against
#          VCFtools reference outputs provided by Dan (2026-06-04).

# SCOPE:   - Two chromosomes: LinJ.36 (diploid) and LinJ.31 (pseudodiploid)
#          - Two populations: AM2 (n=65, all MSL types) and AM9 (n=20)
#          - Windows: 50kb non-overlapping (matching VCFtools exactly)
#          - NOT production logic: no MSL subsetting, no full chromosome loop

# NOTE:    Window step = window size here (non-overlapping). Production script
#          uses 5kb step — do not conflate these results.

# OUTPUT:  Console summary (Pearson r, mean difference per stat)
#          Data/Unified_Popgen/Verification/verification_results.tsv



# 1. CONFIG — all paths and constants defined here, nothing hardcoded below ----

library(vcfR)

BASE <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

# VCF files
VCF_DIPLOID     <- file.path(BASE, "Data",
                             "Linfantum_samples_no_outgroup_biallelic.2025-06-10.noLinJ.31.SNPsonly.vcf.gz")
VCF_PSEUDODIP   <- file.path(BASE, "Data",
                             "tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.diploidised.vcf.gz")

# Verification directory
VERIF_DIR <- file.path(BASE, "Data", "Unified_Popgen", "Verification")

# Dan's reference output files
REF <- list(
  pi_AM2_dip      = file.path(VERIF_DIR, "AM2.samples.txt.windowed.pi"),
  pi_AM9_dip      = file.path(VERIF_DIR, "AM9.samples.txt.windowed.pi"),
  pi_AM2_pdip     = file.path(VERIF_DIR, "AM2.samples.diploidised.txt.windowed.pi"),
  pi_AM9_pdip     = file.path(VERIF_DIR, "AM9.samples.diploidised.txt.windowed.pi"),
  tajd_AM2_dip    = file.path(VERIF_DIR, "AM2.samples.txt.Tajima.D"),
  tajd_AM9_dip    = file.path(VERIF_DIR, "AM9.samples.txt.Tajima.D"),
  tajd_AM2_pdip   = file.path(VERIF_DIR, "AM2.samples.diploidised.txt.Tajima.D"),
  tajd_AM9_pdip   = file.path(VERIF_DIR, "AM9.samples.diploidised.txt.Tajima.D"),
  fst_dip         = file.path(VERIF_DIR, "AM2-AM9.samples.txt.windowed.weir.fst"),
  fst_pdip        = file.path(VERIF_DIR, "AM2-AM9.samples.diploidised.txt.windowed.weir.fst")
)

# Sample list files (column 1 = sample ID)
SAMPLES <- list(
  AM2_dip    = file.path(VERIF_DIR, "AM2.samples.txt"),
  AM9_dip    = file.path(VERIF_DIR, "AM9.samples.txt"),
  AM2_pdip   = file.path(VERIF_DIR, "AM2.samples.diploidised.txt"),
  AM9_pdip   = file.path(VERIF_DIR, "AM9.samples.diploidised.txt")
)

# Window parameters — 50kb non-overlapping to match VCFtools reference
WIN_SIZE <- 50000L
STEP     <- 50000L   # non-overlapping; production uses 5000L

# Ploidy
PLOIDY_DIP   <- 2L   # diploid non-chr31
PLOIDY_PDIP  <- 2L   # pseudodiploid chr31 (tetraploid split into 2 pseudo-haplotypes)

# Chromosomes
CHR_DIPLOID  <- "LinJ.36"
CHR_PSEUDODIP <- "LinJ.31"

# MSL locus on chr31 — excluded from all analyses (hard requirement)
MSL_START <- 1181282L
MSL_END   <- 1192406L



# 2. FORMULA FUNCTIONS
# These are the confirmed production formulas — do not modify here.
# n_gene = n_samples × ploidy throughout.

compute_pi <- function(dos_win, ploidy) {
  # Per-base π: sum of per-SNP heterozygosity divided by window length
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

compute_tajima <- function(dos_win, ploidy) {
  # Tajima's D using gene copies (n_samples × ploidy), not n_samples
  # Correlation with VCFtools ~0.79 on previous 5kb-step run;
  # 50kb non-overlapping windows used here for direct comparison
  n_samples <- nrow(dos_win)
  n_gene    <- n_samples * ploidy
  
  per_snp <- apply(dos_win, 2, function(x) {
    x    <- x[!is.na(x)]
    n_x  <- length(x)
    if (n_x < 2) return(c(seg = 0L, k = 0))
    p <- mean(x) / ploidy
    if (p <= 0 || p >= 1) return(c(seg = 0L, k = 0))
    n_gene_x <- n_x * ploidy
    c(seg = 1L, k = (n_gene_x / (n_gene_x - 1)) * 2 * p * (1 - p))
  })
  
  S <- sum(per_snp["seg", ])
  if (S < 2) return(NA_real_)
  
  k  <- sum(per_snp["k", ])
  a1 <- sum(1 / seq_len(n_gene - 1))
  a2 <- sum(1 / seq_len(n_gene - 1)^2)
  b1 <- (n_gene + 1) / (3 * (n_gene - 1))
  b2 <- 2 * (n_gene^2 + n_gene + 3) / (9 * n_gene * (n_gene - 1))
  c1 <- b1 - 1 / a1
  c2 <- b2 - (n_gene + 2) / (a1 * n_gene) + a2 / a1^2
  e1 <- c1 / a1
  e2 <- c2 / (a1^2 + a2)
  
  var_d <- e1 * S + e2 * S * (S - 1)
  if (var_d <= 0) return(NA_real_)
  (k - S / a1) / sqrt(var_d)
}

compute_fst <- function(dos_win, ids1, ids2, ploidy, MIN_FST_SNPS = 10L) {
  # Hudson's estimator: per-SNP Fst averaged across window
  # Corresponds to VCFtools MEAN_FST (not WEIGHTED_FST)
  p1 <- intersect(ids1, rownames(dos_win))
  p2 <- intersect(ids2, rownames(dos_win))
  if (length(p1) < 5 || length(p2) < 5) return(NA_real_)
  
  mat1 <- dos_win[p1, , drop = FALSE]
  mat2 <- dos_win[p2, , drop = FALSE]
  
  per_snp <- sapply(seq_len(ncol(dos_win)), function(j) {
    a <- mat1[, j] / ploidy
    b <- mat2[, j] / ploidy
    a <- a[!is.na(a)]
    b <- b[!is.na(b)]
    if (length(a) < 5 || length(b) < 5) return(NA_real_)
    f1 <- mean(a); f2 <- mean(b)
    n1 <- length(a); n2 <- length(b)
    pi1 <- 2 * f1 * (1 - f1) * n1 / (n1 - 1)
    pi2 <- 2 * f2 * (1 - f2) * n2 / (n2 - 1)
    pi_w <- (pi1 + pi2) / 2
    pi_b <- f1 * (1 - f2) + f2 * (1 - f1)
    if (pi_b == 0) return(NA_real_)
    (pi_b - pi_w) / pi_b
  })
  
  if (sum(!is.na(per_snp)) < MIN_FST_SNPS) return(NA_real_)
  max(0, min(1, mean(per_snp, na.rm = TRUE)))
}



# 3. VCF LOADING
# Reads a single chromosome from VCF, subsets to target samples,
# removes MSL locus if chr31, returns dosage matrix (samples × SNPs).

load_dosage <- function(vcf_path, chr, sample_ids, ploidy) {
  cat(sprintf("  Loading %s from %s ...\n", chr, basename(vcf_path)))
  
  vcf <- read.vcfR(vcf_path, verbose = FALSE)
  
  # Subset to target chromosome
  keep_chr <- vcf@fix[, "CHROM"] == chr
  vcf <- vcf[keep_chr, ]
  
  # Remove MSL locus if chr31 (hard requirement)
  if (chr == "LinJ.31") {
    pos     <- as.integer(vcf@fix[, "POS"])
    in_msl  <- pos >= MSL_START & pos <= MSL_END
    n_removed <- sum(in_msl)
    vcf     <- vcf[!in_msl, ]
    cat(sprintf("  Removed %d SNPs in MSL locus (chr31:%d-%d)\n",
                n_removed, MSL_START, MSL_END))
  }
  
  # Extract genotype matrix and convert to dosage
  gt  <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  pos <- as.integer(vcf@fix[, "POS"])
  
  # Convert GT strings to dosage (sum of alleles)
  # Handles both phased (|) and unphased (/) separators
  gt_to_dos <- function(g) {
    if (is.na(g) || g == "./." || g == ".|.") return(NA_real_)
    alleles <- as.integer(strsplit(g, "[/|]")[[1]])
    if (any(is.na(alleles))) return(NA_real_)
    sum(alleles)
  }
  
  dos <- apply(gt, c(1, 2), gt_to_dos)
  # dos is SNPs × samples; transpose to samples × SNPs
  dos <- t(dos)
  
  # Subset to requested sample IDs
  present <- intersect(sample_ids, rownames(dos))
  missing <- setdiff(sample_ids, rownames(dos))
  if (length(missing) > 0) {
    warning(sprintf("  %d requested samples not found in VCF: %s",
                    length(missing), paste(head(missing, 5), collapse = ", ")))
  }
  dos <- dos[present, , drop = FALSE]
  colnames(dos) <- pos
  
  cat(sprintf("  Loaded: %d samples, %d SNPs\n", nrow(dos), ncol(dos)))
  dos
}



# 4. WINDOWED COMPUTATION
# Slides 50kb non-overlapping windows across the dosage matrix,
# computing π, Tajima's D, and optionally Fst per window.
# Returns a data.frame with one row per window.

compute_windows <- function(dos, ids1, ids2 = NULL, ploidy,
                            chr, compute_fst_flag = FALSE) {
  pos      <- as.integer(colnames(dos))
  chr_max  <- max(pos)
  starts   <- seq(1L, chr_max, by = STEP)
  
  rows <- lapply(starts, function(ws) {
    we      <- ws + WIN_SIZE - 1L
    in_win  <- pos >= ws & pos <= we
    n_snps  <- sum(in_win)
    
    if (n_snps == 0) {
      row <- data.frame(chrom = chr, win_start = ws, win_end = we,
                        n_snps = 0L, pi = NA_real_, tajima_d = NA_real_)
      if (compute_fst_flag) row$fst <- NA_real_
      return(row)
    }
    
    dos_win <- dos[, in_win, drop = FALSE]
    
    pi_val   <- compute_pi(dos_win, ploidy)
    tajd_val <- compute_tajima(dos_win, ploidy)
    row <- data.frame(chrom = chr, win_start = ws, win_end = we,
                      n_snps = n_snps, pi = pi_val, tajima_d = tajd_val)
    
    if (compute_fst_flag) {
      row$fst <- compute_fst(dos_win, ids1, ids2, ploidy)
    }
    row
  })
  
  do.call(rbind, rows)
}



# 5. COMPARISON FUNCTIONS
# Joins R output to VCFtools reference by window position,
# reports Pearson r and mean absolute difference.

compare_pi <- function(r_df, ref_path, label) {
  ref <- read.table(ref_path, header = TRUE, sep = "\t")
  # VCFtools BIN_START is 1-based; our win_start is also 1-based — join on start
  merged <- merge(r_df, ref, by.x = "win_start", by.y = "BIN_START")
  merged <- merged[!is.na(merged$pi) & !is.na(merged$PI), ]
  
  r   <- cor(merged$pi, merged$PI, use = "complete.obs")
  mad <- mean(abs(merged$pi - merged$PI))
  
  cat(sprintf("\n[π] %s\n", label))
  cat(sprintf("  Windows compared: %d\n", nrow(merged)))
  cat(sprintf("  Pearson r:        %.4f\n", r))
  cat(sprintf("  Mean |diff|:      %.3e\n", mad))
  
  merged$stat        <- "pi"
  merged$label       <- label
  merged$r_value     <- merged$pi
  merged$vcft_value  <- merged$PI
  merged$diff        <- merged$pi - merged$PI
  merged[, c("label", "stat", "win_start", "n_snps", "r_value", "vcft_value", "diff")]
}

compare_tajd <- function(r_df, ref_path, label) {
  ref <- read.table(ref_path, header = TRUE, sep = "\t")
  # VCFtools TajimaD BIN_START is 0-based; add 1 to align with our 1-based starts
  ref$BIN_START_1 <- ref$BIN_START + 1L
  merged <- merge(r_df, ref, by.x = "win_start", by.y = "BIN_START_1")
  merged <- merged[!is.na(merged$tajima_d) & !is.na(merged$TajimaD), ]
  # Exclude VCFtools nan rows (0-SNP windows)
  merged <- merged[merged$TajimaD != "nan", ]
  merged$TajimaD <- as.numeric(merged$TajimaD)
  
  r   <- cor(merged$tajima_d, merged$TajimaD, use = "complete.obs")
  mad <- mean(abs(merged$tajima_d - merged$TajimaD))
  
  cat(sprintf("\n[Tajima's D] %s\n", label))
  cat(sprintf("  Windows compared: %d\n", nrow(merged)))
  cat(sprintf("  Pearson r:        %.4f\n", r))
  cat(sprintf("  Mean |diff|:      %.3e\n", mad))
  
  # S count comparison
  s_match <- sum(merged$n_snps == merged$N_SNPS)
  cat(sprintf("  S count match:    %d / %d windows\n", s_match, nrow(merged)))
  
  merged$stat        <- "tajima_d"
  merged$label       <- label
  merged$r_value     <- merged$tajima_d
  merged$vcft_value  <- merged$TajimaD
  merged$diff        <- merged$tajima_d - merged$TajimaD
  merged[, c("label", "stat", "win_start", "n_snps", "r_value", "vcft_value", "diff")]
}

compare_fst <- function(r_df, ref_path, label) {
  ref <- read.table(ref_path, header = TRUE, sep = "\t")
  merged <- merge(r_df, ref, by.x = "win_start", by.y = "BIN_START")
  merged <- merged[!is.na(merged$fst) & !is.na(merged$MEAN_FST), ]
  
  # Compare against both MEAN_FST and WEIGHTED_FST
  r_mean     <- cor(merged$fst, merged$MEAN_FST,     use = "complete.obs")
  r_weighted <- cor(merged$fst, merged$WEIGHTED_FST, use = "complete.obs")
  mad_mean     <- mean(abs(merged$fst - merged$MEAN_FST))
  mad_weighted <- mean(abs(merged$fst - merged$WEIGHTED_FST))
  
  cat(sprintf("\n[Fst] %s\n", label))
  cat(sprintf("  Windows compared:          %d\n", nrow(merged)))
  cat(sprintf("  vs MEAN_FST     — r: %.4f, mean |diff|: %.4f\n", r_mean,     mad_mean))
  cat(sprintf("  vs WEIGHTED_FST — r: %.4f, mean |diff|: %.4f\n", r_weighted, mad_weighted))
  
  # Report which is closer
  better <- if (r_mean >= r_weighted) "MEAN_FST" else "WEIGHTED_FST"
  cat(sprintf("  Better match:              %s\n", better))
  
  merged$stat        <- "fst"
  merged$label       <- label
  merged$r_value     <- merged$fst
  merged$vcft_value  <- merged$MEAN_FST   # primary comparison
  merged$diff        <- merged$fst - merged$MEAN_FST
  merged[, c("label", "stat", "win_start", "n_snps", "r_value", "vcft_value", "diff")]
}



# 6. EXECUTION
# Load VCFs, compute windows, compare against reference outputs.

cat("=============================================================\n")
cat("Verification: R formulas vs VCFtools reference (Dan 2026-06-04)\n")
cat("Windows: 50kb non-overlapping (NOT production 5kb step)\n")
cat("=============================================================\n")

# --- Read sample ID lists (column 1 only) ---
ids_AM2_dip  <- read.table(SAMPLES$AM2_dip,  header = FALSE, sep = "\t", 
                           fill = TRUE, quote = "")[[1]]
ids_AM9_dip  <- read.table(SAMPLES$AM9_dip,  header = FALSE, sep = "\t", 
                           fill = TRUE, quote = "")[[1]]
ids_AM2_pdip <- read.table(SAMPLES$AM2_pdip, header = FALSE, sep = "\t", 
                           fill = TRUE, quote = "")[[1]]
ids_AM9_pdip <- read.table(SAMPLES$AM9_pdip, header = FALSE, sep = "\t", 
                           fill = TRUE, quote = "")[[1]]

all_results <- list()

# --- LinJ.36 (diploid) ---
cat("\n--- LinJ.36 (diploid) ---\n")
dos36_AM2 <- load_dosage(VCF_DIPLOID, CHR_DIPLOID, ids_AM2_dip, PLOIDY_DIP)
dos36_AM9 <- load_dosage(VCF_DIPLOID, CHR_DIPLOID, ids_AM9_dip, PLOIDY_DIP)

# Compute windows for each population; Fst computed on combined matrix
dos36_both <- load_dosage(VCF_DIPLOID, CHR_DIPLOID,
                          c(ids_AM2_dip, ids_AM9_dip), PLOIDY_DIP)

win36_AM2  <- compute_windows(dos36_AM2,  ids_AM2_dip, ploidy = PLOIDY_DIP, chr = CHR_DIPLOID)
win36_AM9  <- compute_windows(dos36_AM9,  ids_AM9_dip, ploidy = PLOIDY_DIP, chr = CHR_DIPLOID)
win36_fst  <- compute_windows(dos36_both, ids_AM2_dip, ids_AM9_dip,
                              ploidy = PLOIDY_DIP, chr = CHR_DIPLOID,
                              compute_fst_flag = TRUE)

all_results[["pi_AM2_dip"]]   <- compare_pi(win36_AM2,  REF$pi_AM2_dip,   "AM2 LinJ.36 diploid")
all_results[["pi_AM9_dip"]]   <- compare_pi(win36_AM9,  REF$pi_AM9_dip,   "AM9 LinJ.36 diploid")
all_results[["tajd_AM2_dip"]] <- compare_tajd(win36_AM2, REF$tajd_AM2_dip, "AM2 LinJ.36 diploid")
all_results[["tajd_AM9_dip"]] <- compare_tajd(win36_AM9, REF$tajd_AM9_dip, "AM9 LinJ.36 diploid")
all_results[["fst_dip"]]      <- compare_fst(win36_fst,  REF$fst_dip,      "AM2-AM9 LinJ.36 diploid")

# --- LinJ.31 (pseudodiploid chr31) ---
cat("\n--- LinJ.31 (pseudodiploid chr31) ---\n")
dos31_AM2 <- load_dosage(VCF_PSEUDODIP, CHR_PSEUDODIP, ids_AM2_pdip, PLOIDY_PDIP)
dos31_AM9 <- load_dosage(VCF_PSEUDODIP, CHR_PSEUDODIP, ids_AM9_pdip, PLOIDY_PDIP)

dos31_both <- load_dosage(VCF_PSEUDODIP, CHR_PSEUDODIP,
                          c(ids_AM2_pdip, ids_AM9_pdip), PLOIDY_PDIP)

win31_AM2  <- compute_windows(dos31_AM2,  ids_AM2_pdip, ploidy = PLOIDY_PDIP, chr = CHR_PSEUDODIP)
win31_AM9  <- compute_windows(dos31_AM9,  ids_AM9_pdip, ploidy = PLOIDY_PDIP, chr = CHR_PSEUDODIP)
win31_fst  <- compute_windows(dos31_both, ids_AM2_pdip, ids_AM9_pdip,
                              ploidy = PLOIDY_PDIP, chr = CHR_PSEUDODIP,
                              compute_fst_flag = TRUE)

all_results[["pi_AM2_pdip"]]   <- compare_pi(win31_AM2,  REF$pi_AM2_pdip,   "AM2 LinJ.31 pseudodiploid")
all_results[["pi_AM9_pdip"]]   <- compare_pi(win31_AM9,  REF$pi_AM9_pdip,   "AM9 LinJ.31 pseudodiploid")
all_results[["tajd_AM2_pdip"]] <- compare_tajd(win31_AM2, REF$tajd_AM2_pdip, "AM2 LinJ.31 pseudodiploid")
all_results[["tajd_AM9_pdip"]] <- compare_tajd(win31_AM9, REF$tajd_AM9_pdip, "AM9 LinJ.31 pseudodiploid")
all_results[["fst_pdip"]]      <- compare_fst(win31_fst,  REF$fst_pdip,      "AM2-AM9 LinJ.31 pseudodiploid")



# 7. OUTPUT

cat("\n=============================================================\n")
cat("Writing verification_results.tsv ...\n")

combined <- do.call(rbind, all_results)
rownames(combined) <- NULL

out_path <- file.path(VERIF_DIR, "verification_results.tsv")
write.table(combined, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("Written: %s\n", out_path))
cat("Done.\n")



