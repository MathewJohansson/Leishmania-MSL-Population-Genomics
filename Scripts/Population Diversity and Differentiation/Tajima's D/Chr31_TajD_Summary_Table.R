# Generate summary table of per-population Tajima's D results
# Summarises chr31 vs genome-wide median Tajima's D for MSL0 and MSL4




# 1. Setup ---------------------------------------------------------------------

library(tidyverse)

base_dir <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
data_dir <- file.path(base_dir, "Population Diversity and Differentiation Analysis/Per_population")
output_dir <- file.path(base_dir, "Tables")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# q_del values at >= 90% threshold (AM K13), V = AM notation confirmed
q_del <- c(
  V1=1.0000, V2=0.3533, V3=1.0000, V4=0.0417, V5=0.8036,
  V6=1.0000, V7=1.0000, V8=0.0500, V9=1.0000, V10=1.0000,
  V11=1.0000, V12=1.0000, V13=0.0000
)




# 2. Helper functions ----------------------------------------------------------

load_tajd <- function(path, msl_status) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(N_SNPS >= 10, !is.na(CHROM)) %>%
    mutate(MSL_status = msl_status)
}

# Get median Tajima's D for chr31 and genome-wide (excluding chr31)
summarise_tajd <- function(data) {
  list(
    chr31_median  = median(data$TajimaD[data$CHROM == "LinJ.31"], na.rm = TRUE),
    genome_median = median(data$TajimaD[data$CHROM != "LinJ.31"], na.rm = TRUE),
    chr31_n_windows = sum(data$CHROM == "LinJ.31"),
    chr31_is_most_negative = {
      per_chrom <- data %>%
        group_by(CHROM) %>%
        summarise(med = median(TajimaD, na.rm = TRUE), .groups = "drop")
      min_chrom <- per_chrom$CHROM[which.min(per_chrom$med)]
      min_chrom == "LinJ.31"
    }
  )
}




# 3. Main loop -----------------------------------------------------------------

populations <- paste0("V", 1:13)
results <- list()

for (pop in populations) {
  
  cat("Processing", pop, "...\n")
  
  f_msl0_chr31  <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL0_samples.diploidised.tsv.Tajima.D"))
  f_msl0_genome <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL0_samples.tsv.Tajima.D"))
  f_msl4_chr31  <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL4_samples.diploidised.tsv.Tajima.D"))
  f_msl4_genome <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL4_samples.tsv.Tajima.D"))
  
  has_msl0 <- file.exists(f_msl0_chr31) && file.exists(f_msl0_genome)
  has_msl4 <- file.exists(f_msl4_chr31) && file.exists(f_msl4_genome)
  
  row <- tibble(
    Population            = pop,
    AM_notation           = paste0("AM", str_extract(pop, "[0-9]+")),
    q_del                 = q_del[pop],
    MSL0_available        = has_msl0,
    MSL4_available        = has_msl4,
    MSL0_chr31_median     = NA_real_,
    MSL0_genome_median    = NA_real_,
    MSL0_chr31_n_windows  = NA_integer_,
    MSL0_chr31_most_neg   = NA,
    MSL4_chr31_median     = NA_real_,
    MSL4_genome_median    = NA_real_,
    MSL4_chr31_n_windows  = NA_integer_,
    MSL4_chr31_most_neg   = NA
  )
  
  if (has_msl0) {
    msl0 <- bind_rows(load_tajd(f_msl0_chr31, "MSL0"),
                      load_tajd(f_msl0_genome, "MSL0"))
    s <- summarise_tajd(msl0)
    row$MSL0_chr31_median    <- round(s$chr31_median, 3)
    row$MSL0_genome_median   <- round(s$genome_median, 3)
    row$MSL0_chr31_n_windows <- s$chr31_n_windows
    row$MSL0_chr31_most_neg  <- s$chr31_is_most_negative
  }
  
  if (has_msl4) {
    msl4 <- bind_rows(load_tajd(f_msl4_chr31, "MSL4"),
                      load_tajd(f_msl4_genome, "MSL4"))
    s <- summarise_tajd(msl4)
    row$MSL4_chr31_median    <- round(s$chr31_median, 3)
    row$MSL4_genome_median   <- round(s$genome_median, 3)
    row$MSL4_chr31_n_windows <- s$chr31_n_windows
    row$MSL4_chr31_most_neg  <- s$chr31_is_most_negative
  }
  
  results[[pop]] <- row
}

summary_table <- bind_rows(results)




# 4. Print and save ------------------------------------------------------------

cat("\n=== SUMMARY TABLE ===\n")
print(summary_table, n = Inf)

# Save as CSV
out_path <- file.path(output_dir, "Chr31_TajD_Per_Population_Summary.csv")
write.csv(summary_table, out_path, row.names = FALSE)
cat("\nSaved to:", out_path, "\n")

# Print key findings
cat("\n=== KEY FINDINGS ===\n")
cat("Populations where chr31 is most negative (MSL0):\n")
print(summary_table %>%
        filter(MSL0_available, MSL0_chr31_most_neg) %>%
        select(Population, q_del, MSL0_chr31_median, MSL0_genome_median))

cat("\nPopulations where chr31 is most negative (MSL4):\n")
print(summary_table %>%
        filter(MSL4_available, MSL4_chr31_most_neg) %>%
        select(Population, q_del, MSL4_chr31_median, MSL4_genome_median))

cat("\nPopulations with both MSL0 and MSL4 (direct comparison possible):\n")
print(summary_table %>%
        filter(MSL0_available, MSL4_available) %>%
        select(Population, q_del,
               MSL0_chr31_median, MSL0_genome_median,
               MSL4_chr31_median, MSL4_genome_median))

