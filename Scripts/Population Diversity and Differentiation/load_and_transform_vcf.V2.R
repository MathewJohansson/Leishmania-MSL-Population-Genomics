# SUMMARY -----------------------------------------------------------------------
# Purpose: Load VCF files, extract dosage data, and save cleaned datasets
# for StAMPP, adegenet (genlight), and other downstream tools.
#
# Note on outgroups: These VCF files contain European samples and African
# outgroups. They are retained in this base dataset as they are useful for
# rooting phylogenies and global analyses. If an Americas-only focus is needed,
# subset the data in the downstream analysis scripts.
#
# Output files (saved to Data/Fst/):
#   cleaned_tetraploid_vcf_data_2026-04-13.RData
#   cleaned_diploid_vcf_data_2026-04-13.RData




## 1. CLEAN UP AND LOAD LIBRARIES ----------------------------------------------

rm(list = ls())
#setwd("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript")
library(tidyverse)
library(vcfR)
library(adegenet)
library(StAMPP)




## 2. LOAD SAMPLE METADATA -----------------------------------------------------

# VCF files used in this script:
# Chr31, tetraploid calls: Data/tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.vcf.gz
# All chromosomes, diploid calls: Data/Linfantum_samples_with_outgroup_biallelic.2025-06-10.vcf.gz

# Column 16 (population) is renamed to Americas_population here.
# V-notation (V1-V13) is remapped to AM-notation (AM1-AM13).
samples <- read_tsv("sample_info/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
names(samples)[1] <- "sra_run_accession"
names(samples)[16] <- "Americas_population"
names(samples)[17] <- "Americas_population_proportion"

samples <- samples |>
  mutate(Americas_population = replace_na(Americas_population, "Non_Americas")) |>
  mutate(Americas_population = str_replace(Americas_population, "^V", "AM"))

unique(samples$Americas_population)




# PART 1: TETRAPLOID CHR31 VCF PROCESSING


## 3. LOAD TETRAPLOID VCF ------------------------------------------------------

# This VCF contains all samples, chromosome 31 only, with tetraploid calls.
vcf.tet.chr31 <- read.vcfR("Data/tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.vcf.gz")

# Alternative VCFs available locally if needed:
# vcf <- read.vcfR("Data/tetraploid.snps.filter.Americas_only_MSL0and4_samples.vcf.gz")  # Americas MSL0/MSL4 only
# vcf <- read.vcfR("Data/tetraploid_Americas_LinJ.31_filtered_2026-04-07.vcf.gz")         # All Americas samples

# To create further subsets on Viking:
# 1. Make a samples list:
#    grep Americas Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv | cut -f 1 > Americas.ids
# 2. Subset the VCF with bcftools:
#    bcftools view -S Americas.ids -O z -o tetraploid_Americas_LinJ.31_filtered_2026-04-07.vcf.gz \
#    tetraploid_freebayes_regenotyped_sorted_LinJ.31_2025-07-01.filtered_2026-04-07.vcf.gz
# 3. Index the VCF:
#    bcftools index -t tetraploid_Americas_LinJ.31_filtered_2026-04-07.vcf.gz




## 4. PROCESS TETRAPLOID VCF FOR DOSAGE DATA -----------------------------------

# Filter out the MSL region (chr31:1,181,282-1,192,406) to prevent structural/selection bias
chrom <- getCHROM(vcf.tet.chr31)
pos <- as.numeric(getPOS(vcf.tet.chr31))

keep_snps <- !(chrom == "LinJ.31" & pos >= 1181282 & pos <= 1192406)
vcf.no.msl <- vcf.tet.chr31[keep_snps, ]

# Extract raw genotype strings (e.g., "0/0/1/1")
gt_strings <- extract.gt(vcf.no.msl, element = "GT")

# Function to count the alt (1) alleles in each genotype string
count_alt_alleles <- function(x) {
  lengths(regmatches(x, gregexpr("1", x)))
}

# Apply across the matrix and restore row/column names
dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings),
  ncol = ncol(gt_strings)
)
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)

# Set missing data back to NA
dosage_matrix[is.na(gt_strings) | gt_strings == "./." | gt_strings == "./././."] <- NA

# Transpose so samples are rows
dosage_matrix <- t(dosage_matrix)




## 5. IMPUTE MISSING DOSAGES (FOR STAMPP/PCA ONLY) -----------------------------

# PCA and StAMPP require complete data.
# Do NOT use imputed dosages for phylogenetic analyses — missing data is handled
# natively and imputation introduces bias.

dosage_imputed <- apply(dosage_matrix, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})




## 6. REMOVE SITES WITH HIGH MISSING DATA RATES --------------------------------

# Sites with >20% missing data are flagged here.
# This filtered object (dosage_filtered_nas) is not used downstream in this
# script but is retained for inspection.
missing_prop <- colMeans(is.na(dosage_matrix))
keep_snps <- missing_prop <= 0.20
dosage_filtered_nas <- dosage_matrix[, keep_snps]

print(paste("Dropped SNPs with >20% missing data. Remaining:", ncol(dosage_filtered_nas)))




## 7. OUTLIER MANAGEMENT -------------------------------------------------------

# African/European outgroups are retained in this base RData file.
# Outlier removal is better handled within each downstream analysis script.
# Add sample IDs to outliers_to_remove below only if they should be excluded
# from ALL analyses globally.

# outliers_to_remove <- c("ERR205787", "SRR8608748")

outliers_to_remove <- c()  # Empty by default

if (length(outliers_to_remove) > 0) {
  dosage_imputed <- dosage_imputed[!grepl(paste(outliers_to_remove, collapse = "|"), rownames(dosage_imputed)), ]
}
print(paste("Sample count after filtering specified outliers:", nrow(dosage_imputed)))




## 8. REMOVE INVARIANT SITES ---------------------------------------------------

column_vars <- apply(dosage_imputed, 2, var, na.rm = TRUE)
keep_cols <- !is.na(column_vars) & column_vars > 0
dosage_filtered <- dosage_imputed[, keep_cols]

print(paste("TETRAPLOID - Total SNPs initially:", ncol(dosage_imputed)))
print(paste("TETRAPLOID - Remaining SNPs for StAMPP & Genlight:", ncol(dosage_filtered)))




## 9. PREPARE METADATA ---------------------------------------------------------

meta_combined <- data.frame(sra_run_accession = rownames(dosage_filtered)) |>
  left_join(samples, by = "sra_run_accession")




## 10. PREPARE STAMPP FORMAT ---------------------------------------------------

allele_freqs <- dosage_filtered / 4

stampp_df <- as.data.frame(allele_freqs)
stampp_df$Sample <- rownames(stampp_df)

stampp_df <- stampp_df |>
  left_join(meta_combined |> select(sra_run_accession, Americas_population),
            by = c("Sample" = "sra_run_accession"))

stampp_final <- stampp_df |>
  rename(Pop = Americas_population) |>
  mutate(Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs)))

stampp_obj <- stamppConvert(stampp_final, type = "r")
print("Tetraploid StAMPP object created successfully!")




## 11. PREPARE GENLIGHT FORMAT -------------------------------------------------

genlight_obj <- new("genlight", dosage_filtered)
ploidy(genlight_obj) <- 4
pop(genlight_obj) <- meta_combined$Americas_population

print("Tetraploid Genlight object created successfully!")




## 12. SAVE TETRAPLOID WORKSPACE -----------------------------------------------

save(stampp_obj, genlight_obj, meta_combined, dosage_filtered,
     file = "Data/Fst/cleaned_tetraploid_vcf_data_2026-04-13.RData")

print("Tetraploid RData saved to Data/Fst/cleaned_tetraploid_vcf_data_2026-04-13.RData")




# PART 2: DIPLOID WHOLE-GENOME VCF PROCESSING


## 13. LOAD DIPLOID VCF --------------------------------------------------------

vcf.dip <- read.vcfR("Data/Linfantum_samples_with_outgroup_biallelic.2025-06-10.vcf.gz")




## 14. PROCESS DIPLOID VCF FOR DOSAGE DATA -------------------------------------

# Filter out the MSL region
chrom_dip <- getCHROM(vcf.dip)
pos_dip <- as.numeric(getPOS(vcf.dip))

keep_snps_dip <- !(chrom_dip == "LinJ.31" & pos_dip >= 1181282 & pos_dip <= 1192406)
vcf.no.msl_dip <- vcf.dip[keep_snps_dip, ]

# Extract raw genotype strings (e.g., "0/1")
gt_strings_dip <- extract.gt(vcf.no.msl_dip, element = "GT")

dosage_matrix_dip <- matrix(
  count_alt_alleles(gt_strings_dip),
  nrow = nrow(gt_strings_dip),
  ncol = ncol(gt_strings_dip)
)
rownames(dosage_matrix_dip) <- rownames(gt_strings_dip)
colnames(dosage_matrix_dip) <- colnames(gt_strings_dip)

# Set missing data back to NA
dosage_matrix_dip[is.na(gt_strings_dip) | gt_strings_dip == "./."] <- NA

# Transpose so samples are rows
dosage_matrix_dip <- t(dosage_matrix_dip)

# Impute missing values with per-SNP mean
dosage_imputed_dip <- apply(dosage_matrix_dip, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})




## 15. OUTLIER MANAGEMENT (DIPLOID) --------------------------------------------

# ERR205787 (Panama) and SRR8608748 (Honduras) are outliers in Americas analyses.
outliers_to_remove_dip <- c("ERR205787", "SRR8608748")

if (length(outliers_to_remove_dip) > 0) {
  dosage_imputed_dip <- dosage_imputed_dip[!grepl(paste(outliers_to_remove_dip, collapse = "|"), rownames(dosage_imputed_dip)), ]
}
print(paste("Removed specific outliers. New diploid sample count:", nrow(dosage_imputed_dip)))




## 16. REMOVE INVARIANT SITES (DIPLOID) ----------------------------------------

column_vars_dip <- apply(dosage_imputed_dip, 2, var, na.rm = TRUE)
keep_cols_dip <- !is.na(column_vars_dip) & column_vars_dip > 0
dosage_filtered_dip <- dosage_imputed_dip[, keep_cols_dip]

print(paste("DIPLOID - Total SNPs initially:", ncol(dosage_imputed_dip)))
print(paste("DIPLOID - Remaining SNPs for StAMPP & Genlight:", ncol(dosage_filtered_dip)))




## 17. PREPARE METADATA (DIPLOID) ----------------------------------------------

meta_combined_dip <- data.frame(sra_run_accession = rownames(dosage_filtered_dip)) |>
  left_join(samples, by = "sra_run_accession")




## 18. PREPARE STAMPP FORMAT (DIPLOID) -----------------------------------------

# Divide by 2 for diploidy
allele_freqs_dip <- dosage_filtered_dip / 2

stampp_df_dip <- as.data.frame(allele_freqs_dip)
stampp_df_dip$Sample <- rownames(stampp_df_dip)

stampp_df_dip <- stampp_df_dip |>
  left_join(meta_combined_dip |> select(sra_run_accession, Americas_population),
            by = c("Sample" = "sra_run_accession"))

stampp_final_dip <- stampp_df_dip |>
  rename(Pop = Americas_population) |>
  mutate(Ploidy = 2, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs_dip)))

stampp_obj_dip <- stamppConvert(stampp_final_dip, type = "r")
print("Diploid StAMPP object created successfully!")




## 19. PREPARE GENLIGHT FORMAT (DIPLOID) ---------------------------------------

genlight_obj_dip <- new("genlight", dosage_filtered_dip)
ploidy(genlight_obj_dip) <- 2
pop(genlight_obj_dip) <- meta_combined_dip$Americas_population

print("Diploid Genlight object created successfully!")




## 20. SAVE DIPLOID WORKSPACE --------------------------------------------------

save(stampp_obj_dip, genlight_obj_dip, meta_combined_dip, dosage_filtered_dip,
     file = "Data/Fst/cleaned_diploid_vcf_data_2026-04-13.RData")

print("Diploid RData saved to Data/Fst/cleaned_diploid_vcf_data_2026-04-13.RData")
print("All done! Both tetraploid and diploid datasets formatted and saved.")



