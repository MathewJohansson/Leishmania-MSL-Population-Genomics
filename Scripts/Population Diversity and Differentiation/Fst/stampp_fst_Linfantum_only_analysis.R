# SUMMARY -----------------------------------------------------------------------
# FST analysis — L. infantum only (LI subset, K=13)
# Adapted from Dan's stampp_fst_Americas_analysis.R
# Uses tetraploid chr31 calls + diploid whole-genome calls
# Populations: LI1-LI13 (excludes 2 L. donovani outgroup samples)
# Groups by LI population x MSL status (0 or 4)
#
# Input RData files (created by load_and_transform_vcf_V2.R):
#   Data/Fst/cleaned_tetraploid_vcf_data_2026-04-13.RData
#   Data/Fst/cleaned_diploid_vcf_data_2026-04-13.RData




## 1. CLEAN UP AND LOAD LIBRARIES -----------------------------------------------

rm(list = ls())
setwd("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript")
library(vcfR)
library(adegenet)
library(tidyverse)
library(StAMPP)
library(patchwork)




## 2. LOAD DATA -----------------------------------------------------------------

# Chr31 tetraploid data
load("Data/Fst/cleaned_tetraploid_vcf_data_2026-04-13.RData")
# Objects loaded: dosage_filtered, genlight_obj, meta_combined, stampp_obj

# Whole-genome diploid data
load("Data/Fst/cleaned_diploid_vcf_data_2026-04-13.RData")
# Objects loaded: dosage_filtered_dip, genlight_obj_dip, meta_combined_dip, stampp_obj_dip




## 3. PREPARE METADATA FOR FST --------------------------------------------------

# Load LI population assignments and join onto meta_combined.
# The membership table is the single source of truth for LI population labels.
li_pops <- read_tsv("Data/population_structure/population_membership_table_infantum_only_k13_LI_notation.tsv") |>
  select(sample, population) |>
  rename(sra_run_accession = sample, LI_population = population)

meta_combined <- meta_combined |>
  left_join(li_pops, by = "sra_run_accession") |>
  mutate(
    # Samples not in the LI membership table (e.g. L. donovani outgroup) get Non_LI
    LI_population = replace_na(LI_population, "Non_LI"),
    MSL_status = ifelse(is.na(MSL_coverage_rounded), "Unknown", as.character(MSL_coverage_rounded)),
    LI_population_MSL = paste(LI_population, MSL_status, sep = "_")
  )

table(meta_combined$LI_population_MSL)
unique(meta_combined$LI_population_MSL)




## 4. FILTER DATA FOR FST ANALYSIS ----------------------------------------------

# Keep only LI populations and MSL 0 or 4
meta_fst <- meta_combined |>
  filter(LI_population != "Non_LI") |>
  filter(MSL_status %in% c("0", "4"))

print("Populations remaining for FST:")
table(meta_fst$LI_population_MSL)

samples_to_keep <- meta_fst$sra_run_accession
dosage_fst <- dosage_filtered[rownames(dosage_filtered) %in% samples_to_keep, ]

col_vars_fst <- apply(dosage_fst, 2, var, na.rm = TRUE)
keep_cols_fst <- !is.na(col_vars_fst) & col_vars_fst > 0
dosage_fst <- dosage_fst[, keep_cols_fst]

print(paste("Samples remaining for FST:", nrow(dosage_fst)))
print(paste("SNPs remaining for FST:", ncol(dosage_fst)))




## 5. REBUILD STAMPP OBJECT FOR FST ---------------------------------------------

allele_freqs_fst <- dosage_fst / 4

stampp_df_fst <- as.data.frame(allele_freqs_fst)
stampp_df_fst$Sample <- rownames(stampp_df_fst)

stampp_final_fst <- stampp_df_fst |>
  left_join(meta_fst |> select(sra_run_accession, LI_population_MSL),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = LI_population_MSL) |>
  mutate(Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs_fst)))

stampp_obj_fst <- stamppConvert(stampp_final_fst, type = "r")
print("StAMPP object created successfully!")




## 6. FILTER BY MINIMUM SAMPLE SIZE --------------------------------------------

# Groups with fewer than 5 samples (= 20 tetraploid alleles) are dropped.
min_samples <- 5

meta_fst_robust <- meta_fst |>
  group_by(LI_population_MSL) |>
  filter(n() >= min_samples) |>
  ungroup()

print(paste("Dropping groups with fewer than", min_samples, "samples..."))
print("Populations remaining for robust FST:")
table(meta_fst_robust$LI_population_MSL)




## 7. RECALCULATE DOSAGE ON ROBUST SAMPLES --------------------------------------

samples_to_keep <- meta_fst_robust$sra_run_accession
dosage_fst <- dosage_filtered[rownames(dosage_filtered) %in% samples_to_keep, ]

col_vars_fst <- apply(dosage_fst, 2, var, na.rm = TRUE)
keep_cols_fst <- !is.na(col_vars_fst) & col_vars_fst > 0
dosage_fst <- dosage_fst[, keep_cols_fst]

print(paste("Samples remaining for FST:", nrow(dosage_fst)))
print(paste("SNPs remaining for FST:", ncol(dosage_fst)))




## 8. REBUILD STAMPP OBJECT FOR ROBUST FST --------------------------------------

allele_freqs_fst <- dosage_fst / 4

stampp_df_fst <- as.data.frame(allele_freqs_fst)
stampp_df_fst$Sample <- rownames(stampp_df_fst)

stampp_final_fst <- stampp_df_fst |>
  left_join(meta_fst_robust |> select(sra_run_accession, LI_population_MSL),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = LI_population_MSL) |>
  mutate(Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs_fst)))

stampp_obj_fst <- stamppConvert(stampp_final_fst, type = "r")
print("Robust StAMPP object created successfully!")




## 9. CALCULATE FST -------------------------------------------------------------

print("Starting FST calculation with bootstrapping...")

fst_results <- stamppFst(stampp_obj_fst, nboots = 100, percent = 95,
                         nclusters = parallel::detectCores() - 1)

fst_matrix <- fst_results$Fsts
pval_matrix <- fst_results$Pvalues

print("Pairwise FST matrix:")
print(fst_matrix)

save(fst_results, fst_matrix, pval_matrix,
     file = "Data/Fst/tetraploid_fst_LI_results_robust.RData")
print("FST results saved!")




## 10. PLOT FST HEATMAP ---------------------------------------------------------

fst_df <- as.data.frame(fst_matrix)
fst_df$Pop1 <- rownames(fst_df)

fst_long <- fst_df |>
  pivot_longer(cols = -Pop1, names_to = "Pop2", values_to = "FST") |>
  filter(!is.na(FST))

fst_heatmap <- ggplot(fst_long, aes(x = Pop2, y = Pop1, fill = FST)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027",
                       midpoint = median(fst_long$FST, na.rm = TRUE),
                       name = "FST") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  ggtitle("Pairwise FST: L. infantum-only Populations (MSL 0 vs 4)")

print(fst_heatmap)

ggsave("Figures/Fst/fst_heatmap_LI_tetraploid.jpeg",
       plot = fst_heatmap, width = 8, height = 8, dpi = 300)




## 11. PLOT FST BOXPLOT BY COMPARISON TYPE --------------------------------------

fst_boxplot_df <- fst_long |>
  mutate(
    msl1 = str_sub(Pop1, -1),
    msl2 = str_sub(Pop2, -1),
    comparison_type = case_when(
      msl1 == "0" & msl2 == "0" ~ "MSL0-0",
      msl1 == "4" & msl2 == "4" ~ "MSL4-4",
      TRUE ~ "MSL0-4"
    )
  ) |>
  mutate(comparison_type = factor(comparison_type, levels = c("MSL0-0", "MSL0-4", "MSL4-4")))

fst_boxplot <- ggplot(fst_boxplot_df, aes(x = comparison_type, y = FST, fill = comparison_type)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6, color = "black", size = 2) +
  scale_fill_manual(values = c("MSL0-0" = "#4575b4", "MSL0-4" = "#d73027", "MSL4-4" = "#fdae61")) +
  theme_minimal() +
  labs(
    title = "FST by MSL Copy Number",
    subtitle = "Chromosome 31 (Tetraploid, L. infantum-only)",
    x = "Comparison Category",
    y = "Pairwise FST"
  ) +
  theme(
    legend.position = "none",
    axis.text = element_text(size = 12, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 13, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

print(fst_boxplot)

ggsave("Figures/Fst/fst_boxplot_LI_tetraploid.jpeg",
       plot = fst_boxplot, width = 6, height = 6, dpi = 300)




## 12. PLOT LI2_0 vs LI2_4 BAR CHART -------------------------------------------

# LI2 = V2, the large mixed MSL0/MSL4 Brazilian population — key population of interest.

get_fst <- function(popA, popB, mat) {
  val1 <- mat[popA, popB]
  val2 <- mat[popB, popA]
  if (!is.na(val1)) return(val1)
  if (!is.na(val2)) return(val2)
  return(NA)
}

target_pops <- c("LI2_0", "LI2_4")
other_pops <- setdiff(rownames(fst_matrix), target_pops)

li2_df <- data.frame(
  Other_Pop = rep(other_pops, 2),
  LI2_Status = rep(target_pops, each = length(other_pops)),
  FST = NA
)

for (i in 1:nrow(li2_df)) {
  li2_df$FST[i] <- get_fst(li2_df$LI2_Status[i], li2_df$Other_Pop[i], fst_matrix)
}

li2_df <- li2_df |>
  mutate(
    Other_MSL = str_sub(Other_Pop, -1),
    Other_Pop = fct_reorder(Other_Pop, as.numeric(Other_MSL))
  )

li2_barplot <- ggplot(li2_df, aes(x = Other_Pop, y = FST, fill = LI2_Status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c("LI2_0" = "#4575b4", "LI2_4" = "#d73027"),
    name = "Comparing against:"
  ) +
  theme_minimal() +
  labs(
    title = "Divergence of LI2 Sub-populations\nfrom the Rest of L. infantum",
    subtitle = "Comparing FST of LI2_0 vs LI2_4 against all other groups",
    x = "Reference Population",
    y = "Pairwise FST"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
    axis.text.y = element_text(face = "bold", size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 11),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(size = 13, face = "bold"),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  )

print(li2_barplot)

ggsave("Figures/Fst/fst_LI2_comparison_barplot.jpeg",
       plot = li2_barplot, width = 10, height = 5, dpi = 300)




## 13. COMBINE INTO SINGLE FIGURE -----------------------------------------------

combined_figure <- fst_heatmap / (fst_boxplot | li2_barplot) +
  plot_annotation(
    tag_levels = "A",
    title = "Chromosome 31 Population Differentiation — L. infantum-only (Tetraploid)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(combined_figure)

ggsave("Figures/Fst/fst_combined_figure_LI_chr31.jpeg",
       plot = combined_figure, width = 18, height = 16, dpi = 300)
print("LI combined Fst figure saved!")



