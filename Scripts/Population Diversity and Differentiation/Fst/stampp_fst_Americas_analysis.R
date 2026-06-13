# SUMMARY -----------------------------------------------------------------------
# FST analysis — Americas only (AM subset, K=13)
# Dan's original script, adapted for local file paths.
# Uses tetraploid chr31 calls + diploid whole-genome calls
# Populations: AM1-AM13
# Groups by AM population x MSL status (0 or 4)
#
# Input RData files (created by load_and_transform_vcf_V2.R):
#   Data/Fst/cleaned_tetraploid_vcf_data_2026-04-13.RData
#   Data/Fst/cleaned_diploid_vcf_data_2026-04-13.RData




## 1. CLEAN UP AND LOAD LIBRARIES -----------------------------------------------

rm(list = ls())
setwd("/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN")
library(vcfR)
library(adegenet)
library(tidyverse)
library(StAMPP)
library(patchwork)




## 2. LOAD DATA -----------------------------------------------------------------

# Chr31 tetraploid data
load("fst_results/cleaned_tetraploid_vcf_data_2026-04-13.RData")
# Objects loaded: dosage_filtered, genlight_obj, meta_combined, stampp_obj

# Whole-genome diploid data
load("fst_results/cleaned_diploid_vcf_data_2026-04-13.RData")
# Objects loaded: dosage_filtered_dip, genlight_obj_dip, meta_combined_dip, stampp_obj_dip




## 3. PREPARE METADATA FOR FST --------------------------------------------------

# Combine population and MSL status into a single grouping variable (e.g., "AM1_0")
meta_combined <- meta_combined |>
  mutate(
    MSL_status = ifelse(is.na(MSL_coverage_rounded), "Unknown", as.character(MSL_coverage_rounded)),
    Americas_population_MSL = paste(Americas_population, MSL_status, sep = "_")
  )

table(meta_combined$Americas_population_MSL)
unique(meta_combined$Americas_population_MSL)




## 4. FILTER DATA FOR FST ANALYSIS ----------------------------------------------

# Keep only Americas populations and MSL 0 or 4
meta_fst <- meta_combined |>
  filter(Americas_population != "Non_Americas") |>
  filter(MSL_status %in% c("0", "4"))

print("Populations remaining for FST:")
table(meta_fst$Americas_population_MSL)

samples_to_keep <- meta_fst$sra_run_accession
dosage_fst <- dosage_filtered[rownames(dosage_filtered) %in% samples_to_keep, ]

# Recalculate variance to drop SNPs now invariant after subsetting
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
  left_join(meta_fst |> select(sra_run_accession, Americas_population_MSL),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = Americas_population_MSL) |>
  mutate(Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs_fst)))

stampp_obj_fst <- stamppConvert(stampp_final_fst, type = "r")
print("StAMPP object created successfully!")




## 6. FILTER BY MINIMUM SAMPLE SIZE --------------------------------------------

# Groups with fewer than 5 samples (= 20 tetraploid alleles) are dropped.
min_samples <- 5

meta_fst_robust <- meta_fst |>
  group_by(Americas_population_MSL) |>
  filter(n() >= min_samples) |>
  ungroup()

print(paste("Dropping groups with fewer than", min_samples, "samples..."))
print("Populations remaining for robust FST:")
table(meta_fst_robust$Americas_population_MSL)

# AM1_0 AM10_0 AM11_0 AM12_0  AM2_0  AM2_4  AM3_0  AM4_4  AM5_0  AM6_0  AM7_0  AM8_4  AM9_0
#    23     31     48      5     23     24     16     14     11     37     25     48     20




## 7. RECALCULATE DOSAGE ON ROBUST SAMPLES --------------------------------------

samples_to_keep <- meta_fst_robust$sra_run_accession
dosage_fst <- dosage_filtered[rownames(dosage_filtered) %in% samples_to_keep, ]

col_vars_fst <- apply(dosage_fst, 2, var, na.rm = TRUE)
keep_cols_fst <- !is.na(col_vars_fst) & col_vars_fst > 0
dosage_fst <- dosage_fst[, keep_cols_fst]

print(paste("Samples remaining for FST:", nrow(dosage_fst)))
# [1] "Samples remaining for FST: 325"

print(paste("SNPs remaining for FST:", ncol(dosage_fst)))
# [1] "SNPs remaining for FST: 13056"




## 8. REBUILD STAMPP OBJECT FOR ROBUST FST --------------------------------------

allele_freqs_fst <- dosage_fst / 4

stampp_df_fst <- as.data.frame(allele_freqs_fst)
stampp_df_fst$Sample <- rownames(stampp_df_fst)

stampp_final_fst <- stampp_df_fst |>
  left_join(meta_fst_robust |> select(sra_run_accession, Americas_population_MSL),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = Americas_population_MSL) |>
  mutate(Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs_fst)))

stampp_obj_fst <- stamppConvert(stampp_final_fst, type = "r")
print("Robust StAMPP object created successfully!")




## 9. CALCULATE FST -------------------------------------------------------------

print("Starting FST calculation with bootstrapping... this may take a few minutes!")

# nboots = 100 generates p-values via bootstrapping
# nclusters uses multiple CPU cores to speed up calculation
fst_results <- stamppFst(stampp_obj_fst, nboots = 100, percent = 95,
                         nclusters = parallel::detectCores() - 1)

fst_matrix <- fst_results$Fsts
pval_matrix <- fst_results$Pvalues

print("Pairwise FST matrix:")
print(fst_matrix)

save(fst_results, fst_matrix, pval_matrix,
     file = "fst_results/tetraploid_fst_AM_results_robust.RData")
print("FST results saved!")




## 10. PLOT FST HEATMAP ---------------------------------------------------------

fst_df <- as.data.frame(fst_matrix)
fst_df$Pop1 <- rownames(fst_df)

fst_long <- fst_df |>
  pivot_longer(cols = -Pop1, names_to = "Pop2", values_to = "FST") |>
  filter(!is.na(FST))

# Lock in population order so axes match
pop_order <- rownames(fst_matrix)
fst_long$Pop1 <- factor(fst_long$Pop1, levels = pop_order)
fst_long$Pop2 <- factor(fst_long$Pop2, levels = pop_order)

fst_heatmap <- ggplot(fst_long, aes(x = Pop1, y = Pop2, fill = FST)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "#4575b4", mid = "#ffffbf", high = "#d73027",
                       midpoint = median(fst_long$FST, na.rm = TRUE),
                       name = "Pairwise FST") +
  coord_fixed() +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold"),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  ggtitle("Pairwise FST: Americas Populations (MSL 0 vs 4)")

print(fst_heatmap)

ggsave("fst_results/figures/fst_heatmap_AM_tetraploid.jpeg",
       plot = fst_heatmap, width = 8, height = 8, dpi = 300)




## 11. PLOT FST BOXPLOT BY COMPARISON TYPE --------------------------------------

fst_boxplot_df <- fst_long |>
  mutate(
    msl1 = str_sub(Pop1, -1),
    msl2 = str_sub(Pop2, -1),
    # Categorise each pairwise comparison by MSL status combination
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
    subtitle = "Chromosome 31 (Tetraploid, Americas)",
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

ggsave("fst_results/figures/fst_boxplot_AM_tetraploid.jpeg",
       plot = fst_boxplot, width = 6, height = 6, dpi = 300)




## 12. PLOT AM2_0 vs AM2_4 BAR CHART -------------------------------------------

# AM2 = V2, the large mixed MSL0/MSL4 Brazilian population — key population of interest.

get_fst <- function(popA, popB, mat) {
  val1 <- mat[popA, popB]
  val2 <- mat[popB, popA]
  if (!is.na(val1)) return(val1)
  if (!is.na(val2)) return(val2)
  return(NA)
}

target_pops <- c("AM2_0", "AM2_4")
other_pops <- setdiff(rownames(fst_matrix), target_pops)

am2_df <- data.frame(
  Other_Pop = rep(other_pops, 2),
  AM2_Status = rep(target_pops, each = length(other_pops)),
  FST = NA
)

for (i in 1:nrow(am2_df)) {
  am2_df$FST[i] <- get_fst(am2_df$AM2_Status[i], am2_df$Other_Pop[i], fst_matrix)
}

# Sort x-axis so MSL0 populations appear left, MSL4 right
am2_df <- am2_df |>
  mutate(
    Other_MSL = str_sub(Other_Pop, -1),
    Other_Pop = fct_reorder(Other_Pop, as.numeric(Other_MSL))
  )

am2_barplot <- ggplot(am2_df, aes(x = Other_Pop, y = FST, fill = AM2_Status)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(
    values = c("AM2_0" = "#4575b4", "AM2_4" = "#d73027"),
    name = "Comparing against:"
  ) +
  theme_minimal() +
  labs(
    title = "Divergence of AM2 Sub-populations\nfrom the Rest of the Americas",
    subtitle = "AM2_0 vs AM2_4 comparisons",
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

print(am2_barplot)

ggsave("fst_results/figures/fst_AM2_comparison_barplot.jpeg",
       plot = am2_barplot, width = 10, height = 5, dpi = 300)




## 13. COMBINE INTO SINGLE FIGURE -----------------------------------------------

combined_figure <- fst_heatmap / (fst_boxplot | am2_barplot) +
  plot_annotation(
    tag_levels = "A",
    title = "Chromosome 31 Population Differentiation — Americas (Tetraploid)",
    theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
  )

print(combined_figure)

ggsave("fst_results/figures/fst_combined_figure_AM_chr31.jpeg",
       plot = combined_figure, width = 18, height = 16, dpi = 300)
print("AM combined Fst figure saved!")