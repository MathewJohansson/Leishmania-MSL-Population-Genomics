# CALCULATING NEI GENETIC DISTANCE TETRAPLOID SNP CALLS (on CHROMOSOME 31) -----

# UPDATE: 2026-04-10 using filtered tetraploid VCF
# tetraploid.snps.filter.Americas_only_MSL0and4_samples.vcf.gz




## CLEAN UP AND LOAD LIBRARIES -------------------------------------------------

rm(list = ls())
library(vcfR)
library(adegenet)
library(ape)
library(tidyverse)
library(StAMPP)
library(ggtree)
library(viridis)




## LOAD DATA -------------------------------------------------------------------

#sample info
samples <- read_tsv("SAMPLE_INFO/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
names(samples)[1]="sra_run_accession"

# Load tetraploid VCF
# IMPORTANT: Use your original FreeBayes VCF here, NOT the pseudo-diploid one!
vcf <- read.vcfR("PROCESSED_DATA/tetraploid.snps.filter.Americas_only_MSL0and4_samples.vcf.gz")




## PROCESS VCF FOR DOSAGE DATA -------------------------------------------------

# Filter out the MSL region to prevent structural/selection bias
chrom <- getCHROM(vcf)
pos <- as.numeric(getPOS(vcf))

# Create a logical filter: keep SNPs that are NOT in the MSL coordinate window
keep_snps <- !(chrom == "LinJ.31" & pos >= 1181282 & pos <= 1192406)

# Subset the VCF using this filter
vcf_no_msl <- vcf[keep_snps, ]

# 1. Extract raw genotype strings (e.g., "0/0/1/1")
gt_strings <- extract.gt(vcf_no_msl, element = "GT")

# 2. Function to count the '1' alleles in each string
count_alt_alleles <- function(x) {
  lengths(regmatches(x, gregexpr("1", x)))
}

# 3. Apply the count across the matrix
dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings),
  ncol = ncol(gt_strings)
)

# Restore names
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)

# Set actual missing data back to NA 
dosage_matrix[is.na(gt_strings) | gt_strings == "./." | gt_strings == "./././."] <- NA

# 4. Transpose (Samples as rows)
dosage_matrix <- t(dosage_matrix)

# 5a. Impute missing values with the mean dosage of each SNP
dosage_imputed <- apply(dosage_matrix, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

# 5b. DROP OUTLIERS HERE (Before variance filtering!)
outliers_to_remove <- c("ERR205787", "SRR8608748")
dosage_imputed <- dosage_imputed[!rownames(dosage_imputed) %in% outliers_to_remove, ]
print(paste("Removed outliers. New sample count:", nrow(dosage_imputed)))

# 5c. Calculate variance, ignoring NAs, to drop invariant sites
column_vars <- apply(dosage_imputed, 2, var, na.rm = TRUE)

# Create a clean mask: Keep if variance is NOT NA AND variance is GREATER than 0
keep_cols <- !is.na(column_vars) & column_vars > 0
dosage_filtered <- dosage_imputed[, keep_cols]

# Check the final counts
print(paste("Total SNPs initially:", ncol(dosage_imputed)))
print(paste("Remaining SNPs for StAMPP:", ncol(dosage_filtered)))




## PREPARE METADATA (Replacing the old PCA join) -------------------------------

# Create a simple dataframe of the samples we actually kept
meta_combined <- data.frame(sra_run_accession = rownames(dosage_filtered)) |>
  left_join(samples, by = "sra_run_accession")




## CALCULATE NEI'S GENETIC DISTANCE (StAMPP) -----------------------------------

# 1. Convert dosages (0-4) to allele frequencies (0.0-1.0)
allele_freqs <- dosage_filtered / 4

# 2. Build the StAMPP data frame
stampp_df <- as.data.frame(allele_freqs)
stampp_df$Sample <- rownames(stampp_df)

# Safe join using meta_combined
stampp_df <- stampp_df |>
  left_join(meta_combined |> select(sra_run_accession, population), 
            by = c("Sample" = "sra_run_accession"))

# Format the final columns
stampp_final <- stampp_df |>
  rename(Pop = population) |>  
  mutate(
    Ploidy = 4,
    Format = "freq"
  ) |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs)))

# 3. Convert to a StAMPP internal object
stampp_obj <- stamppConvert(stampp_final, type = "r")

# 4. Calculate Nei's Genetic Distance Matrix
nei_dist_matrix <- stamppNeisD(stampp_obj, pop = FALSE)

# 5. Export for Mathew (SplitsTree Network)
stamppPhylip(nei_dist_matrix, file = "PROCESSED_DATA/tetraploid_chr31_Neis_distance.txt")

# 6. Export a Newick Tree for ggtree / iTOL
nj_tree <- nj(as.dist(nei_dist_matrix))
write.tree(nj_tree, file = "PROCESSED_DATA/tetraploid_chr31_Neis_tree.nwk")
print("StAMPP Nei's Genetic Distance calculated and exported successfully!")




## PLOT PHYLOGENY WITH MSL ANNOTATION ------------------------------------------

# 1. Prepare the metadata for ggtree
tree_meta <- meta_combined |>
  select(sra_run_accession, population, MSL_coverage_rounded) |>
  distinct() |>
  mutate(
    # Handle NAs to avoid coloring errors
    MSL_coverage = ifelse(is.na(MSL_coverage_rounded), "Unknown", as.character(MSL_coverage_rounded)),
    MSL_coverage = as.factor(MSL_coverage)
  )

# 2. Build the tree plot
nj_plot.msl <- ggtree(nj_tree, layout = "daylight") %<+% tree_meta +
  # Add points at the tips, colored by MSL status
  geom_tippoint(aes(color = MSL_coverage), size = 2.5, alpha = 0.8) +
  # Apply your trusted viridis color scale
  scale_color_viridis_d(name = "MSL Copies") +
  theme(
    legend.position = "right",
    title = element_text(size = 14, face = "bold")
  ) 

# 3. View the plot in your RStudio viewer
print(nj_plot.msl)

# 4. Save the plot
ggsave("PLOTS/tetraploid_chr31_NJ_tree.msl.jpeg", 
       plot = nj_plot.msl, width = 10, height = 10, dpi = 300)

print("Phylogeny plotted and saved successfully!")




## PLOT PHYLOGENY BY POPULATION ------------------------------------------------

# Make sure these libraries are loaded (you can also move these to the top of your script)
library(Polychrome)
library(gtools)

# 1. Ensure population is a factor and naturally sorted (V1, V2... V10, V11, V12)
tree_meta$population <- factor(
  tree_meta$population, 
  levels = mixedsort(unique(tree_meta$population))
)

# 2. Recreate your trusty 13-color palette
all_cols <- glasbey.colors(15)
my_cols_no_grey <- as.vector(all_cols[c(2:8, 10:15)])

# 3. Build the population tree plot
nj_plot.populations <- suppressWarnings(
  ggtree(nj_tree, layout = "daylight") %<+% tree_meta +
    # Color by population this time
    geom_tippoint(aes(color = population), size = 2.5, alpha = 0.8) +
    # Apply your custom manual palette
    scale_color_manual(values = my_cols_no_grey, name = "Population") +
    theme(
      legend.position = "right",
      title = element_text(size = 14, face = "bold")
    ) +
    # Optional: Makes the legend dots a bit larger to see the colors clearly
    guides(color = guide_legend(override.aes = list(size = 4))) 
)

# 4. View the plot in your RStudio viewer
print(nj_plot.populations)

#ggsave this plot
ggsave("PLOTS/tetraploid_chr31_NJ_tree_population.jpeg", 
       plot = nj_plot.populations, width = 10, height = 10, dpi = 300)


#ggarrange nj_plot.msl and nj_plot.populations and save as a single figure
library(ggpubr)
nj_plot.msl.populations <- ggarrange(nj_plot.msl, nj_plot.populations, 
  ncol = 2, nrow = 1, common.legend = F, legend = "bottom", labels = c("A", "B"), 
  font.label = list(size = 16, face = "bold"))
nj_plot.msl.populations

# 5. Save the new plot
ggsave("PLOTS/tetraploid_chr31_NJ_tree.msl.populations.jpeg", 
       plot = nj_plot.msl.populations, width = 20, height = 10, dpi = 300)




## SAVE WORKSPACE --------------------------------------------------------------
# Save all the distance and tree data
save.image("RDATA/tetraploid_Neis.RData")



