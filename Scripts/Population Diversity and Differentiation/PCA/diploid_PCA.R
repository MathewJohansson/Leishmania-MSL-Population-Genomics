# PCA PLOTS DIPLLOID SNP CALLS (ALL NON CHR31 CHROMOSOMES) ---------------------

# using filtered diploid VCF, for strains from Americas, with MSL0 and MSL4 only
# Linfantum_samples_no_outgroup_biallelic.2025-06-10.noLinJ.31.vcf.gz

## CLEAN UP AND LOAD LIBRARIES ----------------------------------------------------------------
rm(list = ls())
library(vcfR)
library(adegenet)
library(ape)
library(tidyverse)
library(patchwork)  #for plotting panels
library(gtools)     #for mixedsort() to order population names in a natural way
library(Polychrome) #for visual accessibility, many colours in plots

## LOAD DATA --------------------------------------------------------------------

#sample info
samples <- read_tsv("SAMPLE_INFO/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
names(samples)[1]="sra_run_accession"

names(samples)
glimpse(samples)

# Load diploid VCF
vcf <- read.vcfR("PROCESSED_DATA/Americas_only_MSL0and4_samples.noLinJ31.vcf.gz")

## PROCESS VCF FOR PCA DATA ----------------------------------------------------

# Filter regions 
# Extract the chromosome and position arrays from the VCF
chrom <- getCHROM(vcf)
pos <- as.numeric(getPOS(vcf))

# Create a logical filter: keep SNPs that are NOT in the MSL coordinate window
# (Using 1181282 as the start of the first gene and 1192406 as the end of the last)
keep_snps <- !(chrom == "LinJ.31")

# Subset the VCF using this filter
vcf_no_msl <- vcf[keep_snps, ]
dim(vcf_no_msl)
#variants fix_cols  gt_cols 
#86265        8      464 

# 1. Extract raw genotype strings (e.g., "0/1")
gt_strings <- extract.gt(vcf_no_msl, element = "GT")

# 2. Function to count the '1' alleles in each string
# This counts every "1" followed by a separator or end-of-string
count_alt_alleles <- function(x) {
  # Count occurrences of "1"
  # stringr is faster if you have it, but base R works:
  lengths(regmatches(x, gregexpr("1", x)))
}

# 3. Apply the count across the matrix
# We handle NAs by ensuring they stay NA
dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings),
  ncol = ncol(gt_strings)
)

#save all the data
save.image("RDATA/diploid_PCA.RData")
load("RDATA/diploid_PCA.RData")



#dosage_matrix[15:20,1:10]
max(dosage_matrix)

# Restore names
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)

# Set actual missing data back to NA (the regex might count 0 for "./.")
dosage_matrix[is.na(gt_strings) | gt_strings == "./." | gt_strings == "./."] <- NA

# 4. Transpose for PCA (Samples as rows)
dosage_matrix <- t(dosage_matrix)

# --- Check the results ---
print("Dosage distribution:")
table(dosage_matrix, useNA = "always") 
# You should now see counts for 0, 1, 2


# 5a. Impute missing values with the mean dosage of each SNP
# We calculate means column-wise (na.rm = TRUE)
snp_means <- colMeans(dosage_matrix, na.rm = TRUE)

# Replace NAs in the matrix
dosage_imputed <- apply(dosage_matrix, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
})

# 5b.1. Calculate variance, ignoring NAs
column_vars <- apply(dosage_imputed, 2, var, na.rm = TRUE)

# 5b.2. Create a clean mask: 
# Keep if variance is NOT NA AND variance is GREATER than 0
keep_cols <- !is.na(column_vars) & column_vars > 0
dosage_filtered <- dosage_imputed[, keep_cols]

# Check the results again
print(paste("Total SNPs initially:", ncol(dosage_imputed)))
print(paste("Removed (Invariant or all NA):", sum(!keep_cols)))
print(paste("Remaining SNPs for PCA:", ncol(dosage_filtered)))

## RUN PCA ANALYSIS AND CREATE DF FOR PLOTTING  --------------------------------

# Calculate only the top 6 PCs
pca_results <- prcomp(dosage_filtered, center = TRUE, scale. = TRUE, rank. = 6)
save(pca_results, file = "RDATA/diploid_PCA_results.RData")
load("RDATA/diploid_PCA_results.RData")

# Calculate the % variance explained for each PC
# sdev^2 gives the eigenvalues (variance)
pc_var <- pca_results$sdev^2
pc_perc <- round(pc_var / sum(pc_var) * 100, 1)

# Now your plot labels will have the correct numbers!

# 1. Prepare PCA data
pca_df <- as.data.frame(pca_results$x)
pca_df$sra_run_accession <- rownames(pca_df)

# 2. Merge with metadata
# We use left_join so we keep all 466 VCF samples
names(pca_df)
pca_combined <- pca_df %>%
  left_join(samples, by = "sra_run_accession")

# 3. Quick check for that one missing sample
missing_sample <- pca_combined %>% 
  filter(is.na(MSL_coverage_rounded)) %>% 
  pull(sra_run_accession)

print(paste("Sample with no metadata:", missing_sample))
names(pca_combined)


## PLOTTING PREPARATIONS  ------------------------------------------------------

# We treat coverage_rounded as a factor to get distinct colors

# Ensure you have your pc_perc defined
pc_perc <- round(pca_results$sdev^2 / sum(pca_results$sdev^2) * 100, 1)


# Convert population to a factor with 'natural' alphanumeric order
pca_combined$population <- factor(
  pca_combined$population, 
  levels = mixedsort(unique(pca_combined$population))
)

# 1. Prepare your 13 non-grey colors (as we did before)
all_cols <- glasbey.colors(15)
my_cols_no_grey <- as.vector(all_cols[c(2:8, 10:15)])


## PLOTS -----------------------------------------------------------------------

#Plot limits for clarity
# Create a named list of limits for all PCs automatically
# lims$PC1 gives you the min/max for PC1
# lims$PC4 gives you the min/max for PC4
#etc
pc_cols <- grep("^PC", names(pca_combined), value = TRUE) # Finds all columns starting with "PC"

lims <- lapply(pca_combined[pc_cols], function(x) {
  quantile(x, probs = c(0.015, 0.985), na.rm = TRUE)
})

### PCA1.PCA4 ------------------------------------------------------------------

# Define the structure
PC1.PC4.base <- ggplot(pca_combined, aes(x = PC1, y = PC4)) +
  # We leave 'fill' out of the global aes and put it here as a placeholder
  geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", pc_perc[1], "%)"),
    y = paste0("PC4 (", pc_perc[4], "%)")
  ) +
  coord_cartesian(xlim = lims$PC1*c(1,8), ylim = lims$PC4) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    title = element_text(size = 14, face = "bold")
  )
PC1.PC4.base

PCA1.PCA4.population <- PC1.PC4.base + 
  aes(fill = population) + 
  scale_fill_manual(values = my_cols_no_grey, name = "Population") +
  guides(fill = guide_legend(ncol = 5)) +
  ggtitle("all other chromosomes")
PCA1.PCA4.population

PCA1.PCA4.MSL <- PC1.PC4.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_viridis_d(name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list
PCA1.PCA4.MSL

# arrange the two plots side by side
PCA1.PCA4.population.MSL <- (PCA1.PCA4.population | PCA1.PCA4.MSL) 
PCA1.PCA4.population.MSL

# Save the combined plot
ggsave("PLOTS/diploid.nonchr31.PCA1.PCA4.population.MSL.jpeg",
       plot = PCA1.PCA4.population.MSL, width = 14, height = 7, dpi = 300)

### PCA1.PCA1 ------------------------------------------------------------------

# Define the structure
PC1.PC1.base <- ggplot(pca_combined, aes(x = PC1, y = PC2)) +
  # We leave 'fill' out of the global aes and put it here as a placeholder
  geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC1 (", pc_perc[1], "%)"),
    y = paste0("PC2 (", pc_perc[2], "%)")
  ) +
  coord_cartesian(xlim = lims$PC1*c(1,10), ylim = lims$PC2)+
theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    title = element_text(size = 14, face = "bold")
  )
PC1.PC1.base

PCA1.PCA2.population <- PC1.PC1.base + 
  aes(fill = population) + 
  scale_fill_manual(values = my_cols_no_grey, name = "Population") +
  guides(fill = guide_legend(ncol = 5))+
  ggtitle("all other chromosomes")
PCA1.PCA2.population

PCA1.PCA2.MSL <-PC1.PC1.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_viridis_d(name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list
PCA1.PCA2.MSL

# arrange the two plots side by side
PCA1.PCA2.population.MSL <- (PCA1.PCA2.population | PCA1.PCA2.MSL) 
PCA1.PCA2.population.MSL

# Save the combined plot
ggsave("PLOTS/diploid.nonchr31.PCA1.PCA2.population.MSL.jpeg",
       plot = PCA1.PCA2.population.MSL, width = 14, height = 7, dpi = 300)

### PCA2.PCA3 ------------------------------------------------------------------

# Define the structure
PC2.PC3.base <- ggplot(pca_combined, aes(x = PC2, y = PC3)) +
  # We leave 'fill' out of the global aes and put it here as a placeholder
  geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC2 (", pc_perc[2], "%)"),
    y = paste0("PC3 (", pc_perc[3], "%)")
  ) +
  coord_cartesian(xlim = lims$PC2*c(0.1,1), ylim = lims$PC3*c(1,0.1))+
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    title = element_text(size = 14, face = "bold")
  )
PC2.PC3.base

PCA2.PCA3.population <- PC2.PC3.base + 
  aes(fill = population) + 
  scale_fill_manual(values = my_cols_no_grey, name = "Population") +
  guides(fill = guide_legend(ncol = 5))+
  ggtitle("all other chromosomes")
PCA2.PCA3.population

PCA2.PCA3.MSL <-PC2.PC3.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_viridis_d(name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list
PCA2.PCA3.MSL

# arrange the two plots side by side
PCA2.PCA3.population.MSL <- (PCA2.PCA3.population | PCA2.PCA3.MSL) 
PCA2.PCA3.population.MSL

# Save the combined plot
ggsave("PLOTS/diploid.nonchr31.PCA2.PCA3.population.MSL.jpeg",
       plot = PCA2.PCA3.population.MSL, width = 14, height = 7, dpi = 300)



### PCA5.PCA6 ------------------------------------------------------------------

# Define the structure
PC5.PC6.base <- ggplot(pca_combined, aes(x = PC5, y = PC6)) +
  # We leave 'fill' out of the global aes and put it here as a placeholder
  geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC5 (", pc_perc[5], "%)"),
    y = paste0("PC6 (", pc_perc[6], "%)")
  ) +
  coord_cartesian(xlim = lims$PC5, ylim = lims$PC6)+
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    title = element_text(size = 14, face = "bold")
  )
PC5.PC6.base

PC5.PC6.population <- PC5.PC6.base + 
  aes(fill = population) + 
  scale_fill_manual(values = my_cols_no_grey, name = "Population") +
  guides(fill = guide_legend(ncol = 5))+
  ggtitle("all other chromosomes")
PC5.PC6.population

PC5.PC6.MSL <-PC5.PC6.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_viridis_d(name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list
PC5.PC6.MSL

# arrange the two plots side by side
PC5.PC6.population.MSL <- (PC5.PC6.population | PC5.PC6.MSL) 
PC5.PC6.population.MSL


# Save the combined plot
ggsave("PLOTS/diploid.nonchr31.PCA5.PCA6.population.MSL.jpeg",
       plot = PC5.PC6.population.MSL, width = 14, height = 7, dpi = 300)



### PCA2.PCA4 ------------------------------------------------------------------

# Define the structure
PC2.PC4.base <- ggplot(pca_combined, aes(x = PC2, y = PC4)) +
  # We leave 'fill' out of the global aes and put it here as a placeholder
  geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC2 (", pc_perc[2], "%)"),
    y = paste0("PC4 (", pc_perc[4], "%)")
  ) +
  coord_cartesian(xlim = lims$PC2*c(0.1,1), ylim = lims$PC4) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    title = element_text(size = 14, face = "bold")
  )
PC2.PC4.base

PCA2.PCA4.population <- PC2.PC4.base + 
  aes(fill = population) + 
  scale_fill_manual(values = my_cols_no_grey, name = "Population") +
  guides(fill = guide_legend(ncol = 5)) +
  ggtitle("all other chromosomes")
PCA2.PCA4.population

PCA2.PCA4.MSL <- PC2.PC4.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_viridis_d(name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list

# arrange the two plots side by side
PCA2.PCA4.population.MSL <- (PCA2.PCA4.population | PCA2.PCA4.MSL) 
PCA2.PCA4.population.MSL

# Save the combined plot
ggsave("PLOTS/diploid.nonchr31.PCA2.PCA4.population.MSL.jpeg",
       plot = PCA2.PCA4.population.MSL, width = 14, height = 7, dpi = 300)


#save all the data
save.image("RDATA/tetraploid_PCA.RData")
load("RDATA/tetraploid_PCA.RData")

# END OF SCRIPT ----------------------------------------------------------------



