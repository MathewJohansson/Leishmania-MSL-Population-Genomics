# PCA PLOTS DIPLOID SNP CALLS (ALL NON CHR31 CHROMOSOMES) ---------------------

## CLEAN UP AND LOAD LIBRARIES ----------------------------------------------------------------
rm(list = ls())
library(vcfR)
library(adegenet)
library(tidyverse)
library(patchwork) 
library(gtools)     
library(Polychrome)

## LOAD DATA --------------------------------------------------------------------

# 1. Load sample info and immediately fix the population naming
samples <- read_tsv("sample_info/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
names(samples)[1] <- "sra_run_accession"

# Standardize population naming early: Replace 'V' with 'AM'
samples <- samples %>%
  mutate(Americas_population = str_replace(population, "^V", "AM"))

# 2. Load diploid VCF
vcf <- read.vcfR("vcf_files/Americas_only_MSL0and4_samples.noLinJ31.vcf.gz")

## PROCESS VCF FOR PCA DATA (MEMORY OPTIMIZED) ---------------------------------

# Extract info to filter Chr31
chrom <- getCHROM(vcf)
vcf_no_msl <- vcf[chrom != "LinJ.31", ]

# 1. Convert directly to genlight (This is where you save all your RAM)
# vcfR2genlight converts "0/0", "0/1", "1/1" into 0, 1, 2 efficiently
gl_data <- vcfR2genlight(vcf_no_msl)

# 2. Filter monomorphic SNPs (zero variance)
# PCA cannot run if a SNP is the same in every single sample
n_samples <- nInd(gl_data)
snps_to_keep <- !is.na(glSum(gl_data)) & 
  glSum(gl_data) > 0 & 
  glSum(gl_data) < (2 * n_samples)

gl_filtered <- gl_data[, snps_to_keep]

print(paste("Original SNPs:", nLoc(gl_data)))
#[1] "Original SNPs: 78930"
print(paste("SNPs after removing monomorphic sites:", nLoc(gl_filtered)))
#[1] "SNPs after removing monomorphic sites: 41436"

## RUN PCA ANALYSIS ------------------------------------------------------------


# 1. Extract the numeric matrix and handle missing data
pca_input <- tab(gl_filtered, NA.method = "mean")

# 2a. Identify and remove constant columns (variance = 0)
# This is the safety check for scale. = TRUE
col_vars <- apply(pca_input, 2, var)
pca_input_final <- pca_input[, col_vars > 0]

print(paste("Dropped", sum(col_vars == 0), "additional constant columns."))
print(paste("Final SNP count for PCA:", ncol(pca_input_final)))

# 2b. Run PCA using prcomp
pca_results <- prcomp(pca_input_final, center = TRUE, scale. = TRUE, rank. = 6)


# 3. Calculate Variance Explained
# For prcomp, eigenvalues are sdev^2
pc_var <- pca_results$sdev^2
pc_perc <- round(pc_var / sum(pc_var) * 100, 1)

# 4. Prepare Dataframe for Plotting
# Note: prcomp scores are stored in 'x', not 'scores'
pca_df <- as.data.frame(pca_results$x)
pca_df$sra_run_accession <- rownames(pca_df)

# Merge with corrected metadata
pca_combined <- pca_df %>%
  left_join(samples, by = "sra_run_accession") %>%
  mutate(Americas_population = factor(Americas_population, 
                                      levels = mixedsort(unique(Americas_population))))


## PLOTTING PREPARATIONS ------------------------------------------------------

paul_tol_palette <- c(
  "#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", 
  "#AA3377", "#BBBBBB", "#EE7733", "#009988", "#994F00", 
  "#332288", "#9955BB", "#FFAABB"
)
palette_AM <- setNames(paul_tol_palette[1:13], paste0("AM", 1:13))

# Dynamic axis limits (97% of data to ignore extreme outliers)
pc_cols <- grep("^PC", names(pca_combined), value = TRUE)
lims <- lapply(pca_combined[pc_cols], function(x) quantile(x, probs = c(0.015, 0.985), na.rm = TRUE))

## PLOTTING FUNCTION (To save space and avoid repetition) ---------------------

plot_pca <- function(df, x_pc, y_pc, fill_var, palette, title_suffix, pc_perc_vec) {
  x_idx <- as.numeric(gsub("PC", "", x_pc))
  y_idx <- as.numeric(gsub("PC", "", y_pc))
  
  ggplot(df, aes_string(x = x_pc, y = y_pc, fill = fill_var)) +
    geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
    scale_fill_manual(values = palette, name = gsub("_", " ", fill_var)) +
    theme_minimal() +
    labs(
      x = paste0(x_pc, " (", pc_perc_vec[x_idx], "%)"),
      y = paste0(y_pc, " (", pc_perc_vec[y_idx], "%)"),
      title = paste("Chr 1-30:", title_suffix)
    ) +
    coord_cartesian(xlim = lims[[x_pc]], ylim = lims[[y_pc]]) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
}

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


### PCA1.PCA2 ------------------------------------------------------------------

# Define the structure
PC1.PC2.base <- ggplot(pca_combined, aes(x = PC1, y = PC2)) +
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
PC1.PC2.base

PC1.PC2.population <- PC1.PC2.base + 
  aes(fill = Americas_population) + 
  scale_fill_manual(values = palette_AM, name = "Americas Population") +
  guides(fill = guide_legend(ncol = 5))+
  ggtitle("")
PC1.PC2.population

PC1.PC2.MSL <- PC1.PC2.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_manual(values = c("0" = "#d62728", "1" = "#ff7f0e", "2" = "#ffff00", "3" = "#2ca02c", "4" = "#1f77b4"), name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list
PC1.PC2.MSL

# arrange the two plots side by side
PC1.PC2.population.MSL <- (PC1.PC2.population | PC1.PC2.MSL) 
PC1.PC2.population.MSL

# Save the combined plot
ggsave("plots/diploid.nonchr31.PCA1.PCA2.population.MSL.standard_colours.jpeg",
       plot = PCA1.PCA2.population.MSL, width = 14, height = 7, dpi = 300)



### PCA1.PCA2 ------------------------------------------------------------------

# Define the structure
PC3.PC4.base <- ggplot(pca_combined, aes(x = PC3, y = PC4)) +
  # We leave 'fill' out of the global aes and put it here as a placeholder
  geom_point(color = "black", shape = 21, size = 3, stroke = 0.5, alpha = 0.8) +
  theme_minimal() +
  labs(
    x = paste0("PC3 (", pc_perc[3], "%)"),
    y = paste0("PC4 (", pc_perc[4], "%)")
  ) +
  coord_cartesian(xlim = lims$PC3*c(1,0.1), ylim = lims$PC4)+
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    title = element_text(size = 14, face = "bold")
  )
PC3.PC4.base


PC3.PC4.population <- PC3.PC4.base + 
  aes(fill = Americas_population) + 
  scale_fill_manual(values = palette_AM, name = "Americas Population") +
  guides(fill = guide_legend(ncol = 5))+
  ggtitle("")
PC3.PC4.population

PC3.PC4.MSL <-PC3.PC4.base + 
  aes(fill = as.factor(MSL_coverage_rounded)) + 
  scale_fill_manual(values = c("0" = "#d62728", "1" = "#ff7f0e", "2" = "#ffff00", "3" = "#2ca02c", "4" = "#1f77b4"), name = "MSL Copies") +
  guides(fill = guide_legend(ncol = 1)) # Keep this one as a simple list
PC3.PC4.MSL

# arrange the two plots side by side
PC3.PC4.population.MSL <- (PC3.PC4.population | PC3.PC4.MSL) 
PC3.PC4.population.MSL

# Save the combined plot
ggsave("plots/diploid.nonchr31.PCA3.PCA4.population.MSL.standard_colours.jpeg",
       plot = PC3.PC4.population.MSL, width = 14, height = 7, dpi = 300)


# 1. Arrange the four plots into a 2x2 grid
# Top row: PC1 vs PC2 (Pop and MSL)
# Bottom row: PC3 vs PC4 (Pop and MSL)
combined_plot <- (PCA1.PCA2.population + PCA1.PCA2.MSL) / 
  (PC3.PC4.population + PC3.PC4.MSL)

# 2. Collect the legends and format the layout
combined_plot <- combined_plot + 
  plot_layout(guides = "collect") + 
  plot_annotation(tag_levels = "A") & 
  theme(legend.position = "bottom")

# 3. View the result
combined_plot

# 4. Save the high-res 4-panel figure
ggsave("plots/diploid.nonchr31.FullPCA.population.MSL.standard_colours.jpeg",
       plot = combined_plot, width = 14, height = 12, dpi = 300)


#save all the data
save.image("rdata/diploid_PCA.standard_colours.RData")


# END OF SCRIPT ----------------------------------------------------------------



