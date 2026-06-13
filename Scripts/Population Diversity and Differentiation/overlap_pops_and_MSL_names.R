#!/usr/bin/env Rscript

# =============================================================================
# Sample Overlap Analysis: MSL Status by Population
# =============================================================================
# 
# PURPOSE:
#   This script finds overlaps between MSL copy number status samples (MSL0/MSL4) 
#   and ADMIXTURE population assignments, creating population-specific MSL sample lists.
#
# WORKFLOW:
#   1. Load MSL status sample lists (MSL0 = zero copies, MSL4 = four copies)
#   2. Find all population-specific sample files from ADMIXTURE analysis
#   3. For each population, find overlap with MSL0 and MSL4 samples
#   4. Create new sample lists for each population × MSL status combination
#   5. Generate summary statistics and visualizations
#
# INPUT FILES:
#   - MSL status files: sample_info/Linfantum_only.MSL{0,4}.pseudo-diploid.list
#   - Population files: Americas.samples.diploidised_names.txt.populations.K13.vcfspopgen/*.pop*.txt
#
# OUTPUT FILES:
#   - Per-population MSL sample lists: Americas.K13.pop{N}.MSL{0,4}.diploidised.sample_names.txt
#   - Summary table: MSL_population_summary.tsv
#
# AUTHOR: GitHub Copilot
# DATE: July 15, 2025
# VERSION: 1.0
# =============================================================================

library(tidyverse)  # For data manipulation and file I/O

# =============================================================================
# INPUT FILE PATHS AND CONFIGURATION
# =============================================================================

# Define paths to MSL status sample lists
# These files contain sample names (one per line) for samples with specific MSL copy numbers
msl0_file <- "sample_info/Linfantum_only.MSL0.pseudo-diploid.list"  # Zero MSL copies
msl4_file <- "sample_info/Linfantum_only.MSL4.pseudo-diploid.list"  # Four MSL copies

# Define directory containing population-specific sample files from ADMIXTURE analysis
# These files were created by admixture_to_vcfsubset_to_popstats.pl
pop_dir <- "Americas.samples.diploidised_names.txt.populations.K13.vcfspopgen"

# =============================================================================
# LOAD MSL STATUS SAMPLES
# =============================================================================

cat("Loading MSL status sample lists...\n")

# Read MSL0 samples (samples with zero MSL copies)
# These represent samples that have lost the MSL locus entirely
if (file.exists(msl0_file)) {
    # Read all lines from file, trim whitespace, and remove empty lines
    msl0_samples <- read_lines(msl0_file) %>%
        str_trim() %>%                    # Remove leading/trailing whitespace
        .[. != ""]                       # Filter out empty lines using logical indexing
    cat("Loaded", length(msl0_samples), "MSL0 samples\n")
} else {
    stop("MSL0 file not found: ", msl0_file)
}

# Read MSL4 samples (samples with four MSL copies)
# These represent samples with the full complement of MSL loci
if (file.exists(msl4_file)) {
    # Same processing as MSL0 samples
    msl4_samples <- read_lines(msl4_file) %>%
        str_trim() %>%
        .[. != ""]
    cat("Loaded", length(msl4_samples), "MSL4 samples\n")
} else {
    stop("MSL4 file not found: ", msl4_file)
}

# =============================================================================
# FIND AND PROCESS POPULATION FILES
# =============================================================================

cat("Finding population files...\n")

# Create file pattern to match all population files
# Pattern matches: Americas.samples.diploidised_names.populations.K13.pop*.txt
# The * wildcard will match any population number (pop1, pop2, etc.)
pop_pattern <- file.path(pop_dir, "Americas.samples.diploidised_names.populations.K13.pop*.txt")

# Use Sys.glob() to find all files matching the pattern
# This returns a character vector of full file paths
pop_files <- Sys.glob(pop_pattern)

# Check if any population files were found
if (length(pop_files) == 0) {
    stop("No population files found matching pattern: ", pop_pattern)
}

cat("Found", length(pop_files), "population files\n")

# Sort files by population number for consistent processing order
# Extract population numbers using regex, convert to numeric, and sort
# str_extract() finds the first match of the pattern "pop(\\d+)" and extracts group 1 (the digits)
pop_files <- pop_files[order(as.numeric(str_extract(pop_files, "pop(\\d+)", group = 1)))]

# =============================================================================
# PROCESS EACH POPULATION FILE AND FIND OVERLAPS
# =============================================================================

cat("Processing population files and finding overlaps...\n")

# Initialize a summary data frame to track results across all populations
# This will store statistics for each population processed
summary_df <- tibble(
    population = character(),           # Population name (e.g., "pop1", "pop2")
    total_samples = integer(),          # Total samples in this population
    msl0_samples = integer(),          # Number of MSL0 samples in this population
    msl4_samples = integer(),          # Number of MSL4 samples in this population
    msl0_proportion = numeric(),       # Proportion of samples that are MSL0
    msl4_proportion = numeric()        # Proportion of samples that are MSL4
)

# Loop through each population file
for (pop_file in pop_files) {
    
    # Extract population number from filename using regex
    # This captures the digits after "pop" in the filename
    pop_num <- str_extract(basename(pop_file), "pop(\\d+)", group = 1)
    cat("\n--- Processing population", pop_num, "---\n")
    
    # Read population-specific sample list
    if (file.exists(pop_file)) {
        # Load samples for this population (same processing as MSL files)
        pop_samples <- read_lines(pop_file) %>%
            str_trim() %>%
            .[. != ""]  # Remove empty lines
        
        cat("Population", pop_num, "has", length(pop_samples), "samples\n")
    } else {
        # Skip this population if file doesn't exist
        cat("Warning: Population file not found:", pop_file, "\n")
        next  # Continue to next iteration of loop
    }
    
    # =============================================================================
    # FIND OVERLAPS BETWEEN POPULATION AND MSL STATUS
    # =============================================================================
    
    # Find samples that are both in this population AND have MSL0 status
    # intersect() returns elements present in both vectors
    msl0_overlap <- intersect(pop_samples, msl0_samples)
    cat("  MSL0 overlap:", length(msl0_overlap), "samples\n")
    
    # Find samples that are both in this population AND have MSL4 status
    msl4_overlap <- intersect(pop_samples, msl4_samples)
    cat("  MSL4 overlap:", length(msl4_overlap), "samples\n")
    
    # =============================================================================
    # CREATE OUTPUT FILE NAMES
    # =============================================================================
    
    # Define output file paths for this population's MSL-specific sample lists
    # These follow the naming convention specified in the requirements
    msl0_output <- file.path(pop_dir, paste0("Americas.K13.pop", pop_num, ".MSL0.diploidised.sample_names.txt"))
    msl4_output <- file.path(pop_dir, paste0("Americas.K13.pop", pop_num, ".MSL4.diploidised.sample_names.txt"))
    
    # =============================================================================
    # WRITE MSL-SPECIFIC SAMPLE LISTS
    # =============================================================================
    
    # Write MSL0 overlap samples to file
    if (length(msl0_overlap) > 0) {
        # write_lines() writes each element of the vector as a separate line
        write_lines(msl0_overlap, msl0_output)
        cat("  Created MSL0 file:", basename(msl0_output), "\n")
    } else {
        # Create empty file even if no samples found (for consistency)
        cat("  No MSL0 samples found for population", pop_num, "\n")
        write_lines(character(0), msl0_output)  # character(0) creates empty character vector
    }
    
    # Write MSL4 overlap samples to file
    if (length(msl4_overlap) > 0) {
        write_lines(msl4_overlap, msl4_output)
        cat("  Created MSL4 file:", basename(msl4_output), "\n")
    } else {
        cat("  No MSL4 samples found for population", pop_num, "\n")
        write_lines(character(0), msl4_output)
    }
    
    # =============================================================================
    # UPDATE SUMMARY STATISTICS
    # =============================================================================
    
    # Add results for this population to the summary data frame
    summary_df <- summary_df %>%
        add_row(
            population = paste0("pop", pop_num),                                    # Population identifier
            total_samples = length(pop_samples),                                    # Total samples in population
            msl0_samples = length(msl0_overlap),                                   # MSL0 samples in population
            msl4_samples = length(msl4_overlap),                                   # MSL4 samples in population
            msl0_proportion = round(length(msl0_overlap) / length(pop_samples), 3), # MSL0 proportion (3 decimal places)
            msl4_proportion = round(length(msl4_overlap) / length(pop_samples), 3)  # MSL4 proportion (3 decimal places)
        )
}

# =============================================================================
# SUMMARY STATISTICS AND REPORTING
# =============================================================================

cat("\n" , rep("=", 60), "\n")  # Create separator line with 60 equals signs
cat("=== SUMMARY STATISTICS ===\n")
cat(rep("=", 60), "\n\n")

# Print summary table to console
# This shows the distribution of MSL status across all populations
cat("Population-wise MSL distribution:\n")
print(summary_df)

# Write summary table to file for later analysis
summary_file <- file.path(pop_dir, "MSL_population_summary.tsv")
write_tsv(summary_df, summary_file)  # Tab-separated values format
cat("\nSummary table written to:", summary_file, "\n")

# =============================================================================
# CALCULATE OVERALL STATISTICS
# =============================================================================

# Calculate totals across all populations
# sum() adds up values across all rows in the specified column
total_pop_samples <- sum(summary_df$total_samples)      # Total samples across all populations
total_msl0_assigned <- sum(summary_df$msl0_samples)     # Total MSL0 samples assigned
total_msl4_assigned <- sum(summary_df$msl4_samples)     # Total MSL4 samples assigned

cat("\nOverall statistics:\n")
cat("  Total population samples processed:", total_pop_samples, "\n")
cat("  Total MSL0 samples assigned to populations:", total_msl0_assigned, "\n")
cat("  Total MSL4 samples assigned to populations:", total_msl4_assigned, "\n")
cat("  Overall MSL0 proportion:", round(total_msl0_assigned / total_pop_samples, 3), "\n")
cat("  Overall MSL4 proportion:", round(total_msl4_assigned / total_pop_samples, 3), "\n")

# =============================================================================
# IDENTIFY POPULATIONS WITH MSL BIAS
# =============================================================================

# Find populations with significant differences in MSL0 vs MSL4 proportions
# This could indicate population-specific evolutionary pressures on MSL
cat("\nPopulations with notable MSL bias (>20% difference):\n")

summary_df %>%
    # Calculate absolute difference between MSL0 and MSL4 proportions
    mutate(msl_difference = abs(msl0_proportion - msl4_proportion)) %>%
    # Filter for populations with >20% difference and sufficient sample size
    filter(
        msl_difference > 0.2,           # More than 20% difference
        total_samples >= 5              # At least 5 samples (avoid small-sample artifacts)
    ) %>%
    # Sort by difference (largest first)
    arrange(desc(msl_difference)) %>%
    # Select relevant columns for display
    select(population, total_samples, msl0_proportion, msl4_proportion, msl_difference) %>%
    print()

cat("\nAnalysis completed successfully!\n")
cat("Output files created in:", pop_dir, "\n")

# =============================================================================
# FILE VERIFICATION AND FINAL REPORTING
# =============================================================================

cat("\nCreated files:\n")

# Find all output files that were created
# Sys.glob() with pattern matching finds both MSL0 and MSL4 files
created_files <- c(
    Sys.glob(file.path(pop_dir, "Americas.K13.pop*.MSL0.diploidised.sample_names.txt")),
    Sys.glob(file.path(pop_dir, "Americas.K13.pop*.MSL4.diploidised.sample_names.txt"))
)

# Report on each created file
for (file in sort(created_files)) {
    # Count lines in each file to show how many samples it contains
    file_size <- length(read_lines(file))
    cat("  ", basename(file), "(", file_size, "samples )\n")
}

cat("\nTotal files created:", length(created_files), "\n")

# =============================================================================
# ADDITIONAL ANALYSIS SUGGESTIONS
# =============================================================================

# The following analyses could be performed with the generated data:
#
# 1. Statistical tests:
#    - Chi-square test for MSL distribution differences between populations
#    - Fisher's exact test for small sample sizes
#
# 2. Visualizations:
#    - Stacked bar charts showing MSL0/MSL4 proportions per population
#    - Scatter plots of MSL0 vs MSL4 proportions
#    - Heatmaps of MSL distribution patterns
#
# 3. Population genetics analyses:
#    - VCF subsetting using the generated sample lists
#    - Population-specific diversity statistics (π, Tajima's D)
#    - Hardy-Weinberg equilibrium tests within MSL groups
#    - FST calculations between MSL0 and MSL4 groups within populations
#
# 4. Evolutionary analyses:
#    - Selection tests comparing MSL0 vs MSL4 groups
#    - Phylogenetic analysis of MSL status evolution
#    - Geographic pattern analysis of MSL distribution