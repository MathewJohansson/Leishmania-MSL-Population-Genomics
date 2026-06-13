# Hardy-Weinberg Equilibrium Test for Polysomic Inheritance with ADMIXTURE Data - FINAL VERSION
# Using binary notation: 1 = normal allele, 0 = deletion allele
# For gene copy numbers: 0, 1, 2, 3, 4 (5 genotypes)
# Author: Zeynep Sakaoglu
# Date: 20/08/2025


# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(viridis)  # For colorblind-friendly palettes
library(RColorBrewer)

# Set colorblind-friendly palette
colorblind_palette <- c("#0173B2", "#DE8F05", "#029E73", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133")

# Function to convert coverage to binary genotype names
coverage_to_genotype <- function(coverage) {
  genotype_map <- c(
    "0" = "0000",  # 0 copies = homozygous deletion (all 0s)
    "1" = "1000",  # 1 copy = heterozygous with 1 normal allele
    "2" = "1100",  # 2 copies = heterozygous with 2 normal alleles  
    "3" = "1110",  # 3 copies = heterozygous with 3 normal alleles
    "4" = "1111"   # 4 copies = homozygous normal (all 1s)
  )
  return(genotype_map[as.character(coverage)])
}

# Function to extract population number for sorting
extract_pop_number <- function(pop_names) {
  # Extract numbers from population names (e.g., "v1", "v2", "v10")
  numbers <- as.numeric(gsub("[^0-9]", "", pop_names))
  # Handle cases where no number is found
  numbers[is.na(numbers)] <- 999
  return(numbers)
}

# Function to sort populations naturally (v1, v2, v3, ..., v10, v11, etc.)
sort_populations_naturally <- function(pop_names) {
  if(length(pop_names) == 0) return(pop_names)
  
  # Extract the numeric part and prefix
  pop_df <- data.frame(
    original = pop_names,
    prefix = gsub("[0-9]+$", "", pop_names),
    number = extract_pop_number(pop_names),
    stringsAsFactors = FALSE
  )
  
  # Sort by prefix first, then by number
  pop_df <- pop_df[order(pop_df$prefix, pop_df$number), ]
  
  return(pop_df$original)
}

# Function to create genotype count table from merged data
create_genotype_counts <- function(merged_data) {
  # Add genotype column using binary notation
  merged_data$genotype <- coverage_to_genotype(merged_data$coverage_rounded)
  
  # Create count table by population
  genotype_counts <- merged_data %>%
    group_by(population, genotype) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = genotype, values_from = count, values_fill = 0) %>%
    column_to_rownames("population")
  
  # Ensure all genotype columns exist in correct order (binary representation)
  required_genotypes <- c("0000", "1000", "1100", "1110", "1111")
  for(geno in required_genotypes) {
    if(!geno %in% colnames(genotype_counts)) {
      genotype_counts[[geno]] <- 0
    }
  }
  
  # Reorder columns to match expected order (1111 first for highest coverage)
  genotype_counts <- genotype_counts[, c("1111", "1110", "1100", "1000", "0000")]
  
  return(genotype_counts)
}

# Core Hardy-Weinberg test function (binary notation) - FIXED VERSION
polysomic_hwe_test <- function(genotype_counts, population_name = "Population") {
  if(length(genotype_counts) != 5) {
    stop("genotype_counts must have exactly 5 values: 1111, 1110, 1100, 1000, 0000")
  }
  
  names(genotype_counts) <- c("1111", "1110", "1100", "1000", "0000")
  total_n <- sum(genotype_counts)
  
  if(total_n == 0) {
    return(list(
      population = population_name,
      sample_size = 0,
      allele_freq_normal = NA,
      allele_freq_deletion = NA,
      chi_square = NA,
      df = NA,
      p_value = NA,
      significant = NA,
      test_reliable = FALSE,
      results_table = data.frame(),
      valid_genotypes = rep(FALSE, 5),
      test_performed = FALSE
    ))
  }
  
  # Also check if all individuals have the same genotype (no variation)
  non_zero_genotypes <- sum(genotype_counts > 0)
  if(non_zero_genotypes <= 1) {
    return(list(
      population = population_name,
      sample_size = total_n,
      allele_freq_normal = ifelse(genotype_counts["1111"] == total_n, 1, 
                                  ifelse(genotype_counts["0000"] == total_n, 0, NA)),
      allele_freq_deletion = ifelse(genotype_counts["1111"] == total_n, 0, 
                                    ifelse(genotype_counts["0000"] == total_n, 1, NA)),
      chi_square = NA,
      df = NA,
      p_value = NA,
      significant = NA,
      test_reliable = FALSE,
      results_table = data.frame(
        Genotype = names(genotype_counts),
        Coverage = c(4, 3, 2, 1, 0),
        Observed = as.numeric(genotype_counts),
        Expected_Freq = rep(NA, 5),
        Expected_Count = rep(NA, 5),
        Residual = rep(NA, 5),
        Standardized_Residual = rep(NA, 5),
        Valid_for_test = rep(FALSE, 5)
      ),
      valid_genotypes = rep(FALSE, 5),
      test_performed = FALSE
    ))
  }
  
  # Calculate allele frequencies (1 = normal allele, 0 = deletion allele)
  total_alleles <- total_n * 4
  normal_alleles <- 4*genotype_counts["1111"] + 3*genotype_counts["1110"] + 
    2*genotype_counts["1100"] + 1*genotype_counts["1000"]
  
  p <- normal_alleles / total_alleles  # frequency of normal allele (1)
  q <- 1 - p  # frequency of deletion allele (0)
  
  # Expected frequencies under HWE
  expected_freq <- c(
    "1111" = p^4,
    "1110" = 4 * p^3 * q,
    "1100" = 6 * p^2 * q^2,
    "1000" = 4 * p * q^3,
    "0000" = q^4
  )
  
  expected_counts <- expected_freq * total_n
  
  # Determine which genotypes are valid for testing
  if(total_n < 5) {
    # For very small samples, mark as unreliable but still calculate
    valid_genotypes <- rep(TRUE, 5)
    test_reliable <- FALSE
    warning(paste("Population", population_name, ": Very small sample size (n =", total_n, "). Results unreliable."))
  } else if(total_n < 10) {
    # For small samples, use all genotypes but mark as less reliable
    valid_genotypes <- rep(TRUE, 5)
    test_reliable <- FALSE
    warning(paste("Population", population_name, ": Small sample size (n =", total_n, "). Results should be interpreted cautiously."))
  } else {
    # For larger samples, use standard criteria
    valid_genotypes <- expected_counts >= 5
    test_reliable <- sum(valid_genotypes) >= 3
  }
  
  # Calculate chi-square using valid genotypes
  if(sum(valid_genotypes) >= 2 && total_n >= 2) {
    chi_square <- sum((genotype_counts[valid_genotypes] - expected_counts[valid_genotypes])^2 / 
                        pmax(expected_counts[valid_genotypes], 0.1))
    df <- sum(valid_genotypes) - 1 - 1  # -1 for constraint, -1 for estimated parameter
    
    if(df > 0) {
      p_value <- pchisq(chi_square, df, lower.tail = FALSE)
      significant <- !is.na(p_value) && p_value < 0.05
    } else {
      p_value <- NA
      significant <- NA
      warning(paste("Population", population_name, ": Insufficient degrees of freedom for chi-square test"))
    }
  } else {
    chi_square <- NA
    df <- NA
    p_value <- NA
    significant <- NA
    if(total_n >= 2) {
      warning(paste("Population", population_name, ": Cannot perform chi-square test - insufficient expected counts"))
    }
  }
  
  results_df <- data.frame(
    Genotype = names(genotype_counts),
    Coverage = c(4, 3, 2, 1, 0),
    Observed = as.numeric(genotype_counts),
    Expected_Freq = expected_freq,
    Expected_Count = expected_counts,
    Residual = genotype_counts - expected_counts,
    Standardized_Residual = (genotype_counts - expected_counts) / sqrt(pmax(expected_counts, 0.1)),
    Valid_for_test = valid_genotypes
  )
  
  summary_results <- list(
    population = population_name,
    sample_size = total_n,
    allele_freq_normal = p,  # frequency of normal allele (1)
    allele_freq_deletion = q,  # frequency of deletion allele (0)
    chi_square = chi_square,
    df = df,
    p_value = p_value,
    significant = significant,
    test_reliable = test_reliable,
    results_table = results_df,
    valid_genotypes = valid_genotypes
  )
  
  return(summary_results)
}

# FIXED: Function to test all populations from ADMIXTURE data
test_admixture_populations <- function(genotype_counts_df, min_sample_size = 1) {
  cat("=== TESTING ALL POPULATIONS ===\n")
  cat("Input genotype counts table:\n")
  print(genotype_counts_df)
  cat("\n")
  
  # Get all population names
  all_populations <- rownames(genotype_counts_df)
  pop_sizes <- rowSums(genotype_counts_df)
  
  cat("All population sizes:\n")
  print(pop_sizes)
  cat("\n")
  
  # Sort populations naturally
  sorted_populations <- sort_populations_naturally(all_populations)
  cat("Testing populations in order:", paste(sorted_populations, collapse = ", "), "\n\n")
  
  all_results <- list()
  summary_table <- data.frame()
  
  # Test ALL populations, regardless of size
  for(pop_name in sorted_populations) {
    cat("=== Testing", pop_name, "===\n")
    
    # Extract counts for this population
    pop_counts <- as.numeric(genotype_counts_df[pop_name, ])
    names(pop_counts) <- colnames(genotype_counts_df)
    
    result <- polysomic_hwe_test(pop_counts, pop_name)
    all_results[[pop_name]] <- result
    
    cat("Sample size:", result$sample_size, "\n")
    if(result$sample_size > 0) {
      cat("Normal allele frequency (1):", round(result$allele_freq_normal, 4), "\n")
      cat("Deletion allele frequency (0):", round(result$allele_freq_deletion, 4), "\n")
      if(!is.na(result$chi_square)) {
        cat("Chi-square:", round(result$chi_square, 4), "\n")
        cat("Degrees of freedom:", result$df, "\n")
        cat("P-value:", ifelse(is.na(result$p_value), "NA", format(result$p_value, scientific = TRUE)), "\n")
        cat("Significant deviation from HWE:", ifelse(is.na(result$significant), "NA", result$significant), "\n")
      } else {
        cat("Chi-square test: Not applicable\n")
      }
      cat("Test reliable:", result$test_reliable, "\n")
    } else {
      cat("No samples in this population\n")
    }
    cat("\n")
    
    # Add to summary table
    summary_table <- rbind(summary_table, data.frame(
      Population = pop_name,
      N = result$sample_size,
      p_normal = ifelse(is.na(result$allele_freq_normal), NA, round(result$allele_freq_normal, 4)),
      q_deletion = ifelse(is.na(result$allele_freq_deletion), NA, round(result$allele_freq_deletion, 4)),
      Chi_square = ifelse(is.na(result$chi_square), NA, round(result$chi_square, 4)),
      df = result$df,
      P_value = result$p_value,
      Significant = result$significant,
      Test_reliable = result$test_reliable,
      Test_performed = ifelse(is.null(result$test_performed), !is.na(result$p_value), result$test_performed)
    ))
  }
  
  # Sort summary table by population name (naturally)
  summary_table$Population <- factor(summary_table$Population, levels = sorted_populations)
  summary_table <- summary_table[order(summary_table$Population), ]
  
  cat("=== FINAL SUMMARY TABLE ===\n")
  print(summary_table)
  cat("\n")
  
  return(list(detailed_results = all_results, summary = summary_table))
}

# Function to process population dataset 
process_population_dataset <- function(pop_data, samples_data, dataset_name, min_proportion = 0.5) {
  cat("Processing population data for:", dataset_name, "\n")
  cat("Raw population data:", nrow(pop_data), "entries\n")
  
  # Filter samples with sufficient population assignment (above threshold)
  pop_data_filtered <- pop_data %>%
    filter(population_proportion >= min_proportion)
  
  cat("After filtering by proportion >=", min_proportion, ":", nrow(pop_data_filtered), "entries\n")
  
  # Merge with sample data using correct column names
  merged_data <- pop_data_filtered %>%
    left_join(samples_data, by = c("sample" = "sra_run_accession")) %>%
    filter(!is.na(coverage_rounded)) %>%
    filter(coverage_rounded %in% 0:4) # Ensure valid copy numbers
  
  cat("After merging and filtering:", nrow(merged_data), "samples from", length(unique(merged_data$population)), "populations\n")
  
  if(nrow(merged_data) > 0) {
    pop_sizes <- table(merged_data$population)
    cat("Population sizes:\n")
    print(pop_sizes)
  } else {
    cat("WARNING: No samples remaining after merging datasets!\n")
    cat("Check that sample IDs match between population and sample files\n")
  }
  
  return(merged_data)
}

# IMPROVED VISUALIZATION FUNCTIONS - COLORBLIND FRIENDLY

# Function to create individual HWE plots
plot_hwe_results <- function(hwe_result) {
  df <- hwe_result$results_table
  if(nrow(df) == 0 || hwe_result$sample_size == 0) {
    return(NULL)
  }
  
  df$Genotype <- factor(df$Genotype, levels = c("1111", "1110", "1100", "1000", "0000"))
  
  plot_data <- data.frame(
    Genotype = rep(df$Genotype, 2),
    Count = c(df$Observed, df$Expected_Count),
    Type = rep(c("Observed", "Expected"), each = nrow(df))
  )
  
  # Use colorblind-friendly colors
  ggplot(plot_data, aes(x = Genotype, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    labs(title = paste("Hardy-Weinberg Test:", hwe_result$population),
         subtitle = paste("N =", hwe_result$sample_size, 
                          ifelse(is.na(hwe_result$chi_square), "", 
                                 paste(", χ² =", round(hwe_result$chi_square, 3))),
                          ifelse(is.na(hwe_result$p_value), "", 
                                 paste(", p =", format(hwe_result$p_value, digits = 3)))),
         x = "Genotype (Binary representation)", 
         y = "Count",
         caption = "1111 = 4 copies, 1110 = 3 copies, 1100 = 2 copies, 1000 = 1 copy, 0000 = 0 copies\n1 = normal allele, 0 = deletion allele") +
    theme_minimal() +
    scale_fill_manual(values = c("Observed" = colorblind_palette[1], "Expected" = colorblind_palette[2])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 10))
}

# Function to create comprehensive overview plots - IMPROVED
create_population_overview_plots <- function(results_summary, dataset_name) {
  if(nrow(results_summary) == 0) {
    cat("No data to plot for", dataset_name, "\n")
    return(NULL)
  }
  
  # Clean and validate data
  results_summary <- results_summary[!is.na(results_summary$Population), ]
  if(nrow(results_summary) == 0) {
    cat("No valid population data for", dataset_name, "\n")
    return(NULL)
  }
  
  # Sort populations naturally
  sorted_pops <- sort_populations_naturally(as.character(results_summary$Population))
  results_summary$Population <- factor(results_summary$Population, levels = sorted_pops)
  results_summary <- results_summary[order(results_summary$Population), ]
  
  # Ensure numeric columns are properly formatted
  results_summary$N <- as.numeric(results_summary$N)
  results_summary$p_normal <- as.numeric(results_summary$p_normal)
  results_summary$P_value <- as.numeric(results_summary$P_value)
  
  plots_list <- list()
  
  # Plot 1: Sample sizes across all populations
  tryCatch({
    p1 <- ggplot(results_summary, aes(x = Population, y = N)) +
      geom_col(fill = colorblind_palette[1], alpha = 0.8, color = "black", size = 0.3) +
      geom_text(aes(label = N), vjust = -0.3, size = 3) +
      labs(title = paste("Sample Sizes -", dataset_name),
           x = "Population", y = "Sample Size (N)",
           subtitle = paste("Total populations:", nrow(results_summary))) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold")) +
      ylim(0, max(results_summary$N, na.rm = TRUE) * 1.15)
    plots_list$sample_sizes <- p1
  }, error = function(e) {
    cat("Error creating sample sizes plot:", e$message, "\n")
    plots_list$sample_sizes <- NULL
  })
  
  # Plot 2: Normal allele frequencies across all populations - COLORBLIND FRIENDLY
  tryCatch({
    # Create frequency categories for better visualization
    results_summary$freq_category <- cut(results_summary$p_normal, 
                                         breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                                         labels = c("Low (0-0.25)", "Medium-Low (0.25-0.5)", 
                                                    "Medium-High (0.5-0.75)", "High (0.75-1.0)"),
                                         include.lowest = TRUE)
    
    p2 <- ggplot(results_summary, aes(x = Population, y = p_normal, fill = freq_category)) +
      geom_col(alpha = 0.8, color = "black", size = 0.3) +
      geom_text(aes(label = round(p_normal, 3)), vjust = -0.3, size = 2.5) +
      scale_fill_manual(values = colorblind_palette[1:4], name = "Frequency\nCategory") +
      labs(title = paste("Normal Allele Frequencies -", dataset_name),
           x = "Population", y = "Frequency of Normal Allele (1)",
           subtitle = "Higher values = more normal alleles in population") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold")) +
      ylim(0, 1.1)
    plots_list$allele_frequencies <- p2
  }, error = function(e) {
    cat("Error creating allele frequencies plot:", e$message, "\n")
    plots_list$allele_frequencies <- NULL
  })
  
  # Plot 3: Hardy-Weinberg p-values - IMPROVED AND COLORBLIND FRIENDLY
  tryCatch({
    # Include ALL populations, even those without valid p-values
    results_summary$has_pvalue <- !is.na(results_summary$P_value) & 
      is.finite(results_summary$P_value) & 
      results_summary$P_value > 0
    
    # Create categories for display
    results_summary$HWE_Status <- ifelse(results_summary$N == 0, "No Test (Empty)",
                                         ifelse(!results_summary$Test_performed, "No Test (No Variation)",
                                                ifelse(is.na(results_summary$Significant), "No Test (Failed)",
                                                       ifelse(results_summary$Significant, "Deviates from HWE", "Fits HWE"))))
    
    # For plotting, use -log10(p-value) but handle missing values
    results_summary$log10_p <- ifelse(results_summary$has_pvalue, 
                                      -log10(pmax(results_summary$P_value, 1e-15)), 
                                      0)
    
    p3 <- ggplot(results_summary, aes(x = Population, y = log10_p, fill = HWE_Status)) +
      geom_col(alpha = 0.8, color = "black", size = 0.3) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
      geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "black", size = 1) +
      geom_text(aes(label = ifelse(has_pvalue, round(log10_p, 1), "NT")), 
                vjust = -0.3, size = 2.5) +
      scale_fill_manual(values = c("Fits HWE" = colorblind_palette[3], 
                                   "Deviates from HWE" = colorblind_palette[2],
                                   "No Test" = colorblind_palette[7]),
                        name = "HWE Status") +
      labs(title = paste("Hardy-Weinberg P-values -", dataset_name),
           x = "Population", y = "-log10(P-value)",
           subtitle = "All populations shown; NT = No Test possible",
           caption = "Black lines: p=0.05 (dashed), p=0.01 (dotted)\nHigher bars = more significant deviation from HWE") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold")) +
      coord_cartesian(ylim = c(0, max(15, max(results_summary$log10_p, na.rm = TRUE) * 1.1)))
    plots_list$hwe_pvalues <- p3
  }, error = function(e) {
    cat("Error creating p-values plot:", e$message, "\n")
    plots_list$hwe_pvalues <- NULL
  })
  
  # Plot 4: Combined overview - sample size with HWE status
  tryCatch({
    results_summary$Size_Category <- ifelse(results_summary$N == 0, "Empty",
                                            ifelse(results_summary$N < 5, "Very Small (<5)",
                                                   ifelse(results_summary$N < 10, "Small (5-9)",
                                                          ifelse(results_summary$N < 20, "Medium (10-19)", "Large (≥20)"))))
    
    # Debug: Print HWE status distribution
    cat("HWE Status distribution for", dataset_name, ":\n")
    print(table(results_summary$HWE_Status, useNA = "ifany"))
    cat("Sample sizes for each HWE status:\n")
    for(status in unique(results_summary$HWE_Status)) {
      cat(status, ": N =", paste(results_summary$N[results_summary$HWE_Status == status], collapse = ", "), "\n")
    }
    
    # Calculate triangle positions to avoid overlap
    max_height <- max(results_summary$N, na.rm = TRUE)
    triangle_y <- ifelse(max_height > 0, max_height * 0.85, 5)
    
    p4 <- ggplot(results_summary, aes(x = Population, y = N)) +
      geom_col(aes(fill = Size_Category), alpha = 0.7, color = "black", size = 0.3) +
      geom_point(aes(color = HWE_Status), y = triangle_y, size = 5, shape = 17, stroke = 1) +
      scale_fill_manual(values = c("Empty" = colorblind_palette[7],
                                   "Very Small (<5)" = colorblind_palette[8], 
                                   "Small (5-9)" = colorblind_palette[6], 
                                   "Medium (10-19)" = colorblind_palette[4], 
                                   "Large (≥20)" = colorblind_palette[3]),
                        name = "Sample Size",
                        drop = FALSE) +
      scale_color_manual(values = c("Fits HWE" = colorblind_palette[1], 
                                    "Deviates from HWE" = colorblind_palette[2], 
                                    "No Test (Empty)" = colorblind_palette[7],
                                    "No Test (No Variation)" = colorblind_palette[7],
                                    "No Test (Failed)" = colorblind_palette[7]),
                         name = "HWE Test",
                         drop = FALSE) +
      labs(title = paste("Population Overview -", dataset_name),
           x = "Population", y = "Sample Size",
           subtitle = paste("Bars = sample size, Triangles = HWE test results",
                            "\nTriangle counts: Fits HWE =", sum(results_summary$HWE_Status == "Fits HWE", na.rm = TRUE),
                            ", Deviates =", sum(results_summary$HWE_Status == "Deviates from HWE", na.rm = TRUE),
                            ", No Test =", sum(results_summary$HWE_Status == "No Test", na.rm = TRUE))) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold"))
    plots_list$summary <- p4
  }, error = function(e) {
    cat("Error creating summary plot:", e$message, "\n")
    plots_list$summary <- NULL
  })
  
  return(plots_list)
}

# MAIN ANALYSIS WITH IMPROVED DATA PROCESSING

cat("=== HARDY-WEINBERG TEST FOR POLYSOMIC INHERITANCE - FINAL VERSION ===\n")
cat("Using binary notation: 1 = normal allele, 0 = deletion allele\n")
cat("Improvements: All populations tested, colorblind-friendly plots, natural sorting\n\n")

# STEP 1: Load your actual data files
cat("STEP 1: Loading your data files\n")

# Load sample data with coverage information
samples_data <- read_tsv("data/vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv", show_col_types = FALSE)
cat("Loaded samples data:", nrow(samples_data), "samples\n")

# Load population membership tables
pop_membership_Americas_k13 <- read_tsv("data/population_membership_table.Americas_only.k13.tsv", show_col_types = FALSE)
pop_membership_Infantum_k13 <- read_tsv("data/population_membership_table.infantum_only.k13.tsv", show_col_types = FALSE)
pop_membership_All_k12 <- read_tsv("data/population_membership_table.allsamples.k12.tsv", show_col_types = FALSE)
pop_membership_Infantum_k8 <- read_tsv("data/population_membership_table.infantum_only.k8.tsv", show_col_types = FALSE)

cat("Loaded all population membership files\n\n")

# Check the structure of your samples data
cat("=== SAMPLES DATA STRUCTURE ===\n")
cat("Sample data dimensions:", nrow(samples_data), "rows x", ncol(samples_data), "columns\n")
cat("Coverage data distribution:\n")
print(table(samples_data$coverage_rounded, useNA = "ifany"))
cat("\n")

# STEP 2: Process each ADMIXTURE result
cat("STEP 2: Processing ADMIXTURE files with improved analysis\n\n")

# Create a list of population datasets with descriptive names
pop_datasets <- list(
  "Americas_k13" = pop_membership_Americas_k13,
  "Infantum_k13" = pop_membership_Infantum_k13,
  "AllSamples_k12" = pop_membership_All_k12,
  "Infantum_k8" = pop_membership_Infantum_k8
)

# Store results and plots for comparison
all_dataset_results <- list()
all_dataset_plots <- list()

# Process each dataset with improved methodology
for(dataset_name in names(pop_datasets)) {
  cat("\n", rep("=", 60), "\n")
  cat("--- Processing", dataset_name, "---\n")
  
  pop_data <- pop_datasets[[dataset_name]]
  cat("Population data loaded:", nrow(pop_data), "entries\n")
  
  # Process this population dataset
  merged_data <- process_population_dataset(pop_data, samples_data, dataset_name, min_proportion = 0.5)
  
  if(nrow(merged_data) == 0) {
    cat("No valid data found after merging, skipping...\n")
    next
  }
  
  # Create genotype counts
  genotype_counts <- create_genotype_counts(merged_data)
  cat("Genotype count table for", dataset_name, ":\n")
  print(genotype_counts)
  
  # Perform Hardy-Weinberg tests - NOW TESTING ALL POPULATIONS
  cat("\n=== TESTING ALL POPULATIONS (including small ones) ===\n")
  results <- test_admixture_populations(genotype_counts, min_sample_size = 0)  # Changed to 0 to test ALL
  
  # Store results
  all_dataset_results[[dataset_name]] <- results
  
  cat("\n=== SUMMARY TABLE FOR", dataset_name, "===\n")
  if(nrow(results$summary) > 0) {
    print(results$summary)
    
    # Create comprehensive population overview plots with improvements
    cat("Creating improved population plots for", dataset_name, "...\n")
    overview_plots <- create_population_overview_plots(results$summary, dataset_name)
    
    # Store plots for saving later
    all_dataset_plots[[dataset_name]] <- list()
    all_dataset_plots[[dataset_name]]$overview <- overview_plots
    all_dataset_plots[[dataset_name]]$detailed <- list()
    
    if(!is.null(overview_plots)) {
      # Print plots that were successfully created
      if(!is.null(overview_plots$sample_sizes)) {
        print(overview_plots$sample_sizes)
      }
      if(!is.null(overview_plots$allele_frequencies)) {
        print(overview_plots$allele_frequencies)
      }
      if(!is.null(overview_plots$hwe_pvalues)) {
        print(overview_plots$hwe_pvalues)
      }
      if(!is.null(overview_plots$summary)) {
        print(overview_plots$summary)
      }
    }
    
    # Create individual population plots for populations with data
    tryCatch({
      populations_with_data <- results$summary[results$summary$N > 0, "Population"]
      if(length(populations_with_data) > 0) {
        cat("Creating detailed plots for", length(populations_with_data), "populations with data...\n")
        for(pop_name in populations_with_data) {
          if(pop_name %in% names(results$detailed_results)) {
            tryCatch({
              plot_obj <- plot_hwe_results(results$detailed_results[[pop_name]])
              if(!is.null(plot_obj)) {
                print(plot_obj)
                # Store for saving
                all_dataset_plots[[dataset_name]]$detailed[[pop_name]] <- plot_obj
              }
            }, error = function(e) {
              cat("Error creating plot for population", pop_name, ":", e$message, "\n")
            })
          }
        }
      }
    }, error = function(e) {
      cat("Error in individual population plotting:", e$message, "\n")
    })
  } else {
    cat("No populations found in dataset\n")
  }
}

cat("\n", rep("=", 60), "\n")
cat("=== ANALYSIS COMPLETE ===\n")

# Final summary across all datasets
cat("\n=== CROSS-DATASET COMPARISON ===\n")
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    total_populations <- nrow(results$summary)
    populations_with_data <- sum(results$summary$N > 0, na.rm = TRUE)
    tested_populations <- sum(!is.na(results$summary$P_value), na.rm = TRUE)
    significant_count <- sum(results$summary$Significant %in% TRUE, na.rm = TRUE)
    
    cat(dataset_name, ":\n")
    cat("  Total populations: ", total_populations, "\n")
    cat("  Populations with data: ", populations_with_data, "\n")
    cat("  Populations tested: ", tested_populations, "\n")
    cat("  Significant HWE deviations: ", significant_count, "\n\n")
  }
}

# ===== SAVE ALL PLOTS TO HWE_FINAL DIRECTORY =====
cat("\n=== SAVING PLOTS TO HWE_FINAL DIRECTORY ===\n")

# Create output directory if it doesn't exist
if(!dir.exists("HWE_FINAL")) {
  dir.create("HWE_FINAL")
  cat("Created directory: HWE_FINAL/\n")
} else {
  cat("Using existing directory: HWE_FINAL/\n")
}

# Function to safely save plots with improved error handling
safe_ggsave <- function(plot_obj, filename, width = 12, height = 8) {
  if(!is.null(plot_obj)) {
    tryCatch({
      ggsave(filename, plot = plot_obj, width = width, height = height, dpi = 300)
      cat("✓ Saved:", filename, "\n")
      return(TRUE)
    }, error = function(e) {
      cat("✗ Error saving", filename, ":", e$message, "\n")
      return(FALSE)
    })
  } else {
    cat("⚠ Skipped", filename, "(plot object is NULL)\n")
    return(FALSE)
  }
}

# Save plots from stored objects
saved_count <- 0
for(dataset_name in names(all_dataset_plots)) {
  cat("\n--- Saving plots for", dataset_name, "---\n")
  
  dataset_plots <- all_dataset_plots[[dataset_name]]
  
  # Save overview plots
  if(!is.null(dataset_plots$overview)) {
    overview <- dataset_plots$overview
    
    if(safe_ggsave(overview$sample_sizes, 
                   paste0("HWE_FINAL/", dataset_name, "_sample_sizes.png"), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$allele_frequencies, 
                   paste0("HWE_FINAL/", dataset_name, "_allele_frequencies.png"), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$hwe_pvalues, 
                   paste0("HWE_FINAL/", dataset_name, "_hwe_pvalues.png"), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$summary, 
                   paste0("HWE_FINAL/", dataset_name, "_overview.png"), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
  }
  
  # Save individual population plots (limit to first 20 to avoid too many files)
  if(!is.null(dataset_plots$detailed) && length(dataset_plots$detailed) > 0) {
    plot_names <- names(dataset_plots$detailed)
    # Sort population names naturally for consistent ordering
    plot_names <- sort_populations_naturally(plot_names)
    # Limit to first 20 plots to avoid overwhelming number of files
    plot_names <- plot_names[1:min(20, length(plot_names))]
    
    cat("Saving detailed plots for first", length(plot_names), "populations...\n")
    for(pop_name in plot_names) {
      if(safe_ggsave(dataset_plots$detailed[[pop_name]], 
                     paste0("HWE_FINAL/", dataset_name, "_", pop_name, "_detailed.png"),
                     width = 10, height = 6)) {
        saved_count <- saved_count + 1
      }
    }
  }
}

# Save summary tables as CSV files with improved formatting
cat("\n--- Saving summary tables ---\n")
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    tryCatch({
      # Create a clean summary table for export
      clean_summary <- results$summary
      clean_summary$P_value_formatted <- ifelse(is.na(clean_summary$P_value), "NA",
                                                ifelse(clean_summary$P_value < 0.001, 
                                                       format(clean_summary$P_value, scientific = TRUE, digits = 3),
                                                       format(round(clean_summary$P_value, 4), nsmall = 4)))
      
      filename <- paste0("HWE_FINAL/", dataset_name, "_summary_table.csv")
      write.csv(clean_summary, filename, row.names = FALSE)
      cat("✓ Saved:", filename, "\n")
      saved_count <- saved_count + 1
    }, error = function(e) {
      cat("✗ Error saving summary table for", dataset_name, ":", e$message, "\n")
    })
  }
}

# Create a comprehensive combined summary table across all datasets
cat("\n--- Creating comprehensive combined summary ---\n")
combined_summary <- data.frame()
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    temp_summary <- results$summary
    temp_summary$Dataset <- dataset_name
    # Reorder columns to put Dataset first
    temp_summary <- temp_summary[, c("Dataset", setdiff(names(temp_summary), "Dataset"))]
    combined_summary <- rbind(combined_summary, temp_summary)
  }
}

if(nrow(combined_summary) > 0) {
  # Sort combined summary by dataset and then by population (naturally)
  combined_summary$Population <- as.character(combined_summary$Population)
  combined_summary <- combined_summary[order(combined_summary$Dataset, 
                                             extract_pop_number(combined_summary$Population)), ]
  
  # Save combined summary
  tryCatch({
    write.csv(combined_summary, "HWE_FINAL/combined_summary_all_datasets.csv", row.names = FALSE)
    cat("✓ Saved: HWE_FINAL/combined_summary_all_datasets.csv\n")
    saved_count <- saved_count + 1
  }, error = function(e) {
    cat("✗ Error saving combined summary:", e$message, "\n")
  })
  
  # Create final comparison plots across all datasets - COLORBLIND FRIENDLY
  tryCatch({
    # Comparison plot 1: Allele frequencies across datasets
    comparison_plot1 <- ggplot(combined_summary[combined_summary$N > 0, ], 
                               aes(x = Population, y = p_normal, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.2) +
      scale_fill_manual(values = colorblind_palette[1:length(unique(combined_summary$Dataset))]) +
      labs(title = "Normal Allele Frequencies Across All Datasets",
           x = "Population", y = "Frequency of Normal Allele (1)",
           subtitle = "Comparison of all ADMIXTURE results - only populations with data shown") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 16, face = "bold")) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 2)
    
    if(safe_ggsave(comparison_plot1, "HWE_FINAL/all_datasets_allele_frequencies_comparison.png", 
                   width = 18, height = 12)) {
      saved_count <- saved_count + 1
    }
    
    # Comparison plot 2: Sample sizes across datasets
    comparison_plot2 <- ggplot(combined_summary, aes(x = Population, y = N, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.2) +
      scale_fill_manual(values = colorblind_palette[1:length(unique(combined_summary$Dataset))]) +
      labs(title = "Sample Sizes Across All Datasets",
           x = "Population", y = "Sample Size (N)",
           subtitle = "Comparison of population sizes in all ADMIXTURE results") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 16, face = "bold")) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 2)
    
    if(safe_ggsave(comparison_plot2, "HWE_FINAL/all_datasets_sample_sizes_comparison.png", 
                   width = 18, height = 12)) {
      saved_count <- saved_count + 1
    }
    
    # Comparison plot 3: HWE test results across datasets
    hwe_summary <- combined_summary[!is.na(combined_summary$P_value) & combined_summary$N > 0, ]
    if(nrow(hwe_summary) > 0) {
      hwe_summary$log10_p <- -log10(pmax(hwe_summary$P_value, 1e-15))
      hwe_summary$HWE_Status <- ifelse(hwe_summary$Significant, "Deviates from HWE", "Fits HWE")
      
      comparison_plot3 <- ggplot(hwe_summary, aes(x = Population, y = log10_p, fill = HWE_Status)) +
        geom_col(alpha = 0.8, color = "black", size = 0.2) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
        scale_fill_manual(values = c("Fits HWE" = colorblind_palette[3], 
                                     "Deviates from HWE" = colorblind_palette[2])) +
        labs(title = "Hardy-Weinberg Test Results Across All Datasets",
             x = "Population", y = "-log10(P-value)",
             subtitle = "Only populations with valid test results shown",
             caption = "Black dashed line: p = 0.05 significance threshold") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 16, face = "bold")) +
        facet_wrap(~Dataset, scales = "free_x", ncol = 2)
      
      if(safe_ggsave(comparison_plot3, "HWE_FINAL/all_datasets_hwe_results_comparison.png", 
                     width = 18, height = 12)) {
        saved_count <- saved_count + 1
      }
    }
    
  }, error = function(e) {
    cat("✗ Error creating comparison plots:", e$message, "\n")
  })
}

# Create a final statistics summary
cat("\n--- Creating final statistics summary ---\n")
tryCatch({
  stats_summary <- data.frame()
  for(dataset_name in names(all_dataset_results)) {
    results <- all_dataset_results[[dataset_name]]
    if(nrow(results$summary) > 0) {
      total_pops <- nrow(results$summary)
      pops_with_data <- sum(results$summary$N > 0, na.rm = TRUE)
      pops_tested <- sum(!is.na(results$summary$P_value), na.rm = TRUE)
      pops_significant <- sum(results$summary$Significant %in% TRUE, na.rm = TRUE)
      avg_sample_size <- mean(results$summary$N[results$summary$N > 0], na.rm = TRUE)
      avg_normal_freq <- mean(results$summary$p_normal[results$summary$N > 0], na.rm = TRUE)
      
      stats_summary <- rbind(stats_summary, data.frame(
        Dataset = dataset_name,
        Total_Populations = total_pops,
        Populations_With_Data = pops_with_data,
        Populations_Tested = pops_tested,
        Significant_Deviations = pops_significant,
        Percent_Significant = round(100 * pops_significant / pmax(pops_tested, 1), 1),
        Avg_Sample_Size = round(avg_sample_size, 1),
        Avg_Normal_Allele_Freq = round(avg_normal_freq, 3)
      ))
    }
  }
  
  if(nrow(stats_summary) > 0) {
    write.csv(stats_summary, "HWE_FINAL/dataset_statistics_summary.csv", row.names = FALSE)
    cat("✓ Saved: HWE_FINAL/dataset_statistics_summary.csv\n")
    saved_count <- saved_count + 1
    
    cat("\n=== FINAL DATASET STATISTICS ===\n")
    print(stats_summary)
  }
}, error = function(e) {
  cat("✗ Error creating statistics summary:", e$message, "\n")
})

cat("\n=== PLOT SAVING COMPLETE ===\n")
cat("Total files saved:", saved_count, "\n")
cat("All files saved to: HWE_FINAL/ directory\n")

# List all saved files
if(dir.exists("HWE_FINAL")) {
  plot_files <- list.files("HWE_FINAL", full.names = FALSE)
  if(length(plot_files) > 0) {
    cat("\nFiles in HWE_FINAL directory:\n")
    # Group files by type
    png_files <- plot_files[grepl("\\.png$", plot_files)]
    csv_files <- plot_files[grepl("\\.csv$", plot_files)]
    
    if(length(png_files) > 0) {
      cat("\nPlot files (", length(png_files), "):\n")
      for(file in sort(png_files)) {
        cat(" -", file, "\n")
      }
    }
    
    if(length(csv_files) > 0) {
      cat("\nData files (", length(csv_files), "):\n")
      for(file in sort(csv_files)) {
        cat(" -", file, "\n")
      }
    }
  } else {
    cat("\nNo files found in HWE_FINAL directory!\n")
  }
} else {
  cat("\nHWE_FINAL directory not found!\n")
}

#--------------------------PART 2-------------------------------------------------------------
#--------------------------MOCK DATA TESTING--------------------------------------------------

# Mock Data Generator for Hardy-Weinberg Testing in Tetraploid Systems

library(dplyr)
library(readr)

set.seed(12345)  

# Function to generate tetraploid genotype counts given allele frequency
generate_hwe_counts <- function(p, total_n, population_name) {
  # p = frequency of normal allele (1)
  # q = frequency of deletion allele (0)
  q <- 1 - p
  
  # Expected frequencies under HWE for tetraploid
  expected_freq <- c(
    "1111" = p^4,                    # 4 copies normal
    "1110" = 4 * p^3 * q,           # 3 copies normal, 1 deletion
    "1100" = 6 * p^2 * q^2,         # 2 copies normal, 2 deletions
    "1000" = 4 * p * q^3,           # 1 copy normal, 3 deletions
    "0000" = q^4                     # 4 copies deletion
  )
  
  # Generate multinomial counts
  observed_counts <- rmultinom(1, total_n, expected_freq)[,1]
  names(observed_counts) <- c("4", "3", "2", "1", "0")  # coverage levels
  
  return(data.frame(
    population = population_name,
    coverage_4 = observed_counts["4"],
    coverage_3 = observed_counts["3"],
    coverage_2 = observed_counts["2"],
    coverage_1 = observed_counts["1"],
    coverage_0 = observed_counts["0"],
    total_n = total_n,
    expected_p = p,
    status = "HWE"
  ))
}

# Function to generate NON-HWE counts (various deviation patterns)
generate_non_hwe_counts <- function(p, total_n, population_name, deviation_type = "excess_hetero") {
  q <- 1 - p
  
  # Start with HWE expected frequencies
  expected_freq <- c(
    "1111" = p^4,
    "1110" = 4 * p^3 * q,
    "1100" = 6 * p^2 * q^2,
    "1000" = 4 * p * q^3,
    "0000" = q^4
  )
  
  # Apply different types of deviations
  if(deviation_type == "excess_hetero") {
    # Excess heterozygotes (common in population structure)
    # Reduce homozygotes, increase heterozygotes
    modified_freq <- expected_freq
    modified_freq["1111"] <- modified_freq["1111"] * 0.6  # Reduce extreme homozygotes
    modified_freq["0000"] <- modified_freq["0000"] * 0.6
    modified_freq["1110"] <- modified_freq["1110"] * 1.4  # Increase heterozygotes
    modified_freq["1100"] <- modified_freq["1100"] * 1.3
    modified_freq["1000"] <- modified_freq["1000"] * 1.4
    
  } else if(deviation_type == "deficit_hetero") {
    # Deficit of heterozygotes (inbreeding)
    # Increase homozygotes, decrease heterozygotes
    modified_freq <- expected_freq
    modified_freq["1111"] <- modified_freq["1111"] * 1.5
    modified_freq["0000"] <- modified_freq["0000"] * 1.5
    modified_freq["1110"] <- modified_freq["1110"] * 0.5
    modified_freq["1100"] <- modified_freq["1100"] * 0.7
    modified_freq["1000"] <- modified_freq["1000"] * 0.5
    
  } else if(deviation_type == "allele_bias") {
    # Systematic bias toward one allele (technical artifacts)
    modified_freq <- expected_freq
    modified_freq["1111"] <- modified_freq["1111"] * 2.0
    modified_freq["1110"] <- modified_freq["1110"] * 1.3
    modified_freq["0000"] <- modified_freq["0000"] * 0.3
    modified_freq["1000"] <- modified_freq["1000"] * 0.7
  }
  
  # Normalize to sum to 1
  modified_freq <- modified_freq / sum(modified_freq)
  
  # Generate multinomial counts
  observed_counts <- rmultinom(1, total_n, modified_freq)[,1]
  names(observed_counts) <- c("4", "3", "2", "1", "0")
  
  return(data.frame(
    population = population_name,
    coverage_4 = observed_counts["4"],
    coverage_3 = observed_counts["3"],
    coverage_2 = observed_counts["2"],
    coverage_1 = observed_counts["1"],
    coverage_0 = observed_counts["0"],
    total_n = total_n,
    expected_p = p,
    status = paste("NON_HWE", deviation_type, sep = "_")
  ))
}

# Function to generate populations with no variation (monomorphic)
generate_monomorphic_counts <- function(coverage_level, total_n, population_name) {
  counts <- c(0, 0, 0, 0, 0)
  names(counts) <- c("4", "3", "2", "1", "0")
  counts[as.character(coverage_level)] <- total_n
  
  return(data.frame(
    population = population_name,
    coverage_4 = counts["4"],
    coverage_3 = counts["3"],
    coverage_2 = counts["2"],
    coverage_1 = counts["1"],
    coverage_0 = counts["0"],
    total_n = total_n,
    expected_p = ifelse(coverage_level == 4, 1.0, 
                        ifelse(coverage_level == 0, 0.0, coverage_level/4)),
    status = "MONOMORPHIC"
  ))
}

cat("=== GENERATING MOCK DATA FOR HWE TESTING ===\n\n")

# Create comprehensive test dataset
mock_populations <- data.frame()

# 1. POPULATIONS IN HWE (should pass HWE test)
cat("Creating populations in Hardy-Weinberg equilibrium...\n")

# Various allele frequencies and sample sizes
hwe_scenarios <- list(
  list(p = 0.1, n = 50, name = "v1"),   # Low frequency, medium sample
  list(p = 0.3, n = 100, name = "v2"),  # Low-medium frequency, large sample
  list(p = 0.5, n = 80, name = "v3"),   # Balanced frequency, medium sample
  list(p = 0.7, n = 60, name = "v4"),   # High frequency, medium sample
  list(p = 0.9, n = 40, name = "v5"),   # Very high frequency, small sample
  list(p = 0.2, n = 25, name = "v6"),   # Low frequency, small sample
  list(p = 0.6, n = 120, name = "v7"),  # Medium-high frequency, large sample
  list(p = 0.8, n = 35, name = "v8")    # High frequency, small sample
)

for(scenario in hwe_scenarios) {
  pop_data <- generate_hwe_counts(scenario$p, scenario$n, scenario$name)
  mock_populations <- rbind(mock_populations, pop_data)
}

# 2. POPULATIONS NOT IN HWE (should fail HWE test)
cat("Creating populations deviating from Hardy-Weinberg equilibrium...\n")

# Excess heterozygotes scenarios
non_hwe_excess <- list(
  list(p = 0.3, n = 50, name = "v9", type = "excess_hetero"),
  list(p = 0.6, n = 75, name = "v10", type = "excess_hetero"),
  list(p = 0.4, n = 90, name = "v11", type = "excess_hetero")
)

for(scenario in non_hwe_excess) {
  pop_data <- generate_non_hwe_counts(scenario$p, scenario$n, scenario$name, scenario$type)
  mock_populations <- rbind(mock_populations, pop_data)
}

# Deficit heterozygotes scenarios (inbreeding-like)
non_hwe_deficit <- list(
  list(p = 0.5, n = 60, name = "v12", type = "deficit_hetero"),
  list(p = 0.3, n = 45, name = "v13", type = "deficit_hetero"),
  list(p = 0.7, n = 55, name = "v14", type = "deficit_hetero")
)

for(scenario in non_hwe_deficit) {
  pop_data <- generate_non_hwe_counts(scenario$p, scenario$n, scenario$name, scenario$type)
  mock_populations <- rbind(mock_populations, pop_data)
}

# Allele bias scenarios
non_hwe_bias <- list(
  list(p = 0.4, n = 70, name = "v15", type = "allele_bias"),
  list(p = 0.6, n = 40, name = "v16", type = "allele_bias")
)

for(scenario in non_hwe_bias) {
  pop_data <- generate_non_hwe_counts(scenario$p, scenario$n, scenario$name, scenario$type)
  mock_populations <- rbind(mock_populations, pop_data)
}

# 3. EDGE CASES
cat("Creating edge case populations...\n")

# Monomorphic populations (no variation)
edge_cases <- list(
  list(coverage = 4, n = 20, name = "v17"),  # All normal alleles
  list(coverage = 0, n = 15, name = "v18"),  # All deletion alleles
  list(coverage = 2, n = 10, name = "v19")   # All intermediate
)

for(scenario in edge_cases) {
  pop_data <- generate_monomorphic_counts(scenario$coverage, scenario$n, scenario$name)
  mock_populations <- rbind(mock_populations, pop_data)
}

# Very small samples
small_samples <- list(
  list(p = 0.5, n = 3, name = "v20"),
  list(p = 0.3, n = 1, name = "v21"),
  list(p = 0.7, n = 2, name = "v22")
)

for(scenario in small_samples) {
  pop_data <- generate_hwe_counts(scenario$p, scenario$n, scenario$name)
  pop_data$status <- "SMALL_SAMPLE"
  mock_populations <- rbind(mock_populations, pop_data)
}

# Empty population
empty_pop <- data.frame(
  population = "v23",
  coverage_4 = 0, coverage_3 = 0, coverage_2 = 0, coverage_1 = 0, coverage_0 = 0,
  total_n = 0, expected_p = NA, status = "EMPTY"
)
mock_populations <- rbind(mock_populations, empty_pop)

# Display summary of created data
cat("\n=== MOCK DATA SUMMARY ===\n")
cat("Total populations created:", nrow(mock_populations), "\n\n")

status_summary <- table(mock_populations$status)
cat("Population types:\n")
print(status_summary)

cat("\nSample size distribution:\n")
print(table(mock_populations$total_n))

cat("\nFirst few rows of mock data:\n")
print(head(mock_populations, 10))

# Create the individual sample data (expand populations to individual samples)
cat("\n=== CREATING INDIVIDUAL SAMPLE DATA ===\n")

create_individual_samples <- function(pop_row) {
  pop_name <- pop_row$population
  
  # Create individual samples based on genotype counts
  samples <- c()
  
  # Add samples for each coverage level
  if(pop_row$coverage_4 > 0) {
    samples <- c(samples, rep(4, pop_row$coverage_4))
  }
  if(pop_row$coverage_3 > 0) {
    samples <- c(samples, rep(3, pop_row$coverage_3))
  }
  if(pop_row$coverage_2 > 0) {
    samples <- c(samples, rep(2, pop_row$coverage_2))
  }
  if(pop_row$coverage_1 > 0) {
    samples <- c(samples, rep(1, pop_row$coverage_1))
  }
  if(pop_row$coverage_0 > 0) {
    samples <- c(samples, rep(0, pop_row$coverage_0))
  }
  
  # Create sample IDs
  if(length(samples) > 0) {
    sample_ids <- paste0(pop_name, "_sample_", 1:length(samples))
    return(data.frame(
      sra_run_accession = sample_ids,
      coverage_rounded = samples,
      population_source = pop_name,
      expected_status = pop_row$status
    ))
  } else {
    return(data.frame(
      sra_run_accession = character(0),
      coverage_rounded = numeric(0),
      population_source = character(0),
      expected_status = character(0)
    ))
  }
}

# Generate individual sample data
individual_samples <- data.frame()
for(i in 1:nrow(mock_populations)) {
  pop_samples <- create_individual_samples(mock_populations[i, ])
  individual_samples <- rbind(individual_samples, pop_samples)
}

cat("Total individual samples created:", nrow(individual_samples), "\n")
cat("Coverage distribution in samples:\n")
print(table(individual_samples$coverage_rounded))

# Create population membership data for ADMIXTURE format
cat("\n=== CREATING ADMIXTURE-STYLE POPULATION MEMBERSHIP DATA ===\n")

# Create population membership with high confidence (proportion = 0.8)
pop_membership <- data.frame(
  sample = individual_samples$sra_run_accession,
  population = individual_samples$population_source,
  population_proportion = 0.8,  # High confidence assignment
  expected_hwe_status = individual_samples$expected_status
)

cat("Population membership data created:", nrow(pop_membership), "entries\n")
cat("Populations in membership data:", length(unique(pop_membership$population)), "\n")

# Save the mock datasets
cat("\n=== SAVING MOCK DATA FILES ===\n")

# Create mock_data directory
if(!dir.exists("mock_data")) {
  dir.create("mock_data")
  cat("Created directory: mock_data/\n")
}

# Save sample data (equivalent to your vcf_samples file)
write_tsv(individual_samples, "mock_data/mock_vcf_samples_data.tsv")
cat("✓ Saved: mock_data/mock_vcf_samples_data.tsv\n")

# Save population membership (equivalent to your ADMIXTURE files)
write_tsv(pop_membership, "mock_data/mock_population_membership.tsv")
cat("✓ Saved: mock_data/mock_population_membership.tsv\n")

# Save population summary for reference
write_tsv(mock_populations, "mock_data/mock_population_summary.tsv")
cat("✓ Saved: mock_data/mock_population_summary.tsv\n")

# Create a test script to run with your HWE code
test_script <- 

cat("=== TESTING HWE ANALYSIS WITH MOCK DATA ===\\n")

# Load mock data
samples_data <- read_tsv("mock_data/mock_vcf_samples_data.tsv", show_col_types = FALSE)
pop_membership_mock <- read_tsv("mock_data/mock_population_membership.tsv", show_col_types = FALSE)

cat("Loaded mock data:\\n")
cat("- Samples:", nrow(samples_data), "\\n")
cat("- Population memberships:", nrow(pop_membership_mock), "\\n")

# Process mock data using your existing functions
merged_data <- process_population_dataset(pop_membership_mock, samples_data, "MockData", min_proportion = 0.5)

if(nrow(merged_data) > 0) {
  # Create genotype counts
  genotype_counts <- create_genotype_counts(merged_data)
  cat("\\nGenotype count table for mock data:\\n")
  print(genotype_counts)
  
  # Perform Hardy-Weinberg tests
  results <- test_admixture_populations(genotype_counts, min_sample_size = 0)
  
  cat("\\n=== MOCK DATA TEST RESULTS ===\\n")
  if(nrow(results$summary) > 0) {
    # Add expected results for comparison
    expected_results <- read_tsv("mock_data/mock_population_summary.tsv", show_col_types = FALSE)
    comparison <- results$summary %>%
      left_join(expected_results[, c("population", "status")], by = c("Population" = "population"))
    
    cat("Comparison of results with expected:\\n")
    print(comparison[, c("Population", "N", "P_value", "Significant", "status")])
    
    # Check if results match expectations
    hwe_populations <- comparison$Population[comparison$status == "HWE"]
    non_hwe_populations <- comparison$Population[grepl("NON_HWE", comparison$status)]
    
    cat("\\n=== VALIDATION RESULTS ===\\n")
    cat("Populations expected to be in HWE:", length(hwe_populations), "\\n")
    cat("Populations expected to deviate from HWE:", length(non_hwe_populations), "\\n")
    
    if(length(hwe_populations) > 0) {
      hwe_passed <- sum(comparison$Significant[comparison$status == "HWE"] %in% FALSE, na.rm = TRUE)
      cat("HWE populations correctly identified:", hwe_passed, "/", length(hwe_populations), "\\n")
    }
    
    if(length(non_hwe_populations) > 0) {
      non_hwe_passed <- sum(comparison$Significant[grepl("NON_HWE", comparison$status)] %in% TRUE, na.rm = TRUE)
      cat("Non-HWE populations correctly identified:", non_hwe_passed, "/", length(non_hwe_populations), "\\n")
    }
    
    # Create plots
    overview_plots <- create_population_overview_plots(results$summary, "MockData_Test")
    if(!is.null(overview_plots)) {
      print(overview_plots$sample_sizes)
      print(overview_plots$allele_frequencies)
      print(overview_plots$hwe_pvalues)
      print(overview_plots$summary)
    }
  }
} else {
  cat("ERROR: No data found after merging mock datasets!\\n")
}

cat("\\n=== MOCK DATA TESTING COMPLETE ===\\n")


writeLines(test_script, "mock_data/test_hwe_with_mock_data.R")
cat("✓ Saved: mock_data/test_hwe_with_mock_data.R\n")

cat("\n=== EXPECTED RESULTS SUMMARY ===\n")
cat("The following populations should be in HWE (p > 0.05):\n")
hwe_pops <- mock_populations$population[mock_populations$status == "HWE"]
cat(paste(hwe_pops, collapse = ", "), "\n\n")

cat("The following populations should deviate from HWE (p < 0.05):\n")
non_hwe_pops <- mock_populations$population[grepl("NON_HWE", mock_populations$status)]
cat(paste(non_hwe_pops, collapse = ", "), "\n\n")

cat("Edge cases (may not be testable):\n")
edge_pops <- mock_populations$population[mock_populations$status %in% c("MONOMORPHIC", "SMALL_SAMPLE", "EMPTY")]
cat(paste(edge_pops, collapse = ", "), "\n\n")


