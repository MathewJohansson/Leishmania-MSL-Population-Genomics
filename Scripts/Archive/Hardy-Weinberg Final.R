# Hardy-Weinberg Equilibrium Test for Polysomic Inheritance with ADMIXTURE Data
# Using binary notation: 1 = normal allele, 0 = deletion allele
# For gene copy numbers: 0, 1, 2, 3, 4 (5 genotypes)

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)

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

# Core Hardy-Weinberg test function (binary notation)
polysomic_hwe_test <- function(genotype_counts, population_name = "Population") {
  if(length(genotype_counts) != 5) {
    stop("genotype_counts must have exactly 5 values: 1111, 1110, 1100, 1000, 0000")
  }
  
  names(genotype_counts) <- c("1111", "1110", "1100", "1000", "0000")
  total_n <- sum(genotype_counts)
  
  if(total_n == 0) {
    stop("No individuals found in the data")
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
  
  # More lenient approach for small samples
  if(total_n < 20) {
    valid_genotypes <- rep(TRUE, 5)
    test_reliable <- FALSE
    warning(paste("Population", population_name, ": Small sample size (n =", total_n, "). Results should be interpreted cautiously."))
  } else {
    valid_genotypes <- expected_counts >= 5
    test_reliable <- sum(valid_genotypes) >= 3
  }
  
  # Calculate chi-square using valid genotypes
  if(sum(valid_genotypes) >= 2) {
    chi_square <- sum((genotype_counts[valid_genotypes] - expected_counts[valid_genotypes])^2 / 
                        expected_counts[valid_genotypes])
    df <- sum(valid_genotypes) - 1 - 1
    
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
    warning(paste("Population", population_name, ": Cannot perform chi-square test"))
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

# Function to test all populations from ADMIXTURE data
test_admixture_populations <- function(genotype_counts_df, min_sample_size = 10) {
  cat("=== TESTING POPULATIONS ===\n")
  cat("Input genotype counts table:\n")
  print(genotype_counts_df)
  cat("\n")
  
  # Filter populations with sufficient sample size
  pop_sizes <- rowSums(genotype_counts_df)
  cat("Population sizes before filtering:\n")
  print(pop_sizes)
  
  valid_pops <- names(pop_sizes)[pop_sizes >= min_sample_size]
  
  if(length(valid_pops) == 0) {
    cat("ERROR: No populations meet minimum sample size requirement of", min_sample_size, "\n")
    cat("Consider lowering min_sample_size or checking your data\n")
    cat("Available population sizes:", pop_sizes, "\n")
    return(list(detailed_results = list(), summary = data.frame()))
  }
  
  cat("Testing", length(valid_pops), "populations with >=", min_sample_size, "samples:", paste(valid_pops, collapse = ", "), "\n\n")
  
  all_results <- list()
  summary_table <- data.frame()
  
  for(pop_name in valid_pops) {
    cat("=== Testing", pop_name, "===\n")
    
    result <- polysomic_hwe_test(as.numeric(genotype_counts_df[pop_name, ]), pop_name)
    all_results[[pop_name]] <- result
    
    cat("Sample size:", result$sample_size, "\n")
    cat("Normal allele frequency (1):", round(result$allele_freq_normal, 4), "\n")
    cat("Deletion allele frequency (0):", round(result$allele_freq_deletion, 4), "\n")
    cat("Chi-square:", round(result$chi_square, 4), "\n")
    cat("Degrees of freedom:", result$df, "\n")
    cat("P-value:", format(result$p_value, scientific = TRUE), "\n")
    cat("Significant deviation from HWE:", result$significant, "\n\n")
    
    print(result$results_table)
    cat("\n")
    
    summary_table <- rbind(summary_table, data.frame(
      Population = pop_name,
      N = result$sample_size,
      p_normal = round(result$allele_freq_normal, 4),  # frequency of normal allele (1)
      q_deletion = round(result$allele_freq_deletion, 4),  # frequency of deletion allele (0)
      Chi_square = ifelse(is.na(result$chi_square), NA, round(result$chi_square, 4)),
      df = result$df,
      P_value = result$p_value,
      Significant = result$significant,
      Test_reliable = ifelse(is.null(result$test_reliable), TRUE, result$test_reliable)
    ))
  }
  
  cat("Final summary table:\n")
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

# Visualization functions
plot_hwe_results <- function(hwe_result) {
  df <- hwe_result$results_table
  df$Genotype <- factor(df$Genotype, levels = c("1111", "1110", "1100", "1000", "0000"))
  
  plot_data <- data.frame(
    Genotype = rep(df$Genotype, 2),
    Count = c(df$Observed, df$Expected_Count),
    Type = rep(c("Observed", "Expected"), each = nrow(df))
  )
  
  ggplot(plot_data, aes(x = Genotype, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.7) +
    labs(title = paste("Hardy-Weinberg Test:", hwe_result$population),
         subtitle = paste("N =", hwe_result$sample_size, 
                          ", χ² =", round(hwe_result$chi_square, 3), 
                          ", p =", format(hwe_result$p_value, digits = 3)),
         x = "Genotype (Binary representation)", 
         y = "Count",
         caption = "1111 = 4 copies, 1110 = 3 copies, 1100 = 2 copies, 1000 = 1 copy, 0000 = 0 copies\n1 = normal allele, 0 = deletion allele") +
    theme_minimal() +
    scale_fill_manual(values = c("Observed" = "steelblue", "Expected" = "orange")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

plot_summary_results <- function(summary_table) {
  # Remove rows with NA p-values for plotting
  valid_rows <- !is.na(summary_table$P_value) & is.finite(summary_table$P_value) & summary_table$P_value > 0
  
  if(sum(valid_rows) == 0) {
    cat("Warning: No valid p-values for plotting. This often happens with small sample sizes.\n")
    cat("Showing allele frequencies only.\n")
    
    # Just show allele frequencies
    p1 <- ggplot(summary_table, aes(x = Population, y = p_normal)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      geom_text(aes(label = paste("n =", N)), vjust = -0.5, size = 3) +
      labs(title = "Allele Frequencies Across Populations",
           subtitle = "Normal allele (1) frequency by population",
           x = "Population", y = "Frequency of Normal Allele (1)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim(0, 1.1)
    
    return(list(allele_freq_plot = p1, pvalue_plot = NULL))
  }
  
  # Plot 1: Allele frequencies
  p1 <- ggplot(summary_table, aes(x = Population, y = p_normal)) +
    geom_col(fill = "steelblue", alpha = 0.7) +
    geom_text(aes(label = paste("n =", N)), vjust = -0.5, size = 3) +
    labs(title = "Allele Frequencies Across Populations",
         subtitle = "Normal allele (1) frequency by population", 
         x = "Population", y = "Frequency of Normal Allele (1)") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ylim(0, 1.1)
  
  # Plot 2: P-values (only for valid ones with capped extremes)
  valid_summary <- summary_table[valid_rows, ]
  
  if(nrow(valid_summary) > 0) {
    # Cap extremely small p-values for better visualization
    valid_summary$P_value_capped <- pmax(valid_summary$P_value, 1e-10)
    valid_summary$log10_p <- -log10(valid_summary$P_value_capped)
    
    p2 <- ggplot(valid_summary, aes(x = Population, y = log10_p, 
                                    fill = Significant)) +
      geom_col(alpha = 0.7) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
      geom_hline(yintercept = -log10(0.001), linetype = "dotted", color = "darkred", size = 1) +
      geom_text(aes(label = paste("n =", N)), vjust = -0.5, size = 3) +
      labs(title = "Hardy-Weinberg Test Results",
           subtitle = "Higher bars = more significant deviation from HWE",
           x = "Population", y = "-log10(P-value)",
           caption = "Red dashed line: p = 0.05, Red dotted line: p = 0.001\nP-values < 1e-10 are capped for visualization") +
      scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "orange"),
                        name = "Significant\nDeviation",
                        na.value = "gray") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      coord_cartesian(ylim = c(0, 15))  # Limit y-axis to reasonable range
  } else {
    p2 <- NULL
  }
  
  return(list(allele_freq_plot = p1, pvalue_plot = p2))
}

# ===== MAIN ANALYSIS WITH YOUR ACTUAL DATA =====

cat("=== HARDY-WEINBERG TEST FOR POLYSOMIC INHERITANCE ===\n")
cat("Using binary notation: 1 = normal allele, 0 = deletion allele\n\n")

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

# STEP 2: Process each ADMIXTURE result using your data
cat("STEP 2: Processing ADMIXTURE files\n\n")

# Create a list of population datasets with descriptive names
pop_datasets <- list(
  "Americas_k13" = pop_membership_Americas_k13,
  "Infantum_k13" = pop_membership_Infantum_k13,
  "AllSamples_k12" = pop_membership_All_k12,
  "Infantum_k8" = pop_membership_Infantum_k8
)

# Function to create comprehensive population plots (with error handling)
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
  
  # Sort populations by name for consistent ordering
  results_summary <- results_summary[order(results_summary$Population), ]
  results_summary$Population <- factor(results_summary$Population, levels = results_summary$Population)
  
  # Ensure numeric columns are properly formatted
  results_summary$N <- as.numeric(results_summary$N)
  results_summary$p_normal <- as.numeric(results_summary$p_normal)
  results_summary$P_value <- as.numeric(results_summary$P_value)
  
  plots_list <- list()
  
  # Plot 1: Sample sizes across all populations
  tryCatch({
    p1 <- ggplot(results_summary, aes(x = Population, y = N)) +
      geom_col(fill = "lightblue", alpha = 0.7, color = "black") +
      geom_text(aes(label = N), vjust = -0.3, size = 3) +
      labs(title = paste("Sample Sizes -", dataset_name),
           x = "Population", y = "Sample Size (N)",
           subtitle = paste("Total populations:", nrow(results_summary))) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim(0, max(results_summary$N, na.rm = TRUE) * 1.1)
    plots_list$sample_sizes <- p1
  }, error = function(e) {
    cat("Error creating sample sizes plot:", e$message, "\n")
    plots_list$sample_sizes <- NULL
  })
  
  # Plot 2: Normal allele frequencies across all populations
  tryCatch({
    p2 <- ggplot(results_summary, aes(x = Population, y = p_normal)) +
      geom_col(aes(fill = p_normal), alpha = 0.8) +
      geom_text(aes(label = round(p_normal, 3)), vjust = -0.3, size = 2.5) +
      scale_fill_gradient2(low = "red", mid = "yellow", high = "darkgreen", 
                           midpoint = 0.5, name = "Frequency") +
      labs(title = paste("Normal Allele Frequencies -", dataset_name),
           x = "Population", y = "Frequency of Normal Allele (1)",
           subtitle = "Red = low frequency, Green = high frequency") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ylim(0, 1.1)
    plots_list$allele_frequencies <- p2
  }, error = function(e) {
    cat("Error creating allele frequencies plot:", e$message, "\n")
    plots_list$allele_frequencies <- NULL
  })
  
  # Plot 3: Hardy-Weinberg p-values (only for populations with valid p-values)
  tryCatch({
    # Filter for populations with valid p-values
    valid_pvalue_data <- results_summary[!is.na(results_summary$P_value) & 
                                           is.finite(results_summary$P_value) & 
                                           results_summary$P_value > 0, ]
    
    if(nrow(valid_pvalue_data) > 0) {
      valid_pvalue_data$log10_p <- -log10(pmax(valid_pvalue_data$P_value, 1e-15))
      valid_pvalue_data$HWE_Status <- ifelse(valid_pvalue_data$Significant, 
                                             "Deviates from HWE", "Fits HWE")
      
      p3 <- ggplot(valid_pvalue_data, aes(x = Population, y = log10_p, fill = HWE_Status)) +
        geom_col(alpha = 0.8) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", size = 1) +
        geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "darkred", size = 1) +
        geom_text(aes(label = round(log10_p, 1)), vjust = -0.3, size = 2.5) +
        scale_fill_manual(values = c("Fits HWE" = "steelblue", 
                                     "Deviates from HWE" = "orange"),
                          name = "HWE Status") +
        labs(title = paste("Hardy-Weinberg P-values -", dataset_name),
             x = "Population", y = "-log10(P-value)",
             subtitle = "Only populations with valid p-values shown",
             caption = "Red lines: p=0.05 (dashed), p=0.01 (dotted)") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        coord_cartesian(ylim = c(0, 15))
      plots_list$hwe_pvalues <- p3
    } else {
      cat("No valid p-values to plot for", dataset_name, "\n")
      plots_list$hwe_pvalues <- NULL
    }
  }, error = function(e) {
    cat("Error creating p-values plot:", e$message, "\n")
    plots_list$hwe_pvalues <- NULL
  })
  
  # Plot 4: Simple overview combining sample size and significance
  tryCatch({
    results_summary$Size_Category <- ifelse(results_summary$N < 10, "Small (<10)",
                                            ifelse(results_summary$N < 20, "Medium (10-19)", "Large (≥20)"))
    results_summary$HWE_Status <- ifelse(is.na(results_summary$Significant), "No Test",
                                         ifelse(results_summary$Significant, "Significant", "Not Significant"))
    
    p4 <- ggplot(results_summary, aes(x = Population, y = N, fill = Size_Category)) +
      geom_col(alpha = 0.7) +
      geom_point(aes(color = HWE_Status), y = max(results_summary$N, na.rm = TRUE) * 0.9, size = 3) +
      scale_fill_manual(values = c("Small (<10)" = "lightcoral", 
                                   "Medium (10-19)" = "lightyellow", 
                                   "Large (≥20)" = "lightgreen"),
                        name = "Sample Size") +
      scale_color_manual(values = c("Not Significant" = "blue", 
                                    "Significant" = "red", 
                                    "No Test" = "gray"),
                         name = "HWE Test") +
      labs(title = paste("Population Overview -", dataset_name),
           x = "Population", y = "Sample Size",
           subtitle = "Bars = sample size, Points = HWE test results") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    plots_list$summary <- p4
  }, error = function(e) {
    cat("Error creating summary plot:", e$message, "\n")
    plots_list$summary <- NULL
  })
  
  return(plots_list)
}

# Store results and plots for comparison
all_dataset_results <- list()
all_dataset_plots <- list()  # Add this to store plots

# Process each dataset
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
  
  # Perform Hardy-Weinberg tests (use lower min_sample_size for better results)
  results <- test_admixture_populations(genotype_counts, min_sample_size = 5)
  
  # Store results
  all_dataset_results[[dataset_name]] <- results
  
  cat("\n=== SUMMARY TABLE FOR", dataset_name, "===\n")
  if(nrow(results$summary) > 0) {
    print(results$summary)
    
    # Create comprehensive population overview plots
    cat("Creating comprehensive population plots for", dataset_name, "...\n")
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
    
    # Individual population plots (for significant deviations only to save space)
    tryCatch({
      significant_pops <- results$summary[!is.na(results$summary$P_value) & 
                                            !is.na(results$summary$Significant) & 
                                            results$summary$Significant == TRUE, "Population"]
      if(length(significant_pops) > 0) {
        cat("Creating detailed plots for", length(significant_pops), "populations with significant HWE deviations...\n")
        for(pop_name in significant_pops) {
          if(pop_name %in% names(results$detailed_results)) {
            tryCatch({
              plot_obj <- plot_hwe_results(results$detailed_results[[pop_name]])
              print(plot_obj)
              # Store for saving
              all_dataset_plots[[dataset_name]]$detailed[[pop_name]] <- plot_obj
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
    cat("No populations met minimum sample size requirements\n")
  }
}

cat("\n", rep("=", 60), "\n")
cat("=== ANALYSIS COMPLETE ===\n")

# Final summary across all datasets
cat("\n=== CROSS-DATASET COMPARISON ===\n")
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    significant_count <- sum(results$summary$Significant %in% TRUE, na.rm = TRUE)
    total_count <- nrow(results$summary)
    cat(dataset_name, ": ", significant_count, "/", total_count, " populations show significant HWE deviation\n", sep = "")
  }
}

# ===== SAVE ALL PLOTS TO FILES =====
cat("\n=== SAVING PLOTS TO FILES ===\n")

# Create output directory if it doesn't exist
if(!dir.exists("HWE_plots")) {
  dir.create("HWE_plots")
  cat("Created directory: HWE_plots/\n")
} else {
  cat("Using existing directory: HWE_plots/\n")
}

# Function to safely save plots
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
                   paste0("HWE_plots/", dataset_name, "_sample_sizes.png"))) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$allele_frequencies, 
                   paste0("HWE_plots/", dataset_name, "_allele_frequencies.png"))) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$hwe_pvalues, 
                   paste0("HWE_plots/", dataset_name, "_hwe_pvalues.png"))) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$summary, 
                   paste0("HWE_plots/", dataset_name, "_overview.png"))) {
      saved_count <- saved_count + 1
    }
  }
  
  # Save individual population plots
  if(!is.null(dataset_plots$detailed) && length(dataset_plots$detailed) > 0) {
    cat("Saving", length(dataset_plots$detailed), "detailed population plots...\n")
    for(pop_name in names(dataset_plots$detailed)) {
      if(safe_ggsave(dataset_plots$detailed[[pop_name]], 
                     paste0("HWE_plots/", dataset_name, "_", pop_name, "_detailed.png"),
                     width = 10, height = 6)) {
        saved_count <- saved_count + 1
      }
    }
  }
}

# Save summary tables as CSV files
cat("\n--- Saving summary tables ---\n")
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    tryCatch({
      filename <- paste0("HWE_plots/", dataset_name, "_summary_table.csv")
      write.csv(results$summary, filename, row.names = FALSE)
      cat("✓ Saved:", filename, "\n")
      saved_count <- saved_count + 1
    }, error = function(e) {
      cat("✗ Error saving summary table for", dataset_name, ":", e$message, "\n")
    })
  }
}

# Create a combined summary table across all datasets
cat("\n--- Creating combined summary ---\n")
combined_summary <- data.frame()
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    temp_summary <- results$summary
    temp_summary$Dataset <- dataset_name
    combined_summary <- rbind(combined_summary, temp_summary)
  }
}

if(nrow(combined_summary) > 0) {
  # Save combined summary
  tryCatch({
    write.csv(combined_summary, "HWE_plots/combined_summary_all_datasets.csv", row.names = FALSE)
    cat("✓ Saved: HWE_plots/combined_summary_all_datasets.csv\n")
    saved_count <- saved_count + 1
  }, error = function(e) {
    cat("✗ Error saving combined summary:", e$message, "\n")
  })
  
  # Create a final comparison plot across all datasets
  tryCatch({
    comparison_plot <- ggplot(combined_summary, aes(x = Population, y = p_normal, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.7) +
      labs(title = "Normal Allele Frequencies Across All Datasets",
           x = "Population", y = "Frequency of Normal Allele (1)",
           subtitle = "Comparison of all ADMIXTURE results") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      facet_wrap(~Dataset, scales = "free_x")
    
    if(safe_ggsave(comparison_plot, "HWE_plots/all_datasets_comparison.png", width = 16, height = 10)) {
      saved_count <- saved_count + 1
    }
  }, error = function(e) {
    cat("✗ Error creating comparison plot:", e$message, "\n")
  })
}

cat("\n=== PLOT SAVING COMPLETE ===\n")
cat("Total files saved:", saved_count, "\n")
cat("All files saved to: HWE_plots/ directory\n")

# List all saved files
if(dir.exists("HWE_plots")) {
  plot_files <- list.files("HWE_plots", full.names = FALSE)
  if(length(plot_files) > 0) {
    cat("\nFiles in HWE_plots directory:\n")
    for(file in sort(plot_files)) {
      cat(" -", file, "\n")
    }
  } else {
    cat("\nNo files found in HWE_plots directory!\n")
  }
} else {
  cat("\nHWE_plots directory not found!\n")
}
