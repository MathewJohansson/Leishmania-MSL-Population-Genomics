# Complete Hardy-Weinberg Equilibrium Analysis for Tetraploid Systems



# 1. LOAD LIBRARIES ------------------------------------------------------------

# Load required libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(tibble)
library(viridis)
library(RColorBrewer)




# 2. PROJECT DIRECTORY CONFIG --------------------------------------------------

# Set base directory - ONLY CHANGE THIS LINE to run on different machines
BASE_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

# Define subdirectories
METADATA_DIR <- file.path(BASE_DIR, "Metadata")
DATA_DIR     <- file.path(BASE_DIR, "Data/population_structure")
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/HWE")
TABLES_DIR   <- file.path(BASE_DIR, "Tables/HWE")

# Create output directories if needed
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLES_DIR,  recursive = TRUE, showWarnings = FALSE)


# Set colour-blind-friendly palette
colorblind_palette <- c("#0173B2", "#DE8F05", "#029E73", "#CC78BC", "#CA9161", "#FBAFE4", "#949494", "#ECE133")

cat("=== HARDY-WEINBERG EQUILIBRIUM ANALYSIS FOR TETRAPLOID SYSTEMS ===\n")
cat("Analysis started at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")




# 3. UTILITY FUNCTIONS ---------------------------------------------------------

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

# Function to extract population number for natural sorting
extract_pop_number <- function(pop_names) {
  numbers <- as.numeric(gsub("[^0-9]", "", pop_names))
  numbers[is.na(numbers)] <- 999
  return(numbers)
}

# Function to sort populations naturally (v1, v2, v3, ..., v10, v11, etc.)
sort_populations_naturally <- function(pop_names) {
  if(length(pop_names) == 0) return(pop_names)
  
  pop_df <- data.frame(
    original = pop_names,
    prefix = gsub("[0-9]+$", "", pop_names),
    number = extract_pop_number(pop_names),
    stringsAsFactors = FALSE
  )
  
  pop_df <- pop_df[order(pop_df$prefix, pop_df$number), ]
  return(pop_df$original)
}




# 4. DATA PROCESSING FUNCTIONS -------------------------------------------------

# Function to create genotype count table from merged data
create_genotype_counts <- function(merged_data) {
  merged_data$genotype <- coverage_to_genotype(merged_data$coverage_rounded)
  
  genotype_counts <- merged_data %>%
    group_by(population, genotype) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = genotype, values_from = count, values_fill = 0) %>%
    column_to_rownames("population")
  
  # Ensure all genotype columns exist in correct order
  required_genotypes <- c("0000", "1000", "1100", "1110", "1111")
  for(geno in required_genotypes) {
    if(!geno %in% colnames(genotype_counts)) {
      genotype_counts[[geno]] <- 0
    }
  }
  
  # Reorder columns (1111 first for highest coverage)
  genotype_counts <- genotype_counts[, c("1111", "1110", "1100", "1000", "0000")]
  return(genotype_counts)
}

# Function to process population dataset
process_population_dataset <- function(pop_data, samples_data, dataset_name, min_proportion = 0.5) {
  cat("Processing population data for:", dataset_name, "\n")
  cat("Raw population data:", nrow(pop_data), "entries\n")
  
  # Filter samples with sufficient population assignment
  pop_data_filtered <- pop_data %>%
    filter(population_proportion >= min_proportion)
  
  cat("After filtering by proportion >=", min_proportion, ":", nrow(pop_data_filtered), "entries\n")
  
  # Merge with sample data
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




# 5. CORE HARDY-WEINBERG TESTING FUNCTIONS -------------------------------------

# Fixed Hardy-Weinberg test function for tetraploids
polysomic_hwe_test_fixed <- function(genotype_counts, population_name = "Population") {
  if(length(genotype_counts) != 5) {
    stop("genotype_counts must have exactly 5 values: 1111, 1110, 1100, 1000, 0000")
  }
  
  names(genotype_counts) <- c("1111", "1110", "1100", "1000", "0000")
  total_n <- sum(genotype_counts)
  
  # Handle empty populations
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
      test_performed = FALSE,
      test_type = "No test - empty population",
      results_table = data.frame(),
      testable_classes = 0
    ))
  }
  
  # Count number of genotype classes with observations
  observed_classes <- sum(genotype_counts > 0)
  
  # Handle monomorphic populations
  if(observed_classes <= 1) {
    if(genotype_counts["1111"] > 0) {
      p_est <- 1.0; q_est <- 0.0
    } else if(genotype_counts["0000"] > 0) {
      p_est <- 0.0; q_est <- 1.0
    } else if(genotype_counts["1110"] > 0) {
      p_est <- 0.75; q_est <- 0.25
    } else if(genotype_counts["1100"] > 0) {
      p_est <- 0.5; q_est <- 0.5
    } else if(genotype_counts["1000"] > 0) {
      p_est <- 0.25; q_est <- 0.75
    } else {
      p_est <- NA; q_est <- NA
    }
    
    return(list(
      population = population_name,
      sample_size = total_n,
      allele_freq_normal = p_est,
      allele_freq_deletion = q_est,
      chi_square = NA,
      df = NA,
      p_value = NA,
      significant = NA,
      test_reliable = FALSE,
      test_performed = FALSE,
      test_type = "No test - monomorphic",
      results_table = data.frame(
        Genotype = names(genotype_counts),
        Coverage = c(4, 3, 2, 1, 0),
        Observed = as.numeric(genotype_counts),
        Expected = rep(NA, 5)
      ),
      testable_classes = 1
    ))
  }
  
  # Calculate allele frequencies using maximum likelihood estimation
  total_alleles <- total_n * 4
  normal_alleles <- 4*genotype_counts["1111"] + 3*genotype_counts["1110"] + 
    2*genotype_counts["1100"] + 1*genotype_counts["1000"] + 0*genotype_counts["0000"]
  
  p <- normal_alleles / total_alleles  # frequency of normal allele
  q <- 1 - p  # frequency of deletion allele
  
  # Prevent extreme frequencies that cause numerical issues
  if(p <= 0.001) p <- 0.001
  if(p >= 0.999) p <- 0.999
  q <- 1 - p
  
  # Calculate expected frequencies under HWE for tetraploid
  expected_freq <- c(
    "1111" = p^4,                    # AAAA
    "1110" = 4 * p^3 * q,           # AAAa
    "1100" = 6 * p^2 * q^2,         # AAaa
    "1000" = 4 * p * q^3,           # Aaaa
    "0000" = q^4                     # aaaa
  )
  
  expected_counts <- expected_freq * total_n
  
  # More appropriate criteria for tetraploid HWE testing
  min_expected_count <- 1.0  # More lenient than traditional 5
  min_expected_prop <- 0.01  # At least 1% expected frequency
  
  # Determine which classes can be reliably tested
  valid_for_test <- (expected_counts >= min_expected_count) & 
    (expected_freq >= min_expected_prop)
  
  testable_classes <- sum(valid_for_test)
  
  if(testable_classes < 2) {
    return(list(
      population = population_name,
      sample_size = total_n,
      allele_freq_normal = p,
      allele_freq_deletion = q,
      chi_square = NA,
      df = NA,
      p_value = NA,
      significant = NA,
      test_reliable = FALSE,
      test_performed = FALSE,
      test_type = "No test - insufficient expected counts",
      results_table = data.frame(
        Genotype = names(genotype_counts),
        Coverage = c(4, 3, 2, 1, 0),
        Observed = as.numeric(genotype_counts),
        Expected = expected_counts,
        Valid_for_test = valid_for_test
      ),
      testable_classes = testable_classes
    ))
  }
  
  # Perform chi-square test using only valid classes
  observed_valid <- genotype_counts[valid_for_test]
  expected_valid <- expected_counts[valid_for_test]
  
  # Calculate chi-square statistic
  chi_square <- sum((observed_valid - expected_valid)^2 / expected_valid)
  
  # Correct degrees of freedom calculation for tetraploid HWE
  if(testable_classes == 5) {
    df <- 3  # All 5 classes: df = 5 - 1 (constraint) - 1 (estimated parameter) = 3
  } else if(testable_classes == 4) {
    df <- 2
  } else if(testable_classes == 3) {
    df <- 1
  } else {
    df <- max(1, testable_classes - 2)
  }
  
  # Calculate p-value
  if(df > 0) {
    p_value <- pchisq(chi_square, df, lower.tail = FALSE)
    significant <- (p_value < 0.05)
    test_reliable <- (total_n >= 10) && (testable_classes >= 3)
    test_performed <- TRUE
    test_type <- "Chi-square test"
  } else {
    p_value <- NA
    significant <- NA
    test_reliable <- FALSE
    test_performed <- FALSE
    test_type <- "No test - insufficient degrees of freedom"
  }
  
  # For very small samples, mark as unreliable
  if(total_n < 10) {
    test_reliable <- FALSE
    if(test_performed) {
      test_type <- "Chi-square test (small sample - interpret cautiously)"
    }
  }
  
  # Create detailed results table
  results_table <- data.frame(
    Genotype = names(genotype_counts),
    Coverage = c(4, 3, 2, 1, 0),
    Observed = as.numeric(genotype_counts),
    Expected_Freq = expected_freq,
    Expected_Count = expected_counts,
    Residual = genotype_counts - expected_counts,
    Standardized_Residual = (genotype_counts - expected_counts) / sqrt(expected_counts),
    Valid_for_test = valid_for_test
  )
  
  return(list(
    population = population_name,
    sample_size = total_n,
    allele_freq_normal = p,
    allele_freq_deletion = q,
    chi_square = chi_square,
    df = df,
    p_value = p_value,
    significant = significant,
    test_reliable = test_reliable,
    test_performed = test_performed,
    test_type = test_type,
    results_table = results_table,
    testable_classes = testable_classes
  ))
}

# Function to test all populations
test_admixture_populations_fixed <- function(genotype_counts_df, min_sample_size = 0) {
  cat("=== TESTING ALL POPULATIONS ===\n")
  cat("Input genotype counts table:\n")
  print(genotype_counts_df)
  cat("\n")
  
  all_populations <- rownames(genotype_counts_df)
  pop_sizes <- rowSums(genotype_counts_df)
  
  cat("Population sizes:\n")
  print(pop_sizes)
  cat("\n")
  
  sorted_populations <- sort_populations_naturally(all_populations)
  cat("Testing populations in order:", paste(sorted_populations, collapse = ", "), "\n\n")
  
  all_results <- list()
  summary_table <- data.frame()
  
  for(pop_name in sorted_populations) {
    cat("=== Testing", pop_name, "===\n")
    
    pop_counts <- as.numeric(genotype_counts_df[pop_name, ])
    names(pop_counts) <- colnames(genotype_counts_df)
    
    result <- polysomic_hwe_test_fixed(pop_counts, pop_name)
    all_results[[pop_name]] <- result
    
    cat("Sample size:", result$sample_size, "\n")
    if(result$sample_size > 0) {
      cat("Normal allele frequency:", round(result$allele_freq_normal, 4), "\n")
      cat("Deletion allele frequency:", round(result$allele_freq_deletion, 4), "\n")
      cat("Test type:", result$test_type, "\n")
      
      if(result$test_performed && !is.na(result$chi_square)) {
        cat("Chi-square:", round(result$chi_square, 4), "\n")
        cat("Degrees of freedom:", result$df, "\n")
        cat("Testable classes:", result$testable_classes, "\n")
        cat("P-value:", format(result$p_value, scientific = TRUE), "\n")
        cat("Significant deviation from HWE:", result$significant, "\n")
        cat("Test reliable:", result$test_reliable, "\n")
      } else {
        cat("HWE test: Not performed -", result$test_type, "\n")
      }
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
      Test_performed = result$test_performed,
      Test_type = result$test_type,
      Testable_classes = result$testable_classes
    ))
  }
  
  # Sort summary table naturally
  summary_table$Population <- factor(summary_table$Population, levels = sorted_populations)
  summary_table <- summary_table[order(summary_table$Population), ]
  
  cat("=== FINAL SUMMARY TABLE ===\n")
  print(summary_table)
  cat("\n")
  
  return(list(detailed_results = all_results, summary = summary_table))
}




# 6. VISUALISATION FUNCTIONS ------------------------------------------------------

# Function to create individual HWE plots
plot_hwe_results_fixed <- function(hwe_result) {
  if(is.null(hwe_result$results_table) || nrow(hwe_result$results_table) == 0 || hwe_result$sample_size == 0) {
    return(NULL)
  }
  
  df <- hwe_result$results_table
  df$Genotype <- factor(df$Genotype, levels = c("1111", "1110", "1100", "1000", "0000"))
  
  plot_data <- data.frame(
    Genotype = rep(df$Genotype, 2),
    Count = c(df$Observed, df$Expected_Count),
    Type = rep(c("Observed", "Expected"), each = nrow(df))
  )
  
  # Create subtitle with test information
  if(hwe_result$test_performed && !is.na(hwe_result$chi_square)) {
    subtitle <- paste("N =", hwe_result$sample_size, 
                      ", χ² =", round(hwe_result$chi_square, 3),
                      ", df =", hwe_result$df,
                      ", p =", format(hwe_result$p_value, digits = 3),
                      "\nTest:", hwe_result$test_type)
  } else {
    subtitle <- paste("N =", hwe_result$sample_size, "\nTest:", hwe_result$test_type)
  }
  
  ggplot(plot_data, aes(x = Genotype, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    labs(title = paste("Hardy-Weinberg Test:", hwe_result$population),
         subtitle = subtitle,
         x = "Genotype (Binary representation)", 
         y = "Count",
         caption = "1111 = 4 copies, 1110 = 3 copies, 1100 = 2 copies, 1000 = 1 copy, 0000 = 0 copies") +
    theme_minimal() +
    scale_fill_manual(values = c("Observed" = colorblind_palette[1], "Expected" = colorblind_palette[2])) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 12, face = "bold"),
          plot.subtitle = element_text(size = 9))
}

# Function to create comprehensive overview plots
create_population_overview_plots_fixed <- function(results_summary, dataset_name) {
  if(nrow(results_summary) == 0) {
    cat("No data to plot for", dataset_name, "\n")
    return(NULL)
  }
  
  results_summary <- results_summary[!is.na(results_summary$Population), ]
  if(nrow(results_summary) == 0) {
    cat("No valid population data for", dataset_name, "\n")
    return(NULL)
  }
  
  # Sort populations naturally
  sorted_pops <- sort_populations_naturally(as.character(results_summary$Population))
  results_summary$Population <- factor(results_summary$Population, levels = sorted_pops)
  results_summary <- results_summary[order(results_summary$Population), ]
  
  # Apply Benjamini-Hochberg correction across all tested populations
  testable_idx <- which(results_summary$Test_performed & !is.na(results_summary$P_value))
  
  if (length(testable_idx) > 0) {
    p_raw <- results_summary$P_value[testable_idx]
    p_adj <- p.adjust(p_raw, method = "BH")
    results_summary$P_value_BH <- NA_real_
    results_summary$P_value_BH[testable_idx] <- p_adj
    results_summary$Significant <- FALSE
    results_summary$Significant[testable_idx] <- p_adj < 0.05
    
    # BH threshold line: largest raw p-value among those surviving BH correction
    surviving_raw <- p_raw[p_adj < 0.05]
    bh_line <- if (length(surviving_raw) > 0) -log10(max(surviving_raw)) else -log10(0.05)
  } else {
    results_summary$P_value_BH <- NA_real_
    bh_line <- -log10(0.05)
  }
  
  # Create HWE status categories
  results_summary$HWE_Status <- ifelse(results_summary$N == 0, "Empty Population",
                                       ifelse(!results_summary$Test_performed, "No Test Possible",
                                              ifelse(is.na(results_summary$Significant), "Test Failed",
                                                     ifelse(results_summary$Significant, "Deviates from HWE", "Fits HWE"))))
  
  # Create test quality categories
  results_summary$Test_Quality <- ifelse(results_summary$N == 0, "Empty",
                                         ifelse(!results_summary$Test_performed, "No Test",
                                                ifelse(!results_summary$Test_reliable, "Unreliable",
                                                       ifelse(results_summary$Testable_classes >= 4, "High Quality",
                                                              ifelse(results_summary$Testable_classes >= 3, "Good Quality", "Low Quality")))))
  
  plots_list <- list()
  
  # Plot 1: Sample sizes
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
  
  # Plot 2: Allele frequencies with quality indicators
  tryCatch({
    valid_freq_data <- results_summary[!is.na(results_summary$p_normal), ]
    if(nrow(valid_freq_data) > 0) {
      p2 <- ggplot(valid_freq_data, aes(x = Population, y = p_normal, fill = Test_Quality)) +
        geom_col(alpha = 0.8, color = "black", size = 0.3) +
        geom_text(aes(label = round(p_normal, 3)), vjust = -0.3, size = 2.5) +
        scale_fill_manual(values = c("Empty" = colorblind_palette[7],
                                     "No Test" = colorblind_palette[7],
                                     "Unreliable" = colorblind_palette[8],
                                     "Low Quality" = colorblind_palette[6],
                                     "Good Quality" = colorblind_palette[4],
                                     "High Quality" = colorblind_palette[3]),
                          name = "Test Quality") +
        labs(title = paste("Normal Allele Frequencies -", dataset_name),
             x = "Population", y = "Frequency of Normal Allele",
             subtitle = "Colour indicates test quality") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 14, face = "bold")) +
        ylim(0, 1.1)
      plots_list$allele_frequencies <- p2
    } else {
      plots_list$allele_frequencies <- NULL
    }
  }, error = function(e) {
    cat("Error creating allele frequencies plot:", e$message, "\n")
    plots_list$allele_frequencies <- NULL
  })
  
  # Plot 3: HWE test results
  tryCatch({
    # Count categories for subtitle
    status_counts <- table(results_summary$HWE_Status)
    status_text <- paste(names(status_counts), "=", status_counts, collapse = ", ")
    
    # Calculate log10 p-values
    results_summary$log10_p <- ifelse(results_summary$Test_performed & !is.na(results_summary$P_value), 
                                      -log10(pmax(results_summary$P_value, 1e-15)), 
                                      0)
    
    p3 <- ggplot(results_summary, aes(x = Population, y = log10_p, fill = HWE_Status)) +
      geom_col(alpha = 0.8, color = "black", size = 0.3) +
      geom_hline(yintercept = bh_line, linetype = "dashed", color = "black", size = 1) +
      geom_text(aes(label = ifelse(Test_performed & !is.na(P_value), round(log10_p, 1), "")), 
                vjust = -0.3, size = 2.5) +
      scale_fill_manual(values = c("Fits HWE" = colorblind_palette[3], 
                                   "Deviates from HWE" = colorblind_palette[2],
                                   "No Test Possible" = colorblind_palette[7],
                                   "Test Failed" = colorblind_palette[8],
                                   "Empty Population" = colorblind_palette[7]),
                        name = "HWE Status", drop = FALSE) +
      labs(title = paste("Hardy-Weinberg Test Results -", dataset_name),
           x = "Population", y = "-log10(P-value)",
           subtitle = paste("Status counts:", status_text),
           caption = "Black dashed line: BH-corrected significance threshold (FDR < 0.05)") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold")) +
      coord_cartesian(ylim = c(0, max(15, max(results_summary$log10_p, na.rm = TRUE) * 1.1)))
    plots_list$hwe_pvalues <- p3
  }, error = function(e) {
    cat("Error creating HWE results plot:", e$message, "\n")
    plots_list$hwe_pvalues <- NULL
  })
  
  # Plot 4: Combined overview
  tryCatch({
    # Sample size categories
    results_summary$Size_Category <- ifelse(results_summary$N == 0, "Empty",
                                            ifelse(results_summary$N < 5, "Very Small (<5)",
                                                   ifelse(results_summary$N < 10, "Small (5-9)",
                                                          ifelse(results_summary$N < 20, "Medium (10-19)", "Large (≥20)"))))
    
    # Calculate triangle position
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
                        name = "Sample Size", drop = FALSE) +
      scale_color_manual(values = c("Fits HWE" = colorblind_palette[1], 
                                    "Deviates from HWE" = colorblind_palette[2], 
                                    "No Test Possible" = colorblind_palette[7],
                                    "Test Failed" = colorblind_palette[8],
                                    "Empty Population" = colorblind_palette[7]),
                         name = "HWE Test", drop = FALSE) +
      labs(title = paste("Population Overview -", dataset_name),
           x = "Population", y = "Sample Size",
           subtitle = paste("Bars = sample size, Triangles = HWE test results\n", status_text)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 14, face = "bold"))
    plots_list$summary <- p4
  }, error = function(e) {
    cat("Error creating overview plot:", e$message, "\n")
    plots_list$summary <- NULL
  })
  
  return(plots_list)
}

# Function to safely save plots
safe_ggsave <- function(plot_obj, filename, width = 12, height = 8) {
  if(!is.null(plot_obj)) {
    # Save as both PNG and JPG - derive JPG filename from PNG filename
    jpg_filename <- sub("\\.png$", ".jpg", filename)
    tryCatch({
      ggsave(filename,     plot = plot_obj, width = width, height = height, dpi = 300)
      ggsave(jpg_filename, plot = plot_obj, width = width, height = height, dpi = 300)
      cat("✓ Saved:", filename, "\n")
      cat("✓ Saved:", jpg_filename, "\n")
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




# 7. MAIN ANALYSIS EXECUTION ---------------------------------------------------

cat("STEP 1: Loading data files\n")

# Load sample data with coverage information
samples_data <- read_tsv(file.path(METADATA_DIR, "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"), show_col_types = FALSE)
cat("Loaded samples data:", nrow(samples_data), "samples\n")

# Load population membership tables
pop_membership_Americas_k13 <- read_tsv(file.path(DATA_DIR, "population_membership_table_Americas_only_k13_AM_notation.tsv"), show_col_types = FALSE)
pop_membership_Infantum_k13 <- read_tsv(file.path(DATA_DIR, "population_membership_table_infantum_only_k13_LI_notation.tsv"), show_col_types = FALSE)
pop_membership_All_k12 <- read_tsv(file.path(DATA_DIR, "population_membership_table_allsamples_k12_LO_notation.tsv"), show_col_types = FALSE)
pop_membership_Infantum_k8 <- read_tsv(file.path(DATA_DIR, "population_membership_table_infantum_only_k8_LI_notation.tsv"), show_col_types = FALSE)

cat("Loaded all population membership files\n\n")

# Check the structure of samples data
cat("=== SAMPLES DATA STRUCTURE ===\n")
cat("Sample data dimensions:", nrow(samples_data), "rows x", ncol(samples_data), "columns\n")
cat("Coverage data distribution:\n")
print(table(samples_data$coverage_rounded, useNA = "ifany"))
cat("\n")

cat("STEP 2: Processing ADMIXTURE datasets\n\n")

# Create list of population datasets
pop_datasets <- list(
  "Americas_k13" = pop_membership_Americas_k13,
  "Infantum_k13" = pop_membership_Infantum_k13,
  "AllSamples_k12" = pop_membership_All_k12,
  "Infantum_k8" = pop_membership_Infantum_k8
)

# Store results and plots
all_dataset_results <- list()
all_dataset_plots <- list()

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
  cat("\n")
  
  # Perform Hardy-Weinberg tests
  results <- test_admixture_populations_fixed(genotype_counts, min_sample_size = 0)
  
  # Store results
  all_dataset_results[[dataset_name]] <- results
  
  # Apply BH correction immediately so all downstream plots use corrected values
  summ <- all_dataset_results[[dataset_name]]$summary
  testable_idx <- which(summ$Test_performed & !is.na(summ$P_value))
  if(length(testable_idx) > 0) {
    p_adj <- p.adjust(summ$P_value[testable_idx], method = "BH")
    summ$Significant[testable_idx] <- p_adj < 0.05
    all_dataset_results[[dataset_name]]$summary <- summ
    results$summary <- summ
  }
  
  cat("\n=== SUMMARY TABLE FOR", dataset_name, "===\n")
  if(nrow(results$summary) > 0) {
    print(results$summary)
    
    # Create comprehensive population overview plots
    cat("Creating population overview plots for", dataset_name, "...\n")
    overview_plots <- create_population_overview_plots_fixed(results$summary, dataset_name)
    
    # Store plots for saving later
    all_dataset_plots[[dataset_name]] <- list()
    all_dataset_plots[[dataset_name]]$overview <- overview_plots
    all_dataset_plots[[dataset_name]]$detailed <- list()
    
    if(!is.null(overview_plots)) {
      # Display overview plots
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
    
    # Create individual population plots for populations with significant deviations
    populations_with_data <- results$summary[results$summary$N > 0, "Population"]
    significant_pops <- results$summary[!is.na(results$summary$P_value) & 
                                          !is.na(results$summary$Significant) & 
                                          results$summary$Significant == TRUE, "Population"]
    
    if(length(significant_pops) > 0) {
      cat("Creating detailed plots for", length(significant_pops), "populations with significant HWE deviations...\n")
      for(pop_name in significant_pops) {
        if(pop_name %in% names(results$detailed_results)) {
          tryCatch({
            plot_obj <- plot_hwe_results_fixed(results$detailed_results[[pop_name]])
            if(!is.null(plot_obj)) {
              print(plot_obj)
              all_dataset_plots[[dataset_name]]$detailed[[pop_name]] <- plot_obj
            }
          }, error = function(e) {
            cat("Error creating plot for population", pop_name, ":", e$message, "\n")
          })
        }
      }
    }
    
    # Also create plots for a few representative populations with good data
    if(length(populations_with_data) > length(significant_pops)) {
      other_pops <- setdiff(populations_with_data, significant_pops)
      # Select up to 3 representative populations
      sample_pops <- head(other_pops[order(results$summary[results$summary$Population %in% other_pops, "N"], decreasing = TRUE)], 3)
      
      if(length(sample_pops) > 0) {
        cat("Creating detailed plots for", length(sample_pops), "representative populations...\n")
        for(pop_name in sample_pops) {
          if(pop_name %in% names(results$detailed_results)) {
            tryCatch({
              plot_obj <- plot_hwe_results_fixed(results$detailed_results[[pop_name]])
              if(!is.null(plot_obj)) {
                print(plot_obj)
                all_dataset_plots[[dataset_name]]$detailed[[pop_name]] <- plot_obj
              }
            }, error = function(e) {
              cat("Error creating plot for population", pop_name, ":", e$message, "\n")
            })
          }
        }
      }
    }
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
    tested_populations <- sum(results$summary$Test_performed, na.rm = TRUE)
    reliable_tests <- sum(results$summary$Test_reliable, na.rm = TRUE)
    significant_count <- sum(results$summary$Significant %in% TRUE, na.rm = TRUE)
    
    cat(dataset_name, ":\n")
    cat("  Total populations: ", total_populations, "\n")
    cat("  Populations with data: ", populations_with_data, "\n")
    cat("  Populations tested: ", tested_populations, "\n")
    cat("  Reliable tests: ", reliable_tests, "\n")
    cat("  Significant HWE deviations: ", significant_count, "\n\n")
  }
}




# 8. SAVE RESULTS AND PLOTS ----------------------------------------------------

cat("SAVING RESULTS")

# Create output directory
# Output directory already created at top of script (FIGURES_DIR)

saved_count <- 0

# Save plots for each dataset
for(dataset_name in names(all_dataset_plots)) {
  cat("\n--- Saving plots for", dataset_name, "---\n")
  
  dataset_plots <- all_dataset_plots[[dataset_name]]
  
  # Save overview plots
  if(!is.null(dataset_plots$overview)) {
    overview <- dataset_plots$overview
    
    if(safe_ggsave(overview$sample_sizes, 
                   file.path(FIGURES_DIR, paste0(dataset_name, "_sample_sizes.png")), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$allele_frequencies, 
                   file.path(FIGURES_DIR, paste0(dataset_name, "_allele_frequencies.png")), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$hwe_pvalues, 
                   file.path(FIGURES_DIR, paste0(dataset_name, "_hwe_pvalues.png")), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
    
    if(safe_ggsave(overview$summary, 
                   file.path(FIGURES_DIR, paste0(dataset_name, "_overview.png")), 
                   width = 14, height = 8)) {
      saved_count <- saved_count + 1
    }
  }
  
  # Save individual population plots
  if(!is.null(dataset_plots$detailed) && length(dataset_plots$detailed) > 0) {
    plot_names <- names(dataset_plots$detailed)
    plot_names <- sort_populations_naturally(plot_names)
    
    cat("Saving detailed plots for", length(plot_names), "populations...\n")
    for(pop_name in plot_names) {
      if(safe_ggsave(dataset_plots$detailed[[pop_name]], 
                     file.path(FIGURES_DIR, paste0(dataset_name, "_", pop_name, "_detailed.png")),
                     width = 10, height = 6)) {
        saved_count <- saved_count + 1
      }
    }
  }
}

# Save summary tables
cat("\n--- Saving summary tables ---\n")
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    tryCatch({
      # Create clean summary table
      clean_summary <- results$summary
      clean_summary$P_value_formatted <- ifelse(is.na(clean_summary$P_value), "NA",
                                                ifelse(clean_summary$P_value < 0.001, 
                                                       format(clean_summary$P_value, scientific = TRUE, digits = 3),
                                                       format(round(clean_summary$P_value, 4), nsmall = 4)))
      
      filename <- file.path(TABLES_DIR, paste0(dataset_name, "_summary_table.csv"))
      write.csv(clean_summary, filename, row.names = FALSE)
      cat("✓ Saved:", filename, "\n")
      saved_count <- saved_count + 1
    }, error = function(e) {
      cat("✗ Error saving summary table for", dataset_name, ":", e$message, "\n")
    })
  }
}

# Create combined summary
cat("\n--- Creating comprehensive combined summary ---\n")
combined_summary <- data.frame()
for(dataset_name in names(all_dataset_results)) {
  results <- all_dataset_results[[dataset_name]]
  if(nrow(results$summary) > 0) {
    temp_summary <- results$summary
    temp_summary$Dataset <- dataset_name
    temp_summary <- temp_summary[, c("Dataset", setdiff(names(temp_summary), "Dataset"))]
    combined_summary <- rbind(combined_summary, temp_summary)
  }
}

if(nrow(combined_summary) > 0) {
  # Sort combined summary
  combined_summary$Population <- as.character(combined_summary$Population)
  combined_summary <- combined_summary[order(combined_summary$Dataset, 
                                             extract_pop_number(combined_summary$Population)), ]
  
  # Save combined summary
  tryCatch({
    write.csv(combined_summary, file.path(TABLES_DIR, "combined_summary_all_datasets.csv"), row.names = FALSE)
    cat("✓ Saved: combined_summary_all_datasets.csv\n")
    saved_count <- saved_count + 1
  }, error = function(e) {
    cat("✗ Error saving combined summary:", e$message, "\n")
  })
  
  # Create cross-dataset comparison plots
  tryCatch({
    # Comparison plot 1: Allele frequencies
    comparison_plot1 <- ggplot(combined_summary[combined_summary$N > 0 & !is.na(combined_summary$p_normal), ], 
                               aes(x = Population, y = p_normal, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.2) +
      scale_fill_manual(values = colorblind_palette[1:length(unique(combined_summary$Dataset))]) +
      labs(title = "Normal Allele Frequencies Across All Datasets",
           x = "Population", y = "Frequency of Normal Allele",
           subtitle = "Comparison of all ADMIXTURE results") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 16, face = "bold")) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 2)
    
    if(safe_ggsave(comparison_plot1, file.path(FIGURES_DIR, "all_datasets_allele_frequencies_comparison.png"), 
                   width = 18, height = 12)) {
      saved_count <- saved_count + 1
    }
    
    # Comparison plot 2: Sample sizes across datasets
    comparison_plot_sizes <- ggplot(combined_summary, aes(x = Population, y = N, fill = Dataset)) +
      geom_col(position = "dodge", alpha = 0.8, color = "black", size = 0.2) +
      scale_fill_manual(values = colorblind_palette[1:length(unique(combined_summary$Dataset))]) +
      labs(title = "Sample Sizes Across All Datasets",
           x = "Population", y = "Sample Size (N)",
           subtitle = "Comparison of population sizes in all ADMIXTURE results") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(size = 16, face = "bold")) +
      facet_wrap(~Dataset, scales = "free_x", ncol = 2)
    
    if(safe_ggsave(comparison_plot_sizes, file.path(FIGURES_DIR, "all_datasets_sample_sizes_comparison.png"),
                   width = 18, height = 12)) {
      saved_count <- saved_count + 1
    }
    
    # Comparison plot 3: HWE test results
    hwe_summary <- combined_summary[combined_summary$Test_performed & !is.na(combined_summary$P_value) & combined_summary$N > 0, ]
    if(nrow(hwe_summary) > 0) {
      hwe_summary$log10_p <- -log10(pmax(hwe_summary$P_value, 1e-15))
      p_adj_combined <- rep(NA_real_, nrow(hwe_summary))
      testable <- which(!is.na(hwe_summary$P_value))
      p_adj_combined[testable] <- p.adjust(hwe_summary$P_value[testable], method = "BH")
      hwe_summary$HWE_Status <- ifelse(p_adj_combined < 0.05, "Deviates from HWE", "Fits HWE")
      
      comparison_plot2 <- ggplot(hwe_summary, aes(x = Population, y = log10_p, fill = HWE_Status)) +
        geom_col(alpha = 0.8, color = "black", size = 0.2) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 1) +
        scale_fill_manual(values = c("Fits HWE" = colorblind_palette[3], 
                                     "Deviates from HWE" = colorblind_palette[2])) +
        labs(title = "Hardy-Weinberg Test Results Across All Datasets",
             x = "Population", y = "-log10(P-value)",
             subtitle = "Only populations with valid test results shown",
             caption = "Black dashed line: p = 0.05 reference threshold. Bar colours reflect BH-corrected significance (FDR < 0.05).") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(size = 16, face = "bold")) +
        facet_wrap(~Dataset, scales = "free_x", ncol = 2)
      
      if(safe_ggsave(comparison_plot2, file.path(FIGURES_DIR, "all_datasets_hwe_results_comparison.png"), 
                     width = 18, height = 12)) {
        saved_count <- saved_count + 1
      }
    }
  }, error = function(e) {
    cat("✗ Error creating comparison plots:", e$message, "\n")
  })
}

# Create final statistics summary
cat("\n--- Creating final statistics summary ---\n")
tryCatch({
  stats_summary <- data.frame()
  for(dataset_name in names(all_dataset_results)) {
    results <- all_dataset_results[[dataset_name]]
    if(nrow(results$summary) > 0) {
      total_pops <- nrow(results$summary)
      pops_with_data <- sum(results$summary$N > 0, na.rm = TRUE)
      pops_tested <- sum(results$summary$Test_performed, na.rm = TRUE)
      reliable_tests <- sum(results$summary$Test_reliable, na.rm = TRUE)
      pops_significant <- sum(results$summary$Significant %in% TRUE, na.rm = TRUE)
      avg_sample_size <- mean(results$summary$N[results$summary$N > 0], na.rm = TRUE)
      avg_normal_freq <- mean(results$summary$p_normal[results$summary$N > 0], na.rm = TRUE)
      
      stats_summary <- rbind(stats_summary, data.frame(
        Dataset = dataset_name,
        Total_Populations = total_pops,
        Populations_With_Data = pops_with_data,
        Populations_Tested = pops_tested,
        Reliable_Tests = reliable_tests,
        Significant_Deviations = pops_significant,
        Percent_Significant = round(100 * pops_significant / pmax(pops_tested, 1), 1),
        Avg_Sample_Size = round(avg_sample_size, 1),
        Avg_Normal_Allele_Freq = round(avg_normal_freq, 3)
      ))
    }
  }
  
  if(nrow(stats_summary) > 0) {
    write.csv(stats_summary, file.path(TABLES_DIR, "dataset_statistics_summary.csv"), row.names = FALSE)
    cat("✓ Saved: dataset_statistics_summary.csv\n")
    saved_count <- saved_count + 1
    
    cat("\n=== FINAL DATASET STATISTICS ===\n")
    print(stats_summary)
  }
}, error = function(e) {
  cat("✗ Error creating statistics summary:", e$message, "\n")
})

# Save genotype counts (observed and expected) for supplementary table
tryCatch({
  genotype_rows <- data.frame()
  genotypes <- c("1111", "1110", "1100", "1000", "0000")
  for(dataset_name in names(all_dataset_results)) {
    results <- all_dataset_results[[dataset_name]]
    for(pop_name in names(results$detailed_results)) {
      pop_result <- results$detailed_results[[pop_name]]
      if(!is.null(pop_result$results_table) && nrow(pop_result$results_table) > 0) {
        rt <- pop_result$results_table
        row_data <- data.frame(
          Dataset    = dataset_name,
          Population = pop_name,
          N          = pop_result$sample_size,
          stringsAsFactors = FALSE
        )
        for(g in genotypes) {
          g_row <- rt[rt$Genotype == g, ]
          row_data[[paste0("Obs_", g)]] <- if(nrow(g_row) > 0) as.numeric(g_row$Observed) else NA
          row_data[[paste0("Exp_", g)]] <- if(nrow(g_row) > 0 && !is.null(g_row$Expected_Count)) round(as.numeric(g_row$Expected_Count), 2) else NA
        }
        genotype_rows <- rbind(genotype_rows, row_data)
      }
    }
  }
  if(nrow(genotype_rows) > 0) {
    genotype_rows <- genotype_rows[order(genotype_rows$Dataset,
                                         extract_pop_number(genotype_rows$Population)), ]
    write.csv(genotype_rows, file.path(TABLES_DIR, "genotype_counts_observed_expected.csv"), row.names = FALSE)
    cat("✓ Saved: genotype_counts_observed_expected.csv\n")
    saved_count <- saved_count + 1
  }
}, error = function(e) {
  cat("✗ Error saving genotype counts:", e$message, "\n")
})

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Total files saved:", saved_count, "\n")
cat("All files saved to:", FIGURES_DIR, "\n")

# List all saved files
if(dir.exists(FIGURES_DIR)) {
  plot_files <- list.files(FIGURES_DIR, full.names = FALSE)
  if(length(plot_files) > 0) {
    cat("\nFiles in output directory:\n")
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
  }
}

