# Per-population Tajima's D boxplots across all chromosomes
# One plot per AM population (V1-V13 = AM1-AM13)
# Goal: identify which populations show chr31-specific signal in MSL0 strains
# Threshold: admixture proportion >= 0.9




# 1. SETUP ------------------------------------------------------------------------

library(tidyverse)
library(ggplot2)

base_dir <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
data_dir <- file.path(base_dir, "Population Diversity and Differentiation Analysis/Per_population")
output_dir <- file.path(base_dir, "Figures/Population Diversity and Differentiation/Per_population")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Paul Tol palette
# MSL0 = red (deletion), MSL4 = blue (wild-type) — Paul Tol palette
msl_colours <- c("MSL0" = "#BB5566", "MSL4" = "#4477AA")
chr31_highlight <- "#FFDD44"

# Chromosome factor levels
chrom_levels <- c(sprintf("%02d", 1:30), "31", "32", "33", "34", "35", "36")

# q_del values at >= 90% threshold (AM K13) for each population
# V-notation = AM-notation confirmed (V5 = AM5 etc.)
q_del <- c(
  V1  = 1.0000,
  V2  = 0.3533,
  V3  = 1.0000,
  V4  = 0.0417,
  V5  = 0.8036,
  V6  = 1.0000,
  V7  = 1.0000,
  V8  = 0.0500,
  V9  = 1.0000,
  V10 = 1.0000,
  V11 = 1.0000,
  V12 = 1.0000,
  V13 = 0.0000
)




# 2. HELPER FUNCTIONS -------------------------------------------------------------

load_tajd <- function(path, msl_status) {
  read_tsv(path, show_col_types = FALSE) %>%
    filter(N_SNPS >= 10, !is.na(CHROM)) %>%
    mutate(MSL_status = msl_status)
}

clean_chroms <- function(df) {
  df %>%
    mutate(
      chrom_clean = str_replace(CHROM, "LinJ\\.", ""),
      chrom_clean = factor(chrom_clean, levels = chrom_levels)
    ) %>%
    filter(!is.na(chrom_clean))
}

make_plot <- function(data, pop_id, has_msl0, has_msl4) {
  
  q <- q_del[pop_id]
  
  # Build subtitle based on what data is available
  if (has_msl0 && has_msl4) {
    subtitle <- paste0("Chr31 highlighted | q_del = ", round(q, 3),
                       " | MSL0 and MSL4 compared")
  } else if (has_msl0) {
    subtitle <- paste0("Chr31 highlighted | q_del = ", round(q, 3),
                       " | MSL0 only (no MSL4 strains at \u226590% threshold)")
  } else {
    subtitle <- paste0("Chr31 highlighted | q_del = ", round(q, 3),
                       " | MSL4 only (no MSL0 strains at \u226590% threshold)")
  }
  
  p <- ggplot(data, aes(x = chrom_clean, y = TajimaD, fill = MSL_status)) +
    
    # Chr31 highlight
    annotate("rect",
             xmin = which(chrom_levels == "31") - 0.5,
             xmax = which(chrom_levels == "31") + 0.5,
             ymin = -Inf, ymax = Inf,
             fill = chr31_highlight, alpha = 0.3) +
    
    geom_boxplot(
      outlier.size = 0.5,
      outlier.alpha = 0.4,
      linewidth = 0.4,
      position = position_dodge(width = 0.8),
      width = 0.7
    ) +
    
    scale_fill_manual(
      values = msl_colours,
      labels = c("MSL0" = "MSL0 (deletion)", "MSL4" = "MSL4 (wild-type)"),
      name = "MSL status"
    ) +
    
    scale_x_discrete(drop = FALSE) +
    
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    
    labs(
      title = paste0("Tajima's D across all chromosomes — ", pop_id, " (AM K=13)"),
      subtitle = subtitle,
      x = "Chromosome",
      y = "Tajima's D (50kb windows)",
      caption = "Windows with <10 SNPs excluded. Chr31: diploidised tetraploid calls; other chromosomes: diploid calls."
    ) +
    
    theme_bw(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "top",
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(size = 7, colour = "grey50")
    )
  
  return(p)
}




# 3. MAIN LOOP — one plot per population ------------------------------------------

populations <- paste0("V", 1:13)

for (pop in populations) {
  
  cat("Processing", pop, "...\n")
  
  # File paths
  f_msl0_chr31  <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL0_samples.diploidised.tsv.Tajima.D"))
  f_msl0_genome <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL0_samples.tsv.Tajima.D"))
  f_msl4_chr31  <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL4_samples.diploidised.tsv.Tajima.D"))
  f_msl4_genome <- file.path(data_dir, paste0(pop, ".min_proportion0.9_MSL4_samples.tsv.Tajima.D"))
  
  has_msl0 <- file.exists(f_msl0_chr31) && file.exists(f_msl0_genome)
  has_msl4 <- file.exists(f_msl4_chr31) && file.exists(f_msl4_genome)
  
  if (!has_msl0 && !has_msl4) {
    cat("  Skipping", pop, "— no files found\n")
    next
  }
  
  # Load available data
  all_data <- list()
  
  if (has_msl0) {
    msl0 <- bind_rows(
      load_tajd(f_msl0_chr31, "MSL0"),
      load_tajd(f_msl0_genome, "MSL0")
    )
    all_data[["MSL0"]] <- msl0
  }
  
  if (has_msl4) {
    msl4 <- bind_rows(
      load_tajd(f_msl4_chr31, "MSL4"),
      load_tajd(f_msl4_genome, "MSL4")
    )
    all_data[["MSL4"]] <- msl4
  }
  
  combined <- bind_rows(all_data) %>%
    clean_chroms() %>%
    mutate(MSL_status = factor(MSL_status, levels = c("MSL0", "MSL4")))
  
  # Make and save plot
  p <- make_plot(combined, pop, has_msl0, has_msl4)
  
  ggsave(
    filename = file.path(output_dir, paste0(pop, "_tajd_genome_boxplot.jpeg")),
    plot = p, width = 14, height = 6, dpi = 300
  )
  ggsave(
    filename = file.path(output_dir, paste0(pop, "_tajd_genome_boxplot.pdf")),
    plot = p, width = 14, height = 6
  )
  
  cat("  Saved:", pop, "\n")
}

cat("\nAll done. Figures saved to:", output_dir, "\n")

