



# Turkish_Strain_Coverage_Verification.R
# Visualise read coverage across the MSL locus for 5 samples:
#   - 2 MSL0 (complete deletion, 0000) — controls
#   - 2 MSL4 (wild-type, 1111) — controls
#   - 1 MSL3 (partial deletion, 1110) — Turkish strain ERR3550134
# Part of: Evolutionary Forces Shaping the MSL Deletion in Brazilian L. infantum
# Author: Mathew Johansson
# Data: samtools bedcov output (100nt windows across MSL locus, chr31)




# 1. SETUP ------------------------------------------------------------------------

library(tidyverse)

# Directories
base_dir <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"
data_dir <- file.path(base_dir, "Turkish Strain Coverage Verification")
fig_dir  <- file.path(base_dir, "Figures/Turkish Strain Coverage Verification")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# MSL gene coordinates (from MSL-genes-only.bed)
msl_genes <- data.frame(
  gene  = c("Gene 1\n(LINF_310031200)", "Gene 2\n(LINF_310031300)",
            "Gene 3\n(LINF_310031400)", "Gene 4\n(LINF_310031500)"),
  start = c(1181282, 1184205, 1185827, 1191357),
  end   = c(1182328, 1185341, 1188553, 1192406)
)

# Sample metadata
samples <- data.frame(
  file      = c("coverage_SRR20748483.txt", "coverage_SRR20748320.txt",
                "coverage_SRR10478771.txt", "coverage_SRR10478777.txt",
                "coverage_ERR3550134.txt"),
  sample    = c("SRR20748483", "SRR20748320",
                "SRR10478771", "SRR10478777",
                "ERR3550134"),
  label     = c("SRR20748483 (MSL\u2212, Brazil) [control]",
                "SRR20748320 (MSL\u2212, Brazil) [control]",
                "SRR10478771 (MSL+, Brazil) [control]",
                "SRR10478777 (MSL+, Brazil) [control]",
                "ERR3550134 (MSL3/1110, T\u00FCrkiye)"),
  msl_group = c("MSL0", "MSL0", "MSL4", "MSL4", "MSL3"),
  colour    = c("#BB5566", "#CC0000", "#004488", "#0077BB", "#228833"),
  ltype     = c("solid", "dashed", "solid", "dashed", "solid"),
  stringsAsFactors = FALSE
)




# 2. LOAD DATA --------------------------------------------------------------------

load_coverage <- function(path, sample, label, msl_group) {
  read_tsv(path,
           col_names = c("CHROM", "BIN_START", "BIN_END", "TOTAL_BASES"),
           show_col_types = FALSE) %>%
    mutate(
      midpoint   = (BIN_START + BIN_END) / 2,
      mean_depth = TOTAL_BASES / 100,  # divide by window size to get mean per-base depth
      sample     = sample,
      label      = label,
      msl_group  = msl_group
    )
}

cov_data <- pmap_dfr(
  samples[, c("file", "sample", "label", "msl_group")],
  function(file, sample, label, msl_group) {
    load_coverage(file.path(data_dir, file), sample, label, msl_group)
  }
)

# Set factor order for legend
cov_data$label <- factor(cov_data$label, levels = samples$label)




# 3. PLOT -------------------------------------------------------------------------

colour_values <- setNames(samples$colour, samples$label)
ltype_values  <- setNames(samples$ltype,  samples$label)

p <- ggplot(cov_data, aes(x = midpoint, y = mean_depth,
                          colour = label, linetype = label)) +
  
  # MSL gene shading
  geom_rect(data = msl_genes,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE,
            fill = "#DDDDDD", alpha = 0.5) +
  
  # Gene labels
  geom_text(data = msl_genes,
            aes(x = (start + end) / 2, y = Inf, label = gene),
            inherit.aes = FALSE,
            vjust = 1.3, size = 2.8, colour = "grey35", lineheight = 0.85) +
  
  # Raw coverage lines (transparent)
  geom_line(linewidth = 0.4, alpha = 0.3) +
  
  # Loess smoothed lines (prominent)
  geom_smooth(method = "loess", span = 0.15, se = FALSE, linewidth = 1.1) +
  
  scale_colour_manual(name = "Sample", values = colour_values) +
  scale_linetype_manual(name = "Sample", values = ltype_values) +
  
  labs(
    title    = "Read Coverage Across the MSL Locus (chr31)",
    subtitle = "100nt windows | MSL genes shaded grey | ERR3550134: partial deletion (1110), T\u00FCrkiye",
    x        = "Position on chr31 (bp)",
    y        = "Mean read depth per base"
  ) +
  
  theme_bw() +
  theme(
    panel.grid.minor    = element_blank(),
    panel.grid.major.x  = element_blank(),
    axis.text           = element_text(size = 11),
    axis.title          = element_text(size = 12),
    legend.title        = element_text(size = 12, face = "bold"),
    legend.text         = element_text(size = 10),
    legend.position     = "right",
    plot.title          = element_text(size = 13, face = "bold", hjust = 0.5),
    plot.title.position = "plot",
    plot.subtitle       = element_text(size = 10, hjust = 0.5, colour = "grey40")
  )

# Display plot
p




# 4. SAVE -------------------------------------------------------------------------

# Save the plot as both jpeg and pdf images
ggsave(
  filename = file.path(fig_dir, "turkish_strain_MSL_coverage.jpeg"),
  plot     = p,
  width    = 12, height = 6, dpi = 300
)

ggsave(
  filename = file.path(fig_dir, "turkish_strain_MSL_coverage.pdf"),
  plot     = p,
  width    = 12, height = 6
)

message("Figures saved to: ", fig_dir)



