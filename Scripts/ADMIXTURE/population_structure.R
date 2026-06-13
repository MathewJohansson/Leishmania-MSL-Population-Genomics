# ADMIXTURE POPULATION STRUCTURE ANALYSIS
# Leishmania infantum Population Genomics Study

# This script analyses ADMIXTURE results for three datasets:
# 1. All L. infantum + L. donovani outgroup (K=12)
# 2. All L. infantum-only (K=8 and K=13)
# 3. Americas-only L. infantum (K=13)

# Populations use distinct notations: 
# The L.infantum + outgroup subset uses LO notation (e.g., LO1, LO2, etc.).
# The L.infantum-only subset uses LI notation (e.g., LI1, LI2, etc.).
# The Americas-only subset uses AM notation (e.g., AM1, AM2, etc.).
# All population assignments derive from the three master membership tables.




# 1. SETUP AND CONFIGURATION ---------------------------------------------------

# 1.1. Load Required Libraries -------------------------------------------------

library(tidyverse)
library(readr)
library(gtools)
library(ggplot2)
library(stringr)
library(RColorBrewer)
library(colorspace)
library(forcats)
library(ggtext)


# 1.2. Define Project Directories ----------------------------------------------

# Set base directory - ONLY CHANGE THIS LINE to run on different machines
BASE_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

# Define subdirectories
METADATA_DIR <- file.path(BASE_DIR, "Metadata")
DATA_DIR     <- file.path(BASE_DIR, "Data/population_structure")
ADMIXTURE_DIR <- file.path(BASE_DIR, "ADMIXTURE")
SCRIPTS_DIR  <- file.path(BASE_DIR, "Scripts/ADMIXTURE")
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/ADMIXTURE")
OUTPUTS_DIR  <- file.path(ADMIXTURE_DIR, "Analysis Outputs")

# Create output directories if they don't exist
dir.create(file.path(FIGURES_DIR, "All Linfantum With Outgroup"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "Linfantum Only"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(FIGURES_DIR, "Americas Only"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUTS_DIR, "High Confidence Samples"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUTS_DIR, "Defined Populations"),
           recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUTS_DIR, "Population Labels"),
           recursive = TRUE, showWarnings = FALSE)

# Fix for temporary file errors on Linux systems
TEMP_DIR <- file.path(BASE_DIR, "R_temp")
dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
Sys.setenv(TMPDIR = TEMP_DIR)
options(bitmapType = "cairo")
cat("✓ Temporary directory set to:", TEMP_DIR, "\n")


# 1.3. Load Metadata -----------------------------------------------------------

metadata <- read_tsv(file.path(METADATA_DIR,
                               "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"))


# 1.4. Define Colour Palettes --------------------------------------------------

# COLOURBLIND-FRIENDLY PALETTE - Paul Tol's Bright + Muted Scheme
# Source: https://personal.sron.nl/~pault/
# Used consistently across ALL analyses and scripts.
# Exported at the end of this script as paul_tol_colorblind_palette_v1_v13.rds

paul_tol_palette <- c(
  "#4477AA",  # Population 1  - blue
  "#EE6677",  # Population 2  - red
  "#228833",  # Population 3  - green
  "#CCBB44",  # Population 4  - yellow
  "#66CCEE",  # Population 5  - cyan
  "#AA3377",  # Population 6  - purple
  "#BBBBBB",  # Population 7  - grey
  "#EE7733",  # Population 8  - orange
  "#009988",  # Population 9  - teal
  "#994F00",  # Population 10 - brown
  "#332288",  # Population 11 - indigo
  "#9955BB",  # Population 12 - violet
  "#FFAABB"   # Population 13 - light rose
)

# Set palettes for the three subsets - L.i + outgroup, L.i-only, and Americas-only
palette_LO <- setNames(paul_tol_palette[1:12], paste0("LO", 1:12))
palette_LI <- setNames(paul_tol_palette[1:13], paste0("LI", 1:13))
palette_AM <- setNames(paul_tol_palette[1:13], paste0("AM", 1:13))




# 2. ALL L. INFANTUM + OUTGROUP (K=12) -----------------------------------------


# 2.1. Confidence Filtering and Plotting Clusters ------------------------------

# Load K=12 ADMIXTURE results
all_admix_data <- read_tsv(file.path(DATA_DIR,
                                     "population_membership_table_allsamples_k12_LO_notation.tsv"))

# Join metadata
all_admix_data <- all_admix_data %>%
  left_join(metadata, by = c("sample" = "sra_run_accession"))

# Filter by >=90% (high-confidence) assignment
high_conf_all <- all_admix_data %>%
  filter(population_proportion >= 0.9)

# Clean Host_Species column
high_conf_all <- high_conf_all %>%
  mutate(host_species_clean = case_when(
    grepl("Canis familiaris|Dog|Nyctomys procyonoides", host_species, ignore.case = TRUE) ~ "Dog",
    grepl("Homo sapiens|Human", host_species, ignore.case = TRUE) ~ "Human",
    TRUE ~ NA_character_
  ))

# Summarise per population
table_data_all_samples_k12 <- high_conf_all %>%
  group_by(population) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(unique(source_country), collapse = ", "),
    Host_Species = paste(unique(na.omit(host_species_clean)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(Notes = "") %>%
  select(population, N_Samples, Countries_Represented, Host_Species, Notes)

table_data_all_samples_k12

# Save population summary table
write_tsv(table_data_all_samples_k12,
          file.path(OUTPUTS_DIR,
                    "all_linfantum_plus_outgroup_population_summary_k12.tsv"))

# Save high-confidence samples
write_tsv(high_conf_all,
          file.path(OUTPUTS_DIR, "High Confidence Samples",
                    "all_linfantum_plus_outgroup_high_confidence_samples_k12.tsv"))

# Ensure population is a factor with correct notation order
high_conf_all <- high_conf_all %>%
  mutate(population = factor(population, levels = paste0("LO", 1:12)))

# Assignment barplot
ggplot(high_conf_all, aes(x = population, fill = population)) +
  geom_bar() +
  scale_fill_manual(values = palette_LO, drop = FALSE) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE-Inferred Population Assignments (K=12)<br>*L. infantum* Samples with *L. donovani* Outgroup",
    x = NULL,
    y = "Count",
    fill = "Population"
  ) +
  theme(
    plot.title = element_markdown(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save plots
ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "all_linfantum_with_outgroup_population_assignments_k12.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "all_linfantum_with_outgroup_population_assignments_k12.pdf"),
       width = 10, height = 6)


# The above has:
# 1. Read the ADMIXTURE K=12 population table for L. infantum + outgroup.
# 2. Filtered samples with >=90% confidence ancestry to a single population.
# 3. Saved the samples to: all_linfantum_plus_outgroup_high_confidence_samples_k12.tsv
# 4. Counted samples per population.
# 5. Generated and saved a bar plot using notation labels.




# 2.2. K=12 Justification ------------------------------------------------------

# Load the cross-validation data
cv_data_linfantum_with_outgroup <- read_table(
  file.path(DATA_DIR, "Linfantum_samples_with_outgroup_2025-03-12.CV.txt"),
  col_names = c("K", "CV_Error")
)
cv_data_linfantum_with_outgroup

# Plot CV error curve
ggplot(cv_data_linfantum_with_outgroup, aes(x = K, y = CV_Error)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE Cross-Validation Error vs. K",
    subtitle = "L. infantum Samples with L. donovani Outgroup",
    x = "Number of Populations (K)",
    y = "Cross-Validation Error"
  ) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "darkgreen") +
  annotate("text",
           x = 12.5,
           y = min(cv_data_linfantum_with_outgroup$CV_Error),
           label = "K = 12",
           color = "darkgreen",
           hjust = 0) 

# Save plots
ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "linfantum_samples_with_outgroup_cv_error_plot_k1_to_k20.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "linfantum_samples_with_outgroup_cv_error_plot_k1_to_k20.pdf"),
       width = 10, height = 6)




# 2.3. Stacked Ancestry Bar Plot -----------------------------------------------

# Read Q-matrix
q_file <- file.path(DATA_DIR,
                    "Linfantum_samples_with_outgroup_2025-03-12.q30.missing0.8.sampMissing0.25.biallelic.depth.annotate.pruned0.5_renamed.12.Q")
admix_wide <- read_table(q_file, col_names = FALSE)

# Add sample IDs
pop_membership_k12 <- read_tsv(file.path(DATA_DIR,
                                         "population_membership_table_allsamples_k12_LO_notation.tsv"))
stopifnot(nrow(admix_wide) == nrow(pop_membership_k12))

admix_wide <- admix_wide %>%
  mutate(sample = pop_membership_k12$sample) %>%
  relocate(sample)

# Name columns LO1..LO12
colnames(admix_wide)[-1] <- paste0("LO", 1:12)

# Pivot to long format
admix_long <- admix_wide %>%
  pivot_longer(cols = -sample,
               names_to  = "population",
               values_to = "population_proportion")

# Factor with correct ordering
admix_long <- admix_long %>%
  mutate(population = factor(population, levels = paste0("LO", 1:12)))

# Assign dominant population per sample for ordering
dominant_pop <- admix_long %>%
  group_by(sample) %>%
  slice_max(population_proportion, n = 1, with_ties = FALSE) %>%
  select(sample, dominant_population = population)

admix_long <- admix_long %>%
  left_join(dominant_pop, by = "sample")

# Reorder samples by dominant population then by max ancestry proportion
sample_order <- admix_long %>%
  group_by(sample) %>%
  summarise(
    dominant_population = first(dominant_population),
    max_ancestry = max(population_proportion)
  ) %>%
  arrange(dominant_population, desc(max_ancestry)) %>%
  pull(sample)

admix_long <- admix_long %>%
  mutate(sample = factor(sample, levels = sample_order))

# Plot
ggplot(admix_long, aes(x = sample, y = population_proportion, fill = population)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = palette_LO, drop = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_markdown(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  labs(
    title = "All *L. infantum* With Outgroup ADMIXTURE Ancestry Proportions (K=12)",
    x = "Samples Ordered by Population",
    y = "Ancestry Proportion",
    fill  = "Population"
  )

# Save plots
ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "all_linfantum_with_outgroup_ancestry_barplot.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "all_linfantum_with_outgroup_ancestry_barplot.pdf"),
       width = 10, height = 6)


# The above has:
# 1. Loaded ADMIXTURE Q-matrix (K=12) from the .Q file.
# 2. Added sample IDs by matching Q-matrix rows to metadata.
# 3. Pivoted to long format.
# 4. Determined plotting order by dominant ancestry cluster.
# 5. Plotted stacked ancestry barplot using LO notation labels.




# 3. ALL L. INFANTUM ONLY (K=13) -----------------------------------------------


# 3.1. Confidence Filtering and Plotting Clusters ------------------------------

# Load K=13 ADMIXTURE results
linfantum_k13_admix_data <- read_tsv(file.path(DATA_DIR,
                                               "population_membership_table_infantum_only_k13_LI_notation.tsv"))

# Filter by >=90% high-confidence and join metadata
high_conf_linfantum_k13 <- linfantum_k13_admix_data %>%
  filter(population_proportion >= 0.9) %>%
  left_join(metadata, by = c("sample" = "sra_run_accession"))

# Clean Host_Species column
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  mutate(host_species_clean = case_when(
    grepl("Canis familiaris|Dog|Nyctomys procyonoides", host_species, ignore.case = TRUE) ~ "Dog",
    grepl("Homo sapiens|Human", host_species, ignore.case = TRUE) ~ "Human",
    TRUE ~ NA_character_
  ))

# Summarise per population
table_data_linfantum_k13 <- high_conf_linfantum_k13 %>%
  group_by(population) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(unique(source_country), collapse = ", "),
    Host_Species = paste(unique(na.omit(host_species_clean)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(Notes = "") %>%
  select(population, N_Samples, Countries_Represented, Host_Species, Notes)

table_data_linfantum_k13

# Save outputs
write_tsv(high_conf_linfantum_k13,
          file.path(OUTPUTS_DIR, "High Confidence Samples",
                    "linfantum_only_high_confidence_samples_k13.tsv"))

write_tsv(table_data_linfantum_k13,
          file.path(OUTPUTS_DIR, "linfantum_only_population_summary_k13.tsv"))

# Ensure population is a factor with correct LI notation order
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  mutate(population = factor(population, levels = paste0("LI", 1:13)))

# Assignment barplot
ggplot(high_conf_linfantum_k13, aes(x = population, fill = population)) +
  geom_bar() +
  scale_fill_manual(values = palette_LI, drop = FALSE) +
  guides(fill = guide_legend(reverse = FALSE)) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE-Inferred Population Assignments (K=13)\n*L. infantum* Samples Only",
    x = NULL,
    y = "Number of Samples",
    fill  = "Population"
  ) +
  theme(
    plot.title = element_markdown(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save plots
ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_population_assignments_k13.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_population_assignments_k13.pdf"),
       width = 10, height = 6)


# The above has:
# 1. Read the ADMIXTURE K=13 population table for L. infantum only.
# 2. Filtered samples with >=90% confidence ancestry.
# 3. Saved high-confidence samples and population summary.
# 4. Generated and saved bar plot using notation labels.




# 3.2. K=13 Justification ------------------------------------------------------

cv_data_linfantum_only <- read_table(
  file.path(DATA_DIR, "infantum_only_n463.CV.txt"),
  col_names = c("K", "CV_Error")
)
cv_data_linfantum_only

ggplot(cv_data_linfantum_only, aes(x = K, y = CV_Error)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE Cross-Validation Error vs. K",
    subtitle = "All L. infantum Only",
    x = "Number of Populations (K)",
    y = "Cross-Validation Error"
  ) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "purple") +
  annotate("text",
           x = 13.5,
           y = min(cv_data_linfantum_only$CV_Error),
           label = "K = 13",
           color = "purple",
           hjust = 0)

# Save plots
ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_cv_error_plot_k1_to_k20_for_k13_only.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_cv_error_plot_k1_to_k20_for_k13_only.pdf"),
       width = 10, height = 6)




# 3.3. Define Populations ------------------------------------------------------

# Summarise sample counts per population, descending order
high_conf_linfantum_k13 %>%
  group_by(population) %>%
  summarise(num_samples = n()) %>%
  arrange(desc(num_samples))

# Join geographic region info
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  left_join(metadata %>% select(sra_run_accession, source_geographic_region),
            by = c("sample" = "sra_run_accession"))

# Clean up any duplicated columns from joins
high_conf_linfantum_k13 <- high_conf_linfantum_k13 %>%
  mutate(source_geographic_region = coalesce(source_geographic_region.x,
                                             source_geographic_region.y)) %>%
  select(-ends_with(".x"), -ends_with(".y"))

# Summarise by population and region
high_conf_linfantum_k13 %>%
  group_by(population, source_geographic_region) %>%
  summarise(num_samples = n(), .groups = "drop") %>%
  arrange(desc(num_samples))

# Final assignments table: sample, population (notation)
final_population_assignments_linfantum_only <- high_conf_linfantum_k13 %>%
  select(sample, population)
final_population_assignments_linfantum_only

write_tsv(final_population_assignments_linfantum_only,
          file.path(OUTPUTS_DIR, "Defined Populations",
                    "final_population_assignments_linfantum_only.tsv"))


# The above has:
# 1. Summarised sample counts per population.
# 2. Integrated geographic metadata.
# 3. Saved final population assignments using LI notation.




# 3.4. Stacked Ancestry Bar Plot -----------------------------------------------

# Load Q-matrix for K=13
Q_file_k13 <- file.path(DATA_DIR,
                        "infantum_only_n463.q30.missing0.8.sampMissing0.25.biallelic.depth.annotate.pruned0.5_renamed.13.Q")
Qmat_k13 <- read_table2(Q_file_k13, col_names = FALSE)

# Load membership table
pop_membership_k13 <- read_tsv(file.path(DATA_DIR,
                                         "population_membership_table_infantum_only_k13_LI_notation.tsv"))

# Verify sample counts match
stopifnot(nrow(Qmat_k13) == nrow(pop_membership_k13))

# Add sample names and population assignments
Qmat_k13_labeled <- Qmat_k13 %>%
  mutate(
    sample = pop_membership_k13$sample,
    population = pop_membership_k13$population
  ) %>%
  relocate(sample, population)

# Name Q-matrix columns LI1..LI13
colnames(Qmat_k13_labeled)[3:(2 + 13)] <- paste0("LI", 1:13)

# Order samples by population then by max ancestry
Qmat_k13_labeled <- Qmat_k13_labeled %>%
  rowwise() %>%
  mutate(max_ancestry = max(c_across(starts_with("LI")))) %>%
  ungroup() %>%
  arrange(population, desc(max_ancestry))

Qmat_k13_labeled$sample <- factor(Qmat_k13_labeled$sample,
                                  levels = Qmat_k13_labeled$sample)

# Pivot to long format
Qmat_k13_long <- Qmat_k13_labeled %>%
  pivot_longer(cols = starts_with("LI"),
               names_to = "ancestry_component",
               values_to = "proportion") %>%
  mutate(ancestry_component = factor(ancestry_component,
                                     levels = paste0("LI", 1:13)))

# Plot
ggplot(Qmat_k13_long,
       aes(x = sample, y = proportion, fill = ancestry_component)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = palette_LI, drop = FALSE) +
  theme_minimal() +
  labs(
    title = "*L. infantum*-only ADMIXTURE Ancestry Proportions (K=13)",
    x = "Samples Ordered by Population",
    y = "Ancestry Proportion",
    fill = "Population"
  ) +
  theme(
    plot.title = element_markdown(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save plots
ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_stacked_ancestry_barplot_k13.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_stacked_ancestry_barplot_k13.pdf"),
       width = 10, height = 6)


# The above has:
# 1. Loaded ADMIXTURE Q-matrix for K=13.
# 2. Added sample IDs and population assignments from the membership table.
# 3. Pivoted to long format.
# 4. Ordered samples by dominant population.
# 5. Plotted stacked ancestry barplot using LI notation labels.




# 4. AMERICAS ONLY (K=13) ------------------------------------------------------


# 4.1. Confidence Filtering and Plotting Clusters ------------------------------

# Load K=13 ADMIXTURE results for Americas Only
americas_only_k13_admix_data <- read_tsv(file.path(DATA_DIR,
                                                   "population_membership_table_Americas_only_k13_AM_notation.tsv"))

# Filter by >=90% high-confidence and join metadata
high_conf_americas_k13 <- americas_only_k13_admix_data %>%
  filter(population_proportion >= 0.9) %>%
  left_join(
    metadata %>% select(sra_run_accession, source_country,
                        source_state_or_region, host_species),
    by = c("sample" = "sra_run_accession")
  )

# Clean host species
high_conf_americas_k13 <- high_conf_americas_k13 %>%
  mutate(host_species_clean = case_when(
    grepl("Canis familiaris|Dog|Nyctomys procyonoides", host_species, ignore.case = TRUE) ~ "Dog",
    grepl("Homo sapiens|Human", host_species, ignore.case = TRUE) ~ "Human",
    TRUE ~ NA_character_
  ))

# Create population summary table
table_data_americas_k13 <- high_conf_americas_k13 %>%
  group_by(population) %>%
  summarise(
    N_Samples = n(),
    Countries_Represented = paste(sort(unique(source_country)), collapse = ", "),
    Regions_Represented = paste(sort(unique(na.omit(source_state_or_region))), collapse = ", "),
    Host_Species = paste(sort(unique(na.omit(host_species_clean))), collapse = ", "),
    Notes = "",
    .groups = "drop"
  ) %>%
  arrange(population)

# Save population summary table
write_tsv(table_data_americas_k13,
          file.path(OUTPUTS_DIR, "americas_only_population_summary_k13.tsv"))

# Ensure population is a factor with correct V-notation order
high_conf_americas_k13 <- high_conf_americas_k13 %>%
  mutate(population = factor(population, levels = paste0("AM", 1:13)))

# Assignment barplot
ggplot(high_conf_americas_k13, aes(x = population, fill = population)) +
  geom_bar() +
  scale_fill_manual(values = palette_AM,
                    breaks = paste0("AM", 1:13)) +
  labs(
    title = "ADMIXTURE-Inferred Population Assignments (K=13)\nAmericas-Only *L. infantum* Samples",
    x = NULL,
    y = "Number of Samples",
    fill = "Population"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_markdown(hjust = 0.5),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank()
  )

# Save plots
ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_linfantum_population_assignments_k13.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_linfantum_population_assignments_k13.pdf"),
       width = 10, height = 6)


# The above has:
# 1. Read the ADMIXTURE K=13 population table for Americas Only L. infantum.
# 2. Filtered samples with >=90% confidence ancestry.
# 3. Saved population summary table.
# 4. Generated and saved bar plot using AM notation labels.




# 4.2. K=13 Justification ------------------------------------------------------

cv_data_americas_only <- read_table(
  file.path(DATA_DIR, "Americas_samples_n388.CV.txt"),
  col_names = c("K", "CV_Error")
)
cv_data_americas_only

ggplot(cv_data_americas_only, aes(x = K, y = CV_Error)) +
  geom_line(color = "steelblue", size = 1.2) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal() +
  labs(
    title = "ADMIXTURE Cross-Validation Error vs. K",
    subtitle = "L. infantum Samples, Americas Only",
    x = "Number of Populations (K)",
    y = "Cross-Validation Error"
  ) +
  geom_vline(xintercept = 13, linetype = "dashed", color = "darkgreen") +
  annotate("text",
           x = 13.5,
           y = min(cv_data_americas_only$CV_Error),
           label = "K = 13",
           color = "darkgreen",
           hjust = 0)

# Save plots
ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_cv_error_plot_k1_to_k20.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_cv_error_plot_k1_to_k20.pdf"),
       width = 10, height = 6)




# 4.3. Define Populations ------------------------------------------------------

# Summarise sample counts per population
americas_only_k13_admix_data %>%
  group_by(population) %>%
  summarise(num_samples = n()) %>%
  arrange(desc(num_samples))

# Join source_country for geographic breakdown
high_conf_americas_k13 <- high_conf_americas_k13 %>%
  select(-starts_with("source_country")) %>%
  left_join(
    select(metadata, sra_run_accession, source_country),
    by = c("sample" = "sra_run_accession")
  )

# Summarise by population and country
high_conf_americas_k13 %>%
  group_by(population, source_country) %>%
  summarise(num_samples = n(), .groups = "drop") %>%
  arrange(population, desc(num_samples))

# Final assignments table: sample, population (V-notation)
final_population_assignments_americas_only <- high_conf_americas_k13 %>%
  select(sample, population)
final_population_assignments_americas_only

write_tsv(final_population_assignments_americas_only,
          file.path(OUTPUTS_DIR, "Defined Populations",
                    "final_population_assignments_americas_only.tsv"))


# The above has:
# 1. Summarised sample counts per population.
# 2. Geographic breakdown by country.
# 3. Saved final population assignments using AM notation.




# 4.4. Stacked Ancestry Bar Plot -----------------------------------------------

# Read Q-matrix for Americas only K=13
q_file_americas <- file.path(DATA_DIR,
                             "Americas_samples_n388.q30.missing0.8.sampMissing0.25.biallelic.depth.annotate.pruned0.5_renamed.13.Q")
admix_wide_americas <- read_table(q_file_americas, col_names = FALSE)

# Add sample names
admix_wide_americas$sample <- americas_only_k13_admix_data$sample

# Name population columns AM1..AM13
colnames(admix_wide_americas)[1:13] <- paste0("AM", 1:13)

# Pivot to long format
admix_long_americas <- admix_wide_americas %>%
  pivot_longer(cols = starts_with("AM"),
               names_to = "population",
               values_to = "population_proportion") %>%
  mutate(population = factor(population, levels = paste0("AM", 1:13)))

# Assign dominant population for ordering
dominant_pop_americas <- admix_long_americas %>%
  group_by(sample) %>%
  slice_max(population_proportion, n = 1, with_ties = FALSE) %>%
  select(sample, dominant_population = population)

admix_long_americas <- admix_long_americas %>%
  left_join(dominant_pop_americas, by = "sample")

# Reorder samples
sample_order_americas <- admix_long_americas %>%
  group_by(sample) %>%
  summarise(
    dominant_population = first(dominant_population),
    max_ancestry = max(population_proportion),
    .groups = "drop"
  ) %>%
  arrange(dominant_population, desc(max_ancestry)) %>%
  pull(sample)

admix_long_americas <- admix_long_americas %>%
  mutate(sample = factor(sample, levels = sample_order_americas))

# Plot
ancestry_plot_americas <- ggplot(admix_long_americas,
                                 aes(x = sample, y = population_proportion,
                                     fill = population)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values  = palette_AM,
                    breaks  = paste0("AM", 1:13),
                    name = "Population") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", size = 0.5),
    legend.position = "right",
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 11),
    plot.title = element_markdown(size = 14, hjust = 0.5)
  ) +
  labs(
    title = "Americas-only *L. infantum* ADMIXTURE Ancestry Proportions (K=13)",
    x = "Samples Ordered by Population",
    y = "Ancestry Proportion",
    fill = "Population"
  ) +
  guides(fill = guide_legend(ncol = 1))

print(ancestry_plot_americas)

# Save plots
ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_stacked_ancestry_barplot_k13.jpeg"),
       width = 10, height = 6)

ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_stacked_ancestry_barplot_k13.pdf"),
       width = 10, height = 6)




# Save combined palette list for use in geographic map scripts
paul_tol_palette_combined <- list(
  LO = palette_LO,
  LI = palette_LI,
  AM = palette_AM
)
saveRDS(paul_tol_palette_combined,
        file.path(BASE_DIR, "paul_tol_colorblind_palette_lo_li_am.rds"))
cat("Saved: paul_tol_colorblind_palette_lo_li_am.rds\n")




# 5. MSL STATUS COMPANION TRACKS -----------------------------------------------
# Stacked ancestry barplots with MSL coverage track beneath.
# Uses patchwork to combine panels. Sample order matches ancestry plots exactly.
# New figures only — originals above are unchanged.

library(patchwork)

# Shared MSL colour scale
coverage_colors <- c(
  "0" = "#d62728",
  "1" = "#ff7f0e",
  "2" = "#ffff00",
  "3" = "#2ca02c",
  "4" = "#1f77b4"
)

# Helper: build MSL track from a sample order vector
make_msl_track <- function(sample_order_vec, metadata) {
  msl_data <- data.frame(sample = sample_order_vec) %>%
    left_join(
      metadata %>%
        select(sra_run_accession, coverage_rounded) %>%
        mutate(coverage_rounded = as.character(coverage_rounded)),
      by = c("sample" = "sra_run_accession")
    ) %>%
    mutate(
      sample = factor(sample, levels = sample_order_vec),
      coverage_rounded = factor(coverage_rounded, levels = c("0","1","2","3","4"))
    )
  
  ggplot(msl_data, aes(x = sample, y = 1, fill = coverage_rounded)) +
    geom_tile(height = 1) +
    scale_fill_manual(values = coverage_colors, drop = FALSE,
                      name = "MSL copies") +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_minimal() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.title.x     = element_blank(),
      axis.title.y     = element_text(size = 8, angle = 90),
      panel.grid       = element_blank(),
      legend.position  = "right",
      legend.text      = element_text(size = 8),
      legend.title     = element_text(size = 9),
      legend.margin    = margin(t = 20),
      legend.key.size  = unit(0.4, "cm"),
      plot.margin      = margin(0, 5, 2, 5),
      axis.line.x.top  = element_line(colour = "black", linewidth = 0.8),
      legend.justification = "bottom"
    ) +
    labs(y = "MSL")
}


# 5.1. LO (K=12) Ancestry + MSL ------------------------------------------------

p_ancestry_lo <- ggplot(admix_long,
                        aes(x = sample, y = population_proportion,
                            fill = population)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = palette_LO, drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title         = element_markdown(hjust = 0.5),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin        = margin(5, 5, 0, 5),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.8),
    legend.justification = "bottom"
  ) +
  labs(
    title = "All *L. infantum* With Outgroup ADMIXTURE Ancestry Proportions (K=12)",
    x     = NULL,
    y     = "Ancestry Proportion",
    fill  = "Population"
  )

p_msl_lo <- make_msl_track(sample_order, metadata)

combined_lo <- p_ancestry_lo / p_msl_lo +
  plot_layout(heights = c(8, 1), guides = "collect") &
  theme(legend.justification = "bottom")

ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "all_linfantum_with_outgroup_ancestry_barplot_MSL.jpeg"),
       combined_lo, width = 12, height = 7, dpi = 300)

ggsave(file.path(FIGURES_DIR, "All Linfantum With Outgroup",
                 "all_linfantum_with_outgroup_ancestry_barplot_MSL.pdf"),
       combined_lo, width = 12, height = 7)


# 5.2. LI (K=13) Ancestry + MSL ------------------------------------------------

p_ancestry_li <- ggplot(Qmat_k13_long,
                        aes(x = sample, y = proportion,
                            fill = ancestry_component)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values = palette_LI, drop = FALSE) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal() +
  theme(
    plot.title         = element_markdown(hjust = 0.5),
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.margin        = margin(5, 5, 0, 5),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.8),
    legend.justification = "bottom"
  ) +
  labs(
    title = "*L. infantum*-only ADMIXTURE Ancestry Proportions (K=13)",
    x     = NULL,
    y     = "Ancestry Proportion",
    fill  = "Population"
  )

# LI sample order: taken from Qmat_k13_labeled which is already ordered
li_sample_order <- levels(Qmat_k13_long$sample)

p_msl_li <- make_msl_track(li_sample_order, metadata)

combined_li <- p_ancestry_li / p_msl_li +
  plot_layout(heights = c(8, 1), guides = "collect") &
  theme(legend.justification = "bottom")

ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_stacked_ancestry_barplot_k13_MSL.jpeg"),
       combined_li, width = 12, height = 7, dpi = 300)

ggsave(file.path(FIGURES_DIR, "Linfantum Only",
                 "linfantum_only_stacked_ancestry_barplot_k13_MSL.pdf"),
       combined_li, width = 12, height = 7)


# 5.3. AM (K=13) Ancestry + MSL ------------------------------------------------

p_ancestry_am <- ggplot(admix_long_americas,
                        aes(x = sample, y = population_proportion,
                            fill = population)) +
  geom_bar(stat = "identity", width = 1) +
  scale_fill_manual(values  = palette_AM,
                    breaks  = paste0("AM", 1:13),
                    name    = "Population") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    legend.position    = "right",
    legend.text        = element_text(size = 9),
    legend.title       = element_text(size = 11),
    plot.title         = element_markdown(size = 14, hjust = 0.5),
    plot.margin        = margin(5, 5, 0, 5),
    axis.line.x.bottom = element_line(colour = "black", linewidth = 0.8),
    legend.justification = "bottom"
  ) +
  labs(
    title = "Americas-only *L. infantum* ADMIXTURE Ancestry Proportions (K=13)",
    x     = NULL,
    y     = "Ancestry Proportion",
    fill  = "Population"
  ) +
  guides(fill = guide_legend(ncol = 1))

am_sample_order <- levels(admix_long_americas$sample)

p_msl_am <- make_msl_track(am_sample_order, metadata)

combined_am <- p_ancestry_am / p_msl_am +
  plot_layout(heights = c(8, 1), guides = "collect") &
  theme(legend.justification = "bottom")

ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_stacked_ancestry_barplot_k13_MSL.jpeg"),
       combined_am, width = 12, height = 7, dpi = 300)

ggsave(file.path(FIGURES_DIR, "Americas Only",
                 "americas_only_stacked_ancestry_barplot_k13_MSL.pdf"),
       combined_am, width = 12, height = 7)



