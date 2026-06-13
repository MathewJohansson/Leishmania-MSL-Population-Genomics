# GEOGRAPHICAL POPULATION STRUCTURE VISUALISATION
# L.infantum Genomics Study

# This script creates geographical maps showing ADMIXTURE population structure:
# 1. All L.infantum + L.donovani outgroup (K=12)
# 2. All L.infantum only (K=13)
# 3. Americas-only L.infantum (K=13)

# Populations use distinct notations: 
# The L.infantum + outgroup subset uses LO notation (e.g., LO1, LO2, etc.).
# The L.infantum-only subset uses LI notation (e.g., LI1, LI2, etc.).
# The Americas-only subset uses AM notation (e.g., AM1, AM2, etc.).
# All population assignments derive from the three master membership tables.




# 1. SETUP AND CONFIGURATION ---------------------------------------------------

# 1.1. Load Required Libraries -------------------------------------------------

library(tidyverse)    
library(ggplot2)      
library(RColorBrewer) 
library(gridExtra)    
library(sf)           
library(rnaturalearth) 
library(rnaturalearthdata) 
library(cowplot)      
library(readr)        


# 1.2. Define Project Directories ----------------------------------------------

# Set base directory - ONLY CHANGE THIS LINE to run on different machines
BASE_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

# Define subdirectories
METADATA_DIR <- file.path(BASE_DIR, "Metadata")
DATA_DIR <- file.path(BASE_DIR, "Data/population_structure")
FIGURES_DIR <- file.path(BASE_DIR, "Figures/Geographical Work")

# Create output directory if needed
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)


# 1.3. Load Colour-blind-Friendly Palette --------------------------------------

# Load the Paul Tol colour-blind-friendly palette exported by Mat's script
# This ensures consistent colours across all figures for LO/LI/AM populations
paul_tol_palette_path <- file.path(BASE_DIR, "paul_tol_colorblind_palette_lo_li_am.rds")

if (file.exists(paul_tol_palette_path)) {
  paul_tol_palette <- readRDS(paul_tol_palette_path)
  if (is.list(paul_tol_palette)) {
    paul_tol_palette <- unlist(paul_tol_palette, use.names = TRUE)
    names(paul_tol_palette) <- sub("^[^.]+\\.", "", names(paul_tol_palette))
  }
  cat("✓ Loaded Paul Tol colorblind-friendly palette\n")
} else {
  # Fallback: Define palette directly if .rds file not found
  cat("⚠ Warning: paul_tol_colorblind_palette_lo_li_am.rds not found\n")
  cat("  Using fallback palette definition\n")
  tol_colours <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB",
                   "#EE7733", "#009988", "#882255", "#332288", "#DDCC77", "#44BB99")
  paul_tol_palette <- c(
    setNames(tol_colours[1:12], paste0("LO", 1:12)),
    setNames(tol_colours[1:13], paste0("LI", 1:13)),
    setNames(tol_colours[1:13], paste0("AM", 1:13))
  )
}


# 1.4. Define Data Standardisation Mappings ------------------------------------

# Country name standardisation
country_fix <- tribble(
  ~source_country, ~standard_country,
  "Brazil", "Federative Republic of Brazil",
  "Paraguay", "Republic of Paraguay",
  "Uruguay", "Oriental Republic of Uruguay",
  "China", "People's Republic of China",
  "Uzbekistan", "Republic of Uzbekistan",
  "Albania", "Republic of Albania",
  "France", "French Republic",
  "Greece", "Hellenic Republic",
  "Italy", "Italian Republic",
  "Portugal", "Portuguese Republic",
  "Spain", "Kingdom of Spain",
  "Yugoslavia", "Republic of Serbia",
  "Israel", "State of Israel",
  "Turkiye", "Republic of Turkey",
  "Yemen", "Republic of Yemen",
  "Egypt", "Arab Republic of Egypt",
  "Morocco", "Kingdom of Morocco",
  "Tunisia", "Tunisian Republic",
  "Honduras", "Republic of Honduras",
  "Panama", "Republic of Panama",
  "Palestine", "Palestine",
  "Ethiopia", "Federal Democratic Republic of Ethiopia",
  "Sudan", "Republic of the Sudan",
  "Georgia", "Georgia",
  "Syria", "Syria",
  "Iran", "Iran"
)

# Brazil state mapping
state_mapping <- tribble(
  ~full_name, ~code,
  "Acre", "AC", "Alagoas", "AL", "Amapá", "AP", "Amazonas", "AM", "Bahia", "BA",
  "Ceara", "CE", "Ceará", "CE", "Distrito Federal", "DF", "Espírito Santo", "ES", 
  "Goiás", "GO", "Maranhão", "MA", "Mato Grosso", "MT", "Mato Grosso do Sul", "MS", 
  "Minas Gerais", "MG", "Pará", "PA", "Paraíba", "PB", "Paraná", "PR", 
  "Pernambuco", "PE", "Piauí", "PI", "Rio de Janeiro", "RJ", "Rio Grande do Norte", "RN", 
  "Rio Grande do Sul", "RS", "Rondônia", "RO", "Roraima", "RR", "Santa Catarina", "SC", 
  "São Paulo", "SP", "Sergipe", "SE", "Tocantins", "TO"
)


# 1.5. Load Data Files ---------------------------------------------------------

# Load sample metadata
samples_data <- read_tsv(file.path(METADATA_DIR, "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"), 
                         show_col_types = FALSE) %>%
  rename(sample = sra_run_accession)

# Load population membership files
pop_membership_Americas_k13 <- read_tsv(file.path(DATA_DIR, "population_membership_table_Americas_only_k13_AM_notation.tsv"), 
                                        show_col_types = FALSE)
pop_membership_Infantum_k13 <- read_tsv(file.path(DATA_DIR, "population_membership_table_infantum_only_k13_LI_notation.tsv"), 
                                        show_col_types = FALSE)
pop_membership_All_k12 <- read_tsv(file.path(DATA_DIR, "population_membership_table_allsamples_k12_LO_notation.tsv"), 
                                   show_col_types = FALSE)
pop_membership_Infantum_k8 <- read_tsv(file.path(DATA_DIR, "population_membership_table_infantum_only_k8_LI_notation.tsv"), 
                                       show_col_types = FALSE)

# Create a named list with loaded population dataframes
pop_data_list <- list(
  "Americas_only_k13" = pop_membership_Americas_k13,
  "infantum_only_k13" = pop_membership_Infantum_k13,
  "allsamples_k12" = pop_membership_All_k12,
  "infantum_only_k8" = pop_membership_Infantum_k8
)




# 2. HELPER FUNCTIONS ----------------------------------------------------------

# 2.1. Data Processing Function ------------------------------------------------

# Function to process population membership data and merge with sample locations
process_population_data <- function(pop_data, dataset_name, samples_data) {
  cat("\n=== Processing", dataset_name, "===\n")
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("No population data provided for processing (", dataset_name, ").\n")
    return(NULL)
  }
  
  # Filter for high confidence assignment (>=90%)
  # Consistent with ADMIXTURE barplots in population_structure.R
  pop_filtered <- pop_data %>% filter(population_proportion >= 0.9)
  cat("Samples with >=90% population proportion:", nrow(pop_filtered), "\n")
  
  if (nrow(pop_filtered) == 0) {
    cat("No samples with >=90% population proportion found after filtering for", dataset_name, ".\n")
    return(NULL)
  }
  
  # Clean samples data
  samples_clean <- samples_data %>%
    left_join(country_fix, by = "source_country", relationship = "many-to-many") %>%
    mutate(source_country = coalesce(standard_country, source_country)) %>%
    select(-standard_country) %>%
    left_join(state_mapping, by = c("source_state_or_region" = "full_name"), 
              relationship = "many-to-many") %>%
    mutate(source_state_or_region = if_else(!is.na(code), code, source_state_or_region)) %>%
    select(-code)
  
  # Join population data with sample location data
  combined_data <- pop_filtered %>%
    left_join(samples_clean, by = "sample", relationship = "many-to-many")
  
  if (nrow(combined_data) == 0) {
    cat("WARNING: Join resulted in 0 rows for", dataset_name, ".\n")
    return(NULL)
  }
  
  # Filter for Brazil states
  combined_data_brazil <- combined_data %>%
    filter(source_country == "Federative Republic of Brazil" | source_country == "Brazil") %>%
    filter(!is.na(source_state_or_region))
  
  cat("Brazil samples with population data:", nrow(combined_data_brazil), "\n")
  
  if (nrow(combined_data_brazil) == 0) {
    cat("No Brazil samples found with population data.\n")
    return(NULL)
  }
  
  # Create population summary by state
  state_populations <- combined_data_brazil %>%
    group_by(source_state_or_region, population) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(source_state_or_region) %>%
    mutate(
      total_samples = sum(count),
      proportion = count / total_samples
    ) %>%
    ungroup() %>%
    filter(total_samples >= 2) %>%
    arrange(desc(total_samples))
  
  cat("States with sufficient data:", length(unique(state_populations$source_state_or_region)), "\n")
  
  return(state_populations)
}




# 3. MAP VISUALISATION ---------------------------------------------------------

# 3.1. Map Creation Function ---------------------------------------------------

# Function to create map with properly sized uniform bar charts
create_map_with_bar_charts <- function(state_pop_data, dataset_name, k_value) {
  # Load Brazil map data using scale 50 (medium resolution, included in base package)
  # This avoids the rnaturalearthhires package requirement
  brazil_sf <- ne_download(scale = 50, type = 'states', category = 'cultural', returnclass = "sf")
  brazil_states <- brazil_sf[brazil_sf$admin == "Brazil", ]
  
  # Standardise state codes in the map data
  brazil_states <- brazil_states %>%
    mutate(code = case_when(
      name == "Acre" ~ "AC", name == "Alagoas" ~ "AL", name == "Amapá" ~ "AP",
      name == "Amazonas" ~ "AM", name == "Bahia" ~ "BA", name == "Ceará" ~ "CE",
      name == "Distrito Federal" ~ "DF", name == "Espírito Santo" ~ "ES", name == "Goiás" ~ "GO",
      name == "Maranhão" ~ "MA", name == "Mato Grosso" ~ "MT", name == "Mato Grosso do Sul" ~ "MS",
      name == "Minas Gerais" ~ "MG", name == "Pará" ~ "PA", name == "Paraíba" ~ "PB",
      name == "Paraná" ~ "PR", name == "Pernambuco" ~ "PE", name == "Piauí" ~ "PI",
      name == "Rio de Janeiro" ~ "RJ", name == "Rio Grande do Norte" ~ "RN",
      name == "Rio Grande do Sul" ~ "RS", name == "Rondônia" ~ "RO", name == "Roraima" ~ "RR",
      name == "Santa Catarina" ~ "SC", name == "São Paulo" ~ "SP", name == "Sergipe" ~ "SE",
      name == "Tocantins" ~ "TO",
      TRUE ~ as.character(abbrev)
    )) %>%
    select(name, code, geometry) %>%
    st_make_valid()
  
  if (is.null(state_pop_data) || nrow(state_pop_data) == 0) {
    cat("DEBUG: state_pop_data is NULL or empty.\n")
    return(NULL)
  }
  
  cat("Creating map for", dataset_name, "with", nrow(state_pop_data), "data points\n")
  
  # Calculate centroids
  brazil_states_filtered <- brazil_states %>%
    filter(code %in% unique(state_pop_data$source_state_or_region))
  
  if (nrow(brazil_states_filtered) == 0) {
    cat("No matching states found in map data.\n")
    return(NULL)
  }
  
  brazil_states_centroids <- brazil_states_filtered %>%
    st_centroid() %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(lon = X, lat = Y) %>%
    bind_cols(brazil_states_filtered %>% st_drop_geometry()) %>%
    select(name, code, lon, lat)
  
  # Merge population data with centroids
  plot_data <- state_pop_data %>%
    left_join(brazil_states_centroids, by = c("source_state_or_region" = "code")) %>%
    filter(!is.na(lon))
  
  if (nrow(plot_data) == 0) {
    cat("No plot data after joining with centroids.\n")
    return(NULL)
  }
  
  # Create colour palette using Paul Tol colour-blind-friendly scheme
  # Get unique populations and sort them NUMERICALLY (LO1, LO2... / LI1... / AM1...)
  all_populations <- unique(plot_data$population)
  # Sort numerically by extracting the trailing number from LO1, LI1, AM1, etc.
  all_populations <- all_populations[order(as.numeric(gsub("^[A-Za-z]+", "", all_populations)))]
  n_pops <- length(all_populations)
  
  # DEBUG: Print the ordering
  cat("Populations in numeric order:", paste(all_populations, collapse = ", "), "\n")
  
  # Convert population to factor with numeric ordering BEFORE plotting
  plot_data <- plot_data %>%
    mutate(population = factor(population, levels = all_populations))
  
  # DEBUG: Check factor levels
  cat("Factor levels:", paste(levels(plot_data$population), collapse = ", "), "\n")
  
  # Use Paul Tol palette for LO/LI/AM populations
  colors <- paul_tol_palette[all_populations]
  names(colors) <- all_populations  # Explicitly set names to preserve order
  
  # Verify all populations got colours
  if (any(is.na(colors))) {
    cat("WARNING: Some populations not in palette, using fallback colours\n")
    missing_pops <- all_populations[is.na(colors)]
    cat("Missing populations:", paste(missing_pops, collapse = ", "), "\n")
    # Fallback for any unrecognised populations
    colors[is.na(colors)] <- rainbow(sum(is.na(colors)))
  }
  
  # Base map
  p_map <- ggplot() +
    geom_sf(data = brazil_states, fill = "lightgray", color = "white", linewidth = 0.2) +
    coord_sf() +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold", 
                                margin = margin(t = 5, b = 3)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 5)),
      legend.position = "none",
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5)  
    ) +
    labs(
      title = "Population Structure in Brazilian States",
      subtitle = paste(dataset_name, "- K =", k_value, "| Samples with >=90% assignment confidence")
    )
  
  # Medium-sized uniform bar charts - same size for all states
  unique_states <- unique(plot_data$source_state_or_region)
  cat("Processing", length(unique_states), "states for bar charts\n")
  
  # Fixed bar chart dimensions
  fixed_width <- 1.4   
  fixed_height <- 2.0  
  
  for (state_code in unique_states) {
    state_df <- plot_data %>% filter(source_state_or_region == state_code)
    
    if (nrow(state_df) == 0) next
    
    state_df <- state_df %>% arrange(population)
    total_samples <- unique(state_df$total_samples)
    
    cat("State:", state_code, "Samples:", total_samples, 
        "Fixed size:", fixed_width, "x", fixed_height, "\n")
    
    # Create medium-sized bar chart
    state_bar_chart <- ggplot(state_df, aes(x = 1, y = proportion, fill = population)) +
      geom_col(width = 1, color = "black", linewidth = 0.15) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.4)
      )
    
    # Get centroid coordinates
    centroid_lon <- unique(state_df$lon)
    centroid_lat <- unique(state_df$lat)
    
    # Add bar chart to map with fixed dimensions
    p_map <- p_map +
      annotation_custom(
        grob = ggplotGrob(state_bar_chart),
        xmin = centroid_lon - fixed_width/2,
        xmax = centroid_lon + fixed_width/2,
        ymin = centroid_lat - fixed_height/2,
        ymax = centroid_lat + fixed_height/2
      )
  }
  
  # Create separate legend with NUMERIC ordering
  dummy_data <- tibble(
    population = factor(all_populations, levels = all_populations),  # Force numeric order
    y = 1
  )
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = colors, 
      name = "Population", 
      breaks = all_populations,
      limits = all_populations,  # Force these exact levels
      drop = FALSE
    ) +
    guides(fill = guide_legend(override.aes = list(size = 1))) +  # Explicit guide control
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.5, "cm"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11, face = "bold"),
      legend.margin = margin(5, 5, 5, 5)
    )
  
  # Extract legend - use ggplotGrob instead of get_legend to preserve factor order
  legend_grob <- ggplotGrob(legend_plot)
  
  # Find the legend in the grob
  legend_index <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  legend <- legend_grob$grobs[[legend_index]]
  
  # Combine map and legend with proper spacing
  final_plot <- plot_grid(
    p_map, 
    legend, 
    ncol = 2, 
    rel_widths = c(0.75, 0.25),  
    align = "h"
  ) +
    theme(plot.margin = margin(t = 5, r = 5, b = 5, l = 5))  
  
  return(final_plot)
}




# 4. MAIN EXECUTION ------------------------------------------------------------

cat("BRAZIL POPULATION STRUCTURE MAP ANALYSIS")

# Check if samples_data is loaded correctly
if (!exists("samples_data") || !("sample" %in% names(samples_data))) {
  cat("ERROR: samples_data not found or missing 'sample' column.\n")
  stop("Please check data loading.")
}

cat("Samples data loaded with", nrow(samples_data), "rows\n")

# Debug: Check each dataset before processing
cat("\n=== DEBUGGING: Checking all datasets ===\n")
for (i in 1:length(pop_data_list)) {
  file_name <- names(pop_data_list)[i]
  pop_data <- pop_data_list[[i]]
  
  cat("Dataset", i, ":", file_name, "\n")
  if (is.null(pop_data)) {
    cat("  - STATUS: NULL (file not loaded or empty)\n")
  } else if (nrow(pop_data) == 0) {
    cat("  - STATUS: Empty (0 rows)\n")
  } else {
    cat("  - STATUS: OK -", nrow(pop_data), "rows loaded\n")
    cat("  - Columns:", paste(names(pop_data), collapse = ", "), "\n")
    if ("population_proportion" %in% names(pop_data)) {
      high_confidence <- sum(pop_data$population_proportion >= 0.9, na.rm = TRUE)
      cat("  - Samples with >=90% confidence:", high_confidence, "\n")
    }
  }
  cat("\n")
}

# Process each dataset and create plots
plots_created <- 0
failed_datasets <- c()

for (i in 1:length(pop_data_list)) {
  file_name <- names(pop_data_list)[i]
  pop_data <- pop_data_list[[i]]
  
  cat("\n", rep("=", 60), "\n")
  cat("PROCESSING DATASET", i, "of", length(pop_data_list), ":", file_name, "\n")
  cat(rep("=", 60), "\n")
  
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("SKIPPING", file_name, "- no data loaded\n")
    failed_datasets <- c(failed_datasets, file_name)
    next
  }
  
  # Extract dataset info
  k_value <- case_when(
    str_detect(file_name, "k13") ~ "13",
    str_detect(file_name, "k12") ~ "12", 
    str_detect(file_name, "k8") ~ "8",
    TRUE ~ "unknown"
  )
  
  dataset_type <- case_when(
    str_detect(file_name, "Americas_only") ~ "Americas Only",
    str_detect(file_name, "infantum_only") ~ "L. infantum Only", 
    str_detect(file_name, "allsamples") ~ "All Samples",
    TRUE ~ "Dataset"
  )
  
  cat("Dataset type:", dataset_type, "\n")
  cat("K value:", k_value, "\n")
  cat("Raw data rows:", nrow(pop_data), "\n")
  
  # Process data
  state_pop_data <- process_population_data(pop_data, dataset_type, samples_data)
  
  if (is.null(state_pop_data)) {
    cat("FAILED: No processed data for", dataset_type, "\n")
    failed_datasets <- c(failed_datasets, file_name)
    next
  }
  
  # Create map
  map_plot <- create_map_with_bar_charts(state_pop_data, dataset_type, k_value)
  
  if (!is.null(map_plot)) {
    # Display the plot
    print(map_plot)
    
    # Save plot as JPEG
    safe_name <- str_replace_all(tolower(dataset_type), "[^a-z0-9]+", "_") 
    filename <- paste0("brazil_population_map_paper_", safe_name, "_k", k_value, ".jpg")
    filepath <- file.path(FIGURES_DIR, filename)
    ggsave(filepath, map_plot, width = 14, height = 11, dpi = 300)
    cat("SUCCESS: Plot saved as", filepath, "\n")
    plots_created <- plots_created + 1
    
    # Add a pause to see each plot
    if (interactive()) {
      readline(prompt = "Press [Enter] to continue to next plot...")
    }
  } else {
    cat("FAILED: Could not create plot for", dataset_type, "\n")
    failed_datasets <- c(failed_datasets, file_name)
  }
}




# 5. SUMMARY -------------------------------------------------------------------

cat("\n=== FINAL SUMMARY ===\n")
cat("Expected datasets: 4\n")
cat("Total plots successfully created:", plots_created, "\n")
cat("Failed datasets:", length(failed_datasets), "\n")
if (length(failed_datasets) > 0) {
  cat("Failed dataset names:\n")
  for (failed in failed_datasets) {
    cat("  -", failed, "\n")
  }
  cat("\nTROUBLESHOOTING SUGGESTIONS:\n")
  cat("1. Check that all 4 .tsv files exist in the 'data/' directory\n")
  cat("2. Ensure file names match exactly (case sensitive)\n")
  cat("3. Check file permissions and that files aren't empty\n")
  cat("4. Verify file paths are correct\n")
}
cat("Script completed!\n")

save.image("map_bar_plotting_data.Rda")




# 6. EURASIA/AFRICA POPULATION STRUCTURE MAPS ----------------------------------
# Generates Europe, Asia & Africa focused maps for all datasets.
# Uses the same paul_tol_palette and FIGURES_DIR as Brazil maps above.

# 6.1 Process population data worldwide (by country/location, not Brazilian state)
process_population_data_worldwide <- function(pop_data, dataset_name, samples_data) {
  cat("\n=== Processing WORLDWIDE", dataset_name, "===\n")
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("No population data provided for processing (", dataset_name, ").\n")
    return(NULL)
  }
  
  # Filter for high confidence assignment (>0.8)
  pop_filtered <- pop_data %>% filter(population_proportion >= 0.9)
  cat("Samples with >=90% population proportion:", nrow(pop_filtered), "\n")
  
  if (nrow(pop_filtered) == 0) {
    cat("No samples with >=90% population proportion found after filtering for", dataset_name, ".\n")
    return(NULL)
  }
  
  # Clean samples data with existing country fixes
  samples_clean <- samples_data %>%
    left_join(country_fix, by = "source_country", relationship = "many-to-many") %>%
    mutate(source_country = coalesce(standard_country, source_country)) %>%
    select(-standard_country)
  
  # Join population data with sample location data
  combined_data <- pop_filtered %>%
    left_join(samples_clean, by = "sample", relationship = "many-to-many")
  
  if (nrow(combined_data) == 0) {
    cat("WARNING: Join resulted in 0 rows for", dataset_name, ".\n")
    return(NULL)
  }
  
  # Worldwide: don't filter for Brazil only
  combined_data_worldwide <- combined_data %>%
    filter(!is.na(source_country)) %>%
    mutate(
      location_id = case_when(
        (source_country == "Federative Republic of Brazil" | source_country == "Brazil") & 
          !is.na(source_state_or_region) ~ paste0("BR-", source_state_or_region),
        TRUE ~ source_country
      ),
      location_name = case_when(
        str_starts(location_id, "BR-") ~ paste0("Brazil (", source_state_or_region, ")"),
        TRUE ~ source_country
      )
    )
  
  cat("Worldwide samples with population data:", nrow(combined_data_worldwide), "\n")
  
  if (nrow(combined_data_worldwide) == 0) {
    cat("No worldwide samples found with population data.\n")
    return(NULL)
  }
  
  # Create population summary by location
  # For non-Brazil locations, group by country only (ignore sub-country variation)
  location_populations <- combined_data_worldwide %>%
    group_by(location_id, location_name, population, source_country) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(location_id) %>%
    mutate(
      total_samples = sum(count),
      proportion = count / total_samples
    ) %>%
    ungroup() %>%
    mutate(source_state_or_region = NA_character_) %>%
    filter(total_samples >= 1) %>%
    arrange(desc(total_samples))
  
  cat("Locations with sufficient data:", length(unique(location_populations$location_id)), "\n")
  cat("Countries represented:", length(unique(location_populations$source_country)), "\n")
  
  return(location_populations)
}


# 6.2 Create Europe/Asia/Africa focused map
create_europe_asia_focused_map <- function(location_pop_data, dataset_name, k_value) {
  
  european_countries <- c("Spain", "Kingdom of Spain", "Portuguese Republic", "Portugal", 
                          "French Republic", "France", "Italian Republic", "Italy", 
                          "Hellenic Republic", "Greece", "Republic of Albania", "Albania",
                          "Republic of Serbia", "Yugoslavia", "Serbia")
  
  asian_countries <- c("People's Republic of China", "China", "Republic of Turkey", "Turkiye", 
                       "Turkey", "State of Israel", "Israel", "Republic of Yemen", "Yemen",
                       "Republic of Uzbekistan", "Uzbekistan")
  
  african_countries <- c("Arab Republic of Egypt", "Egypt", "Kingdom of Morocco", "Morocco",
                         "Tunisian Republic", "Tunisia", "Federal Democratic Republic of Ethiopia", 
                         "Ethiopia", "Republic of the Sudan", "Sudan")
  
  eurasia_africa_data <- location_pop_data %>%
    filter(source_country %in% c(european_countries, asian_countries, african_countries) |
             str_detect(source_country, "Spain|Portugal|France|Italy|Greece|Albania|Serbia|China|Turkey|Israel|Yemen|Uzbekistan|Egypt|Morocco|Tunisia|Ethiopia|Sudan"))
  
  cat("Europe/Asia/Africa samples found:", nrow(eurasia_africa_data), "\n")
  cat("Unique locations:", length(unique(eurasia_africa_data$location_id)), "\n")
  
  if (nrow(eurasia_africa_data) == 0) {
    cat("No Europe/Asia/Africa data found\n")
    return(NULL)
  }
  
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Direct coordinate lookup using short country names
  # Avoids the broken pattern-matching fallback that misplaces countries
  short_name_lookup <- tribble(
    ~long_name,                                          ~short_name,
    "Tunisian Republic",                                 "Tunisia",
    "State of Israel",                                   "Israel",
    "French Republic",                                   "France",
    "Kingdom of Spain",                                  "Spain",
    "Republic of Uzbekistan",                            "Uzbekistan",
    "Italian Republic",                                  "Italy",
    "Kingdom of Morocco",                                "Morocco",
    "People's Republic of China",                       "China",
    "Portuguese Republic",                               "Portugal",
    "Arab Republic of Egypt",                            "Egypt",
    "Hellenic Republic",                                 "Greece",
    "Republic of Albania",                               "Albania",
    "Republic of Turkey",                                "Turkey",
    "Republic of Yemen",                                 "Yemen",
    "Republic of Serbia",                                "Serbia",
    "Federal Democratic Republic of Ethiopia",           "Ethiopia",
    "Republic of the Sudan",                             "Sudan",
    "Georgia",                                           "Georgia",
    "Syria",                                             "Syria",
    "Palestine",                                         "Palestine",
    "Iran",                                              "Iran"
  )
  
  # Pre-compute world centroids once
  world_centroids <- world_map %>%
    st_centroid() %>%
    mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    select(sovereignt, name, lon, lat)
  
  unique_locations <- unique(eurasia_africa_data$location_id)
  cat("Processing coordinates for", length(unique_locations), "unique locations\n")
  
  # Manual overrides for countries whose rnaturalearth centroids are wrong
  # (e.g. France includes overseas territories pulling centroid southwest)
  manual_coords <- tribble(
    ~location_id,              ~lon,   ~lat,
    "French Republic",          2.2,   46.2,
    "Portuguese Republic",     -8.0,   39.6
  )
  
  location_coords <- eurasia_africa_data %>%
    group_by(location_id, location_name, source_country) %>%
    slice(1) %>%
    ungroup() %>%
    left_join(short_name_lookup, by = c("source_country" = "long_name")) %>%
    mutate(lookup_name = coalesce(short_name, source_country)) %>%
    left_join(world_centroids, by = c("lookup_name" = "name")) %>%
    left_join(manual_coords, by = "location_id", suffix = c("", "_manual")) %>%
    mutate(
      lon = coalesce(lon_manual, lon),
      lat = coalesce(lat_manual, lat)
    ) %>%
    select(-lon_manual, -lat_manual) %>%
    filter(!is.na(lon)) %>%
    select(location_id, location_name, lon, lat)
  
  cat("Successfully processed", nrow(location_coords), "location coordinates\n")
  # Print each for verification
  for (i in seq_len(nrow(location_coords))) {
    cat(" ", location_coords$location_id[i], "-> LON:", round(location_coords$lon[i],1),
        "LAT:", round(location_coords$lat[i],1), "\n")
  }
  
  cat("Successfully processed", nrow(location_coords), "location coordinates\n")
  
  plot_data <- eurasia_africa_data %>%
    left_join(location_coords, by = c("location_id", "location_name")) %>%
    filter(!is.na(lon), !is.na(lat))
  
  if (nrow(plot_data) == 0) {
    cat("No plot data after joining with coordinates.\n")
    return(NULL)
  }
  
  cat("Successfully mapped", length(unique(plot_data$location_id)), "locations for plotting\n")
  
  # Use Paul Tol palette (consistent with Brazil maps)
  all_populations <- sort(unique(plot_data$population))
  pop_numbers <- as.integer(str_extract(all_populations, "[0-9]+"))
  all_populations <- all_populations[order(pop_numbers)]
  colors <- paul_tol_palette[all_populations]
  
  # Convert population to factor with numeric ordering BEFORE plotting
  plot_data <- plot_data %>%
    mutate(population = factor(population, levels = all_populations))
  
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "lightgray", color = "white", linewidth = 0.2) +
    coord_sf(xlim = c(-20, 150), ylim = c(-40, 70)) +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(t = 10, b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
    ) +
    labs(
      title = "Population Structure: Europe, Asia & Africa",
      subtitle = paste(dataset_name, "- K =", k_value, "| Samples with >=90% assignment confidence")
    )
  
  unique_locations_plot <- unique(plot_data$location_id)
  cat("Adding", length(unique_locations_plot), "bar charts to map\n")
  
  fixed_width <- 3.0
  fixed_height <- 4.0
  
  for (location in unique_locations_plot) {
    location_df <- plot_data %>% filter(location_id == location)
    
    if (nrow(location_df) == 0) next
    
    location_df <- location_df %>% arrange(population)
    
    location_bar_chart <- ggplot(location_df, aes(x = 1, y = proportion, fill = population)) +
      geom_col(width = 1, color = "black", linewidth = 0.3) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6)
      )
    
    centroid_lon <- unique(location_df$lon)
    centroid_lat <- unique(location_df$lat)
    
    p_map <- p_map +
      annotation_custom(
        grob = ggplotGrob(location_bar_chart),
        xmin = centroid_lon - fixed_width/2,
        xmax = centroid_lon + fixed_width/2,
        ymin = centroid_lat - fixed_height/2,
        ymax = centroid_lat + fixed_height/2
      )
  }
  
  # Legend using Paul Tol palette
  dummy_data <- tibble(population = all_populations, y = 1)
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(
      values = colors,
      name = "Population",
      breaks = all_populations,
      limits = all_populations,
      drop = FALSE
    ) +
    guides(fill = guide_legend(override.aes = list(size = 1))) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.margin = margin(10, 10, 10, 10)
    )
  
  # Extract legend using ggplotGrob to preserve factor order (consistent with Brazil maps)
  legend_grob <- ggplotGrob(legend_plot)
  legend_index <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  legend <- legend_grob$grobs[[legend_index]]
  
  final_plot <- plot_grid(
    p_map,
    legend,
    ncol = 2,
    rel_widths = c(0.8, 0.2),
    align = "h"
  ) +
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  return(final_plot)
}


# 6.3 Run Eurasia map generation for all datasets
cat("CREATING EUROPE/ASIA/AFRICA POPULATION STRUCTURE MAPS")

eurasia_plots_created <- 0
eurasia_failed <- c()

for (i in 1:length(pop_data_list)) {
  file_name <- names(pop_data_list)[i]
  pop_data <- pop_data_list[[i]]
  
  cat("\n", rep("=", 60), "\n")
  cat("PROCESSING EURASIA DATASET", i, ":", file_name, "\n")
  cat(rep("=", 60), "\n")
  
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("SKIPPING", file_name, "- no data loaded\n")
    next
  }
  
  k_value <- case_when(
    str_detect(file_name, "k13") ~ "13",
    str_detect(file_name, "k12") ~ "12",
    str_detect(file_name, "k8")  ~ "8",
    TRUE ~ "unknown"
  )
  
  dataset_type <- case_when(
    str_detect(file_name, "Americas_only") ~ "Americas Only",
    str_detect(file_name, "infantum_only") ~ "L. infantum Only",
    str_detect(file_name, "allsamples")    ~ "All Samples",
    TRUE ~ "Dataset"
  )
  
  location_pop_data <- process_population_data_worldwide(pop_data, dataset_type, samples_data)
  
  if (is.null(location_pop_data)) {
    cat("FAILED: No worldwide data for", dataset_type, "\n")
    eurasia_failed <- c(eurasia_failed, file_name)
    next
  }
  
  eurasia_map_plot <- create_europe_asia_focused_map(location_pop_data, dataset_type, k_value)
  
  if (!is.null(eurasia_map_plot)) {
    print(eurasia_map_plot)
    safe_name <- str_replace_all(tolower(dataset_type), "[^a-z0-9]+", "_")
    filename <- paste0("eurasia_africa_population_map_", safe_name, "_k", k_value, ".jpg")
    filepath <- file.path(FIGURES_DIR, filename)
    ggsave(filepath, eurasia_map_plot, width = 16, height = 12, dpi = 300)
    cat("SUCCESS: Saved", filepath, "\n")
    eurasia_plots_created <- eurasia_plots_created + 1
  } else {
    cat("FAILED: Could not create eurasia map for", dataset_type, "\n")
    eurasia_failed <- c(eurasia_failed, file_name)
  }
  
  if (interactive()) readline(prompt = "Press [Enter] to continue...")
}

cat("EURASIA MAPS SUMMARY")
cat("Total eurasia maps created:", eurasia_plots_created, "\n")
cat("Failed:", length(eurasia_failed), "\n")




# 7. EUROPEAN ORIGIN MAP FOR BRAZILIAN SAMPLES ---------------------------------
# Shows LI subset (infantum_only_k13) only.
# Two-panel layout: Europe (left) + South America (right), shared legend.
# Answers the question: which European populations are ancestral to Brazilian samples?

# 7.1 Create Europe + South America focused map --------------------------------

create_europe_southamerica_map <- function(location_pop_data, dataset_name, k_value) {
  
  # Countries to include in each panel
  european_countries <- c(
    "Kingdom of Spain", "Spain", "Portuguese Republic", "Portugal",
    "French Republic", "France", "Italian Republic", "Italy",
    "Hellenic Republic", "Greece", "Republic of Albania", "Albania",
    "Republic of Serbia", "Yugoslavia", "Serbia", "Georgia", "Syria",
    "State of Israel", "Israel", "Republic of Turkey", "Turkey", "Turkiye",
    "Tunisian Republic", "Tunisia", "Kingdom of Morocco", "Morocco",
    "Palestine", "Iran"
  )
  
  southamerican_countries <- c(
    "Federative Republic of Brazil", "Brazil",
    "Republic of Paraguay", "Paraguay",
    "Oriental Republic of Uruguay", "Uruguay"
  )
  
  # Filter to Europe and South America only
  europe_data <- location_pop_data %>%
    filter(source_country %in% european_countries)
  
  southam_data <- location_pop_data %>%
    filter(source_country %in% southamerican_countries |
             str_starts(location_id, "BR-"))
  
  cat("European locations found:", length(unique(europe_data$location_id)), "\n")
  cat("South American locations found:", length(unique(southam_data$location_id)), "\n")
  
  if (nrow(europe_data) == 0 && nrow(southam_data) == 0) {
    cat("No data for either region — skipping.\n")
    return(NULL)
  }
  
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Coordinate lookup (reuse same short_name_lookup pattern as Eurasia map)
  short_name_lookup <- tribble(
    ~long_name,                                          ~short_name,
    "Tunisian Republic",                                 "Tunisia",
    "State of Israel",                                   "Israel",
    "French Republic",                                   "France",
    "Kingdom of Spain",                                  "Spain",
    "Italian Republic",                                  "Italy",
    "Kingdom of Morocco",                                "Morocco",
    "Portuguese Republic",                               "Portugal",
    "Hellenic Republic",                                 "Greece",
    "Republic of Albania",                               "Albania",
    "Republic of Turkey",                                "Turkey",
    "Republic of Serbia",                                "Serbia",
    "Georgia",                                           "Georgia",
    "Syria",                                             "Syria",
    "Palestine",                                         "Palestine",
    "Iran",                                              "Iran",
    "Federative Republic of Brazil",                     "Brazil",
    "Republic of Paraguay",                              "Paraguay",
    "Oriental Republic of Uruguay",                      "Uruguay"
  )
  
  manual_coords <- tribble(
    ~location_id,                    ~lon,    ~lat,
    "French Republic",                2.2,    46.2,
    "Portuguese Republic",           -8.0,    39.6,
    "Palestine",                     34.5,    31.2,
    "State of Israel",               35.8,    32.5
  )
  
  world_centroids <- world_map %>%
    st_centroid() %>%
    mutate(lon = st_coordinates(.)[,1], lat = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    select(sovereignt, name, lon, lat)
  
  # Helper: resolve coordinates for a data subset
  resolve_coords <- function(df) {
    df %>%
      group_by(location_id, location_name, source_country) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(short_name_lookup, by = c("source_country" = "long_name")) %>%
      mutate(lookup_name = coalesce(short_name, source_country)) %>%
      left_join(world_centroids, by = c("lookup_name" = "name")) %>%
      left_join(manual_coords, by = "location_id", suffix = c("", "_manual")) %>%
      mutate(
        lon = coalesce(lon_manual, lon),
        lat = coalesce(lat_manual, lat)
      ) %>%
      select(-lon_manual, -lat_manual) %>%
      filter(!is.na(lon)) %>%
      select(location_id, location_name, lon, lat)
  }
  
  # For South America, use Brazilian state centroids where available
  brazil_sf    <- ne_download(scale = 50, type = "states", category = "cultural", returnclass = "sf")
  brazil_states <- brazil_sf[brazil_sf$admin == "Brazil", ] %>%
    mutate(code = case_when(
      name == "Acre" ~ "AC", name == "Alagoas" ~ "AL", name == "Amapá" ~ "AP",
      name == "Amazonas" ~ "AM", name == "Bahia" ~ "BA", name == "Ceará" ~ "CE",
      name == "Distrito Federal" ~ "DF", name == "Espírito Santo" ~ "ES",
      name == "Goiás" ~ "GO", name == "Maranhão" ~ "MA", name == "Mato Grosso" ~ "MT",
      name == "Mato Grosso do Sul" ~ "MS", name == "Minas Gerais" ~ "MG",
      name == "Pará" ~ "PA", name == "Paraíba" ~ "PB", name == "Paraná" ~ "PR",
      name == "Pernambuco" ~ "PE", name == "Piauí" ~ "PI", name == "Rio de Janeiro" ~ "RJ",
      name == "Rio Grande do Norte" ~ "RN", name == "Rio Grande do Sul" ~ "RS",
      name == "Rondônia" ~ "RO", name == "Roraima" ~ "RR", name == "Santa Catarina" ~ "SC",
      name == "São Paulo" ~ "SP", name == "Sergipe" ~ "SE", name == "Tocantins" ~ "TO",
      TRUE ~ as.character(abbrev)
    )) %>%
    select(name, code, geometry) %>%
    st_make_valid()
  
  state_centroids <- brazil_states %>%
    st_centroid() %>%
    st_coordinates() %>%
    as_tibble() %>%
    rename(lon = X, lat = Y) %>%
    bind_cols(brazil_states %>% st_drop_geometry()) %>%
    select(code, lon, lat)
  
  # Resolve Europe coordinates
  europe_coords <- if (nrow(europe_data) > 0) resolve_coords(europe_data) else NULL
  
  # Resolve South America coordinates (state-level for Brazil, country-level otherwise)
  southam_state_data <- southam_data %>% filter(str_starts(location_id, "BR-"))
  southam_country_data <- southam_data %>% filter(!str_starts(location_id, "BR-"))
  
  southam_state_coords <- if (nrow(southam_state_data) > 0) {
    state_centroids_named <- brazil_states %>%
      st_drop_geometry() %>%
      bind_cols(
        brazil_states %>% st_centroid() %>% st_coordinates() %>%
          as_tibble() %>% rename(lon = X, lat = Y)
      ) %>%
      select(name, lon, lat)
    southam_state_data %>%
      mutate(state_name = str_remove(location_id, "BR-")) %>%
      group_by(location_id, location_name, state_name) %>%
      slice(1) %>%
      ungroup() %>%
      left_join(state_centroids_named, by = c("state_name" = "name")) %>%
      filter(!is.na(lon)) %>%
      select(location_id, location_name, lon, lat)
  } else NULL
  
  southam_country_coords <- if (nrow(southam_country_data) > 0) {
    resolve_coords(southam_country_data)
  } else NULL
  
  southam_coords <- bind_rows(southam_state_coords, southam_country_coords)
  
  # Join coords back to data
  join_coords <- function(df, coords) {
    if (is.null(coords) || nrow(coords) == 0) return(NULL)
    df %>%
      left_join(coords, by = c("location_id", "location_name")) %>%
      filter(!is.na(lon), !is.na(lat))
  }
  
  europe_plot_data   <- join_coords(europe_data, europe_coords)
  southam_plot_data  <- join_coords(southam_data, southam_coords)
  
  # Colour palette — LI populations only
  all_pops <- sort(unique(location_pop_data$population),
                   decreasing = FALSE)
  pop_nums <- as.integer(str_extract(all_pops, "[0-9]+"))
  all_pops <- all_pops[order(pop_nums)]
  colors   <- paul_tol_palette[all_pops]
  
  for (df in list(europe_plot_data, southam_plot_data)) {
    if (!is.null(df)) {
      df <- df %>% mutate(population = factor(population, levels = all_pops))
    }
  }
  if (!is.null(europe_plot_data))
    europe_plot_data <- europe_plot_data %>%
    mutate(population = factor(population, levels = all_pops))
  if (!is.null(southam_plot_data))
    southam_plot_data <- southam_plot_data %>%
    mutate(population = factor(population, levels = all_pops))
  
  # Helper: build one panel map
  build_panel <- function(plot_df, xlim, ylim, title, fixed_width, fixed_height) {
    if (is.null(plot_df) || nrow(plot_df) == 0) {
      return(ggplot() + theme_void() +
               labs(title = title) +
               theme(plot.title = element_text(hjust = 0.5, size = 13, face = "bold")))
    }
    
    p <- ggplot() +
      geom_sf(data = world_map, fill = "lightgray", color = "white", linewidth = 0.2) +
      coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
      theme_minimal() +
      theme(
        panel.grid.major  = element_blank(),
        panel.grid.minor  = element_blank(),
        axis.title        = element_blank(),
        axis.text         = element_blank(),
        axis.ticks        = element_blank(),
        plot.title        = element_text(hjust = 0.5, size = 13, face = "bold",
                                         margin = margin(t = 8, b = 4)),
        legend.position   = "none",
        plot.margin       = margin(t = 10, r = 8, b = 8, l = 8)
      ) +
      labs(title = title)
    
    for (loc in unique(plot_df$location_id)) {
      loc_df <- plot_df %>% filter(location_id == loc) %>% arrange(population)
      if (nrow(loc_df) == 0) next
      bar <- ggplot(loc_df, aes(x = 1, y = proportion, fill = population)) +
        geom_col(width = 1, color = "black", linewidth = 0.3) +
        scale_fill_manual(values = colors, drop = FALSE) +
        scale_y_continuous(expand = c(0, 0)) +
        theme_void() +
        theme(
          legend.position    = "none",
          plot.margin        = unit(c(0, 0, 0, 0), "cm"),
          panel.background   = element_rect(fill = "white", color = "black", linewidth = 0.6)
        )
      clon <- unique(loc_df$lon)
      clat <- unique(loc_df$lat)
      p <- p + annotation_custom(
        grob = ggplotGrob(bar),
        xmin = clon - fixed_width / 2, xmax = clon + fixed_width / 2,
        ymin = clat - fixed_height / 2, ymax = clat + fixed_height / 2
      )
    }
    p
  }
  
  # Build panels — South America left, Europe right, uniform bar sizes
  panel_southam <- build_panel(
    southam_plot_data,
    xlim = c(-82, -34), ylim = c(-56, 15),
    title = "South America",
    fixed_width = 1.6, fixed_height = 2.2
  )
  
  panel_europe  <- build_panel(
    europe_plot_data,
    xlim = c(-12, 45), ylim = c(22, 72),
    title = "Europe & Middle East",
    fixed_width = 1.6, fixed_height = 2.2
  )
  
  # Shared legend
  dummy_data  <- tibble(population = factor(all_pops, levels = all_pops), y = 1)
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Population",
                      breaks = all_pops, limits = all_pops, drop = FALSE) +
    guides(fill = guide_legend(override.aes = list(size = 1))) +
    theme_void() +
    theme(
      legend.position    = "right",
      legend.key.size    = unit(0.6, "cm"),
      legend.text        = element_text(size = 11),
      legend.title       = element_text(size = 12, face = "bold"),
      legend.margin      = margin(10, 10, 10, 10)
    )
  
  legend_grob  <- ggplotGrob(legend_plot)
  legend_index <- which(sapply(legend_grob$grobs, function(x) x$name) == "guide-box")
  legend       <- legend_grob$grobs[[legend_index]]
  
  # Title row
  title_text <- paste0(
    "European Origins of Brazilian L. infantum\n",
    dataset_name, " - K = ", k_value,
    " | Samples with \u226590% assignment confidence"
  )
  
  title_grob <- ggdraw() +
    draw_label(title_text, fontface = "bold", size = 14, hjust = 0.5)
  
  two_panels <- plot_grid(panel_southam, panel_europe,
                          ncol = 2, rel_widths = c(0.85, 1), align = "h")
  
  final_plot <- plot_grid(
    title_grob,
    plot_grid(two_panels, legend, ncol = 2, rel_widths = c(0.82, 0.18), align = "h"),
    ncol = 1, rel_heights = c(0.08, 0.92)
  )
  
  return(final_plot)
}


# 7.2 Run European origin map — LI subset only ---------------------------------

cat("\n", rep("=", 60), "\n")
cat("CREATING EUROPEAN ORIGIN MAP FOR BRAZILIAN SAMPLES\n")
cat("Dataset: infantum_only_k13 (LI subset)\n")
cat(rep("=", 60), "\n")

li_pop_data <- pop_data_list[["infantum_only_k13"]]

if (!is.null(li_pop_data) && nrow(li_pop_data) > 0) {
  
  li_worldwide <- process_population_data_worldwide(li_pop_data, "L. infantum Only", samples_data)
  
  if (!is.null(li_worldwide)) {
    euro_southam_plot <- create_europe_southamerica_map(li_worldwide, "L. infantum Only", "13")
    
    if (!is.null(euro_southam_plot)) {
      print(euro_southam_plot)
      filepath <- file.path(FIGURES_DIR, "european_origin_brazil_l_infantum_only_k13.jpg")
      ggsave(filepath, euro_southam_plot, width = 18, height = 11, dpi = 300)
      cat("SUCCESS: Saved", filepath, "\n")
    } else {
      cat("FAILED: Could not create European origin map.\n")
    }
  } else {
    cat("FAILED: No worldwide data returned for LI subset.\n")
  }
  
} else {
  cat("FAILED: infantum_only_k13 data not loaded.\n")
}