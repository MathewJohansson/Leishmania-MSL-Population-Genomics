# WORLDWIDE POPULATION STRUCTURE MAP - UPGRADED VERSION
# Modified from your original Brazil-focused code
# IMPROVEMENTS: Consistent population naming (v1, v2, v3...) and colorblind-friendly colors

# Load additional libraries for worldwide mapping
library(tidyverse)    
library(ggplot2)      
library(RColorBrewer) 
library(gridExtra)    
library(sf)           
library(rnaturalearth) 
library(rnaturalearthdata) 
library(cowplot)      
library(readr)
library(viridis)      # For colorblind-friendly colors

# Create UPGRADE directory if it doesn't exist
if (!dir.exists("UPGRADE")) {
  dir.create("UPGRADE")
  cat("Created UPGRADE directory\n")
}

# --- TRUE COLORBLIND-FRIENDLY PALETTE FUNCTION (AVOIDS RED-GREEN) ---
get_colorblind_palette <- function(n) {
  if (n <= 8) {
    # Paul Tol's colorblind-safe palette - avoids red-green confusion
    colors <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", 
                "#66CCEE", "#AA3377", "#BBBBBB", "#000000")
    return(colors[1:n])
  } else if (n <= 12) {
    # Extended colorblind-safe palette
    colors <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", 
                "#66CCEE", "#AA3377", "#BBBBBB", "#332288",
                "#88CCEE", "#44AA99", "#117733", "#999933")
    return(colors[1:n])
  } else {
    # For many populations, use viridis plasma (purple-pink-yellow, no red-green)
    return(viridis::plasma(n))
  }
}

# --- POPULATION STANDARDIZATION FUNCTION (CONSISTENT ACROSS ALL SCRIPTS) ---
standardize_population_names <- function(pop_data) {
  # Get unique populations and sort them consistently
  unique_pops <- sort(unique(pop_data$population))
  
  # Create mapping from original names to v1, v2, v3, etc.
  pop_mapping <- tibble(
    original_population = unique_pops,
    standardized_population = paste0("v", 1:length(unique_pops))
  )
  
  cat("Population mapping:\n")
  for (i in 1:nrow(pop_mapping)) {
    cat("  ", pop_mapping$original_population[i], " -> ", pop_mapping$standardized_population[i], "\n")
  }
  
  # Apply mapping to data
  standardized_data <- pop_data %>%
    left_join(pop_mapping, by = c("population" = "original_population")) %>%
    select(-population) %>%
    rename(population = standardized_population)
  
  return(list(data = standardized_data, mapping = pop_mapping))
}

# --- WORLDWIDE VERSION: Process population data globally with standardization ---

process_population_data_worldwide <- function(pop_data, dataset_name, samples_data) {
  cat("\n=== Processing WORLDWIDE", dataset_name, "===\n")
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("No population data provided for processing (", dataset_name, ").\n")
    return(NULL)
  }
  
  # Filter for high confidence assignment (>0.8)
  pop_filtered <- pop_data %>% filter(population_proportion > 0.8)
  cat("Samples with >0.8 population proportion:", nrow(pop_filtered), "\n")
  
  if (nrow(pop_filtered) == 0) {
    cat("No samples with >0.8 population proportion found after filtering for", dataset_name, ".\n")
    return(NULL)
  }
  
  # STANDARDIZE POPULATION NAMES TO v1, v2, v3, etc.
  standardized_result <- standardize_population_names(pop_filtered)
  pop_filtered <- standardized_result$data
  pop_mapping <- standardized_result$mapping
  
  # Clean samples data with your existing country fixes
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
  
  # WORLDWIDE VERSION: Don't filter for Brazil only
  combined_data_worldwide <- combined_data %>%
    filter(!is.na(source_country)) %>%
    # Create a location identifier combining country and state/region
    mutate(
      location_id = case_when(
        # For Brazil, use state codes
        (source_country == "Federative Republic of Brazil" | source_country == "Brazil") & 
          !is.na(source_state_or_region) ~ paste0("BR-", source_state_or_region),
        # For other countries, use country name
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
  location_populations <- combined_data_worldwide %>%
    group_by(location_id, location_name, population, source_country, source_state_or_region) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(location_id) %>%
    mutate(
      total_samples = sum(count),
      proportion = count / total_samples
    ) %>%
    ungroup() %>%
    filter(total_samples >= 1) %>%  # Include locations with 1+ samples to show Europe/Asia
    arrange(desc(total_samples))
  
  cat("Locations with sufficient data:", length(unique(location_populations$location_id)), "\n")
  cat("Countries represented:", length(unique(location_populations$source_country)), "\n")
  
  # Print detailed breakdown by region for debugging
  cat("\n=== LOCATION BREAKDOWN ===\n")
  location_summary <- location_populations %>%
    group_by(location_id, location_name, source_country) %>%
    summarise(total_samples = first(total_samples), .groups = 'drop') %>%
    arrange(desc(total_samples))
  
  for (i in 1:min(20, nrow(location_summary))) {
    cat(sprintf("  %s (%s): %d samples\n", 
                location_summary$location_name[i], 
                location_summary$source_country[i], 
                location_summary$total_samples[i]))
  }
  
  return(list(data = location_populations, mapping = pop_mapping))
}

# --- UPGRADED WORLDWIDE VERSION: Create map function ---

create_worldwide_map_with_bar_charts_UPGRADED <- function(processed_result, dataset_name, k_value) {
  if (is.null(processed_result) || is.null(processed_result$data) || nrow(processed_result$data) == 0) {
    cat("DEBUG: processed_result is NULL or empty.\n")
    return(NULL)
  }
  
  location_pop_data <- processed_result$data
  pop_mapping <- processed_result$mapping
  
  cat("Creating UPGRADED worldwide map for", dataset_name, "with", nrow(location_pop_data), "data points\n")
  
  # Get world map data
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Initialize location coordinates tibble with proper structure
  location_coords <- tibble(
    location_id = character(0),
    location_name = character(0),
    lon = numeric(0),
    lat = numeric(0)
  )
  
  # Process each unique location
  unique_locations <- unique(location_pop_data$location_id)
  cat("Processing coordinates for", length(unique_locations), "unique locations\n")
  
  for (location in unique_locations) {
    location_data <- location_pop_data %>% 
      filter(location_id == location) %>% 
      slice(1)  # Get first row for location info
    
    cat("Processing location:", location, "\n")
    
    if (str_starts(location, "BR-")) {
      # Brazil state - use your existing state mapping
      state_code <- str_remove(location, "BR-")
      cat("  Brazil state code:", state_code, "\n")
      
      # Check if brazil_states exists
      if (exists("brazil_states")) {
        brazil_state <- brazil_states %>% filter(code == state_code)
        
        if (nrow(brazil_state) > 0) {
          centroid <- brazil_state %>%
            st_centroid() %>%
            st_coordinates()
          
          location_coords <- bind_rows(location_coords, tibble(
            location_id = location,
            location_name = location_data$location_name,
            lon = as.numeric(centroid[1]),
            lat = as.numeric(centroid[2])
          ))
          cat("  Added Brazil state coordinates\n")
        } else {
          cat("  Warning: Brazil state not found in map data\n")
        }
      } else {
        cat("  Warning: brazil_states object not found\n")
      }
    } else {
      # Other countries
      cat("  Country:", location_data$source_country, "\n")
      country_map <- world_map %>% 
        filter(name == location_data$source_country | 
                 name_long == location_data$source_country |
                 sovereignt == location_data$source_country |
                 admin == location_data$source_country)
      
      if (nrow(country_map) > 0) {
        centroid <- country_map %>%
          st_centroid() %>%
          st_coordinates()
        
        location_coords <- bind_rows(location_coords, tibble(
          location_id = location,
          location_name = location_data$location_name,
          lon = as.numeric(centroid[1]),
          lat = as.numeric(centroid[2])
        ))
        cat("  Added country coordinates\n")
      } else {
        cat("  Warning: Country not found in world map data\n")
        # Try some common alternative country name patterns
        alt_patterns <- c(
          paste0("United States", ".*"),
          paste0(".*", str_remove_all(location_data$source_country, "Republic of |Kingdom of |State of "), ".*")
        )
        
        for (pattern in alt_patterns) {
          country_map_alt <- world_map %>% 
            filter(str_detect(name, regex(pattern, ignore_case = TRUE)) |
                     str_detect(name_long, regex(pattern, ignore_case = TRUE)) |
                     str_detect(admin, regex(pattern, ignore_case = TRUE)))
          
          if (nrow(country_map_alt) > 0) {
            centroid <- country_map_alt %>%
              slice(1) %>%  # Take first match
              st_centroid() %>%
              st_coordinates()
            
            location_coords <- bind_rows(location_coords, tibble(
              location_id = location,
              location_name = location_data$location_name,
              lon = as.numeric(centroid[1]),
              lat = as.numeric(centroid[2])
            ))
            cat("  Added country coordinates using alternative matching\n")
            break
          }
        }
      }
    }
  }
  
  cat("Successfully processed", nrow(location_coords), "location coordinates\n")
  
  # Merge population data with coordinates
  plot_data <- location_pop_data %>%
    left_join(location_coords, by = c("location_id", "location_name")) %>%
    filter(!is.na(lon), !is.na(lat))
  
  if (nrow(plot_data) == 0) {
    cat("No plot data after joining with coordinates.\n")
    return(NULL)
  }
  
  cat("Successfully mapped", length(unique(plot_data$location_id)), "locations\n")
  
  # Create COLORBLIND-FRIENDLY palette with consistent ordering
  all_populations <- sort(unique(plot_data$population))  # v1, v2, v3, etc. will sort correctly
  n_pops <- length(all_populations)
  
  colors <- get_colorblind_palette(n_pops)
  names(colors) <- all_populations
  
  cat("Using", n_pops, "colorblind-friendly colors for populations:", paste(all_populations, collapse = ", "), "\n")
  
  # Base world map with improved styling
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "lightgray", color = "white", linewidth = 0.2) +
    coord_sf(crs = st_crs(4326)) +  # Use WGS84 projection
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(t = 15, b = 8)),
      plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
      legend.position = "none",
      plot.margin = margin(t = 25, r = 15, b = 15, l = 15)  
    ) +
    labs(
      title = "Worldwide Population Structure",
      subtitle = paste(dataset_name, "- K =", k_value, "| Standardized populations (v1, v2, v3...) | >80% assignment confidence")
    )
  
  # Add bar charts for each location
  unique_locations <- unique(plot_data$location_id)
  cat("Processing", length(unique_locations), "locations for bar charts\n")
  
  # Adjust bar chart size for worldwide view
  fixed_width <- 2.8   # Slightly larger for worldwide view
  fixed_height <- 3.5  
  
  for (location in unique_locations) {
    location_df <- plot_data %>% filter(location_id == location)
    
    if (nrow(location_df) == 0) next
    
    # Ensure consistent ordering by population name (v1, v2, v3, etc.)
    location_df <- location_df %>% arrange(population)
    total_samples <- unique(location_df$total_samples)
    
    # Create bar chart with colorblind-friendly colors
    location_bar_chart <- ggplot(location_df, aes(x = 1, y = proportion, fill = population)) +
      geom_col(width = 1, color = "black", linewidth = 0.25) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.6)
      )
    
    # Get coordinates
    centroid_lon <- unique(location_df$lon)
    centroid_lat <- unique(location_df$lat)
    
    # Add bar chart to map
    p_map <- p_map +
      annotation_custom(
        grob = ggplotGrob(location_bar_chart),
        xmin = centroid_lon - fixed_width/2,
        xmax = centroid_lon + fixed_width/2,
        ymin = centroid_lat - fixed_height/2,
        ymax = centroid_lat + fixed_height/2
      )
  }
  
  # Create legend with improved styling
  dummy_data <- tibble(
    population = all_populations,
    y = 1
  )
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Population\n(Standardized)") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.7, "cm"),
      legend.text = element_text(size = 12, family = "sans"),
      legend.title = element_text(size = 13, face = "bold", family = "sans"),
      legend.margin = margin(15, 15, 15, 15)
    )
  
  # Extract legend
  legend <- get_legend(legend_plot)
  
  # Combine map and legend
  final_plot <- plot_grid(
    p_map, 
    legend, 
    ncol = 2, 
    rel_widths = c(0.8, 0.2),  
    align = "h"
  ) +
    theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15))
  
  # Save population mapping information
  mapping_info <- list(
    dataset_name = dataset_name,
    k_value = k_value,
    population_mapping = pop_mapping,
    total_populations = n_pops,
    colors_used = colors,
    map_type = "worldwide"
  )
  
  return(list(plot = final_plot, mapping_info = mapping_info))
}

# --- UPGRADED EUROPE/ASIA FOCUSED VERSION ---

create_europe_asia_focused_map_UPGRADED <- function(processed_result, dataset_name, k_value) {
  if (is.null(processed_result) || is.null(processed_result$data) || nrow(processed_result$data) == 0) {
    cat("DEBUG: processed_result is NULL or empty.\n")
    return(NULL)
  }
  
  location_pop_data <- processed_result$data
  pop_mapping <- processed_result$mapping
  
  # Define European and Asian countries for better matching
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
  
  # Filter for Europe, Asia, and Africa (excluding Americas)
  eurasia_africa_data <- location_pop_data %>%
    filter(source_country %in% c(european_countries, asian_countries, african_countries) |
             str_detect(source_country, "Spain|Portugal|France|Italy|Greece|Albania|Serbia|China|Turkey|Israel|Yemen|Uzbekistan|Egypt|Morocco|Tunisia|Ethiopia|Sudan"))
  
  cat("Europe/Asia/Africa samples found:", nrow(eurasia_africa_data), "\n")
  cat("Unique locations:", length(unique(eurasia_africa_data$location_id)), "\n")
  
  if (nrow(eurasia_africa_data) == 0) {
    cat("No Europe/Asia/Africa data found\n")
    return(NULL)
  }
  
  # Get world map data
  world_map <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Initialize location coordinates tibble with proper structure
  location_coords <- tibble(
    location_id = character(0),
    location_name = character(0),
    lon = numeric(0),
    lat = numeric(0)
  )
  
  # Process each unique location in Europe/Asia/Africa
  unique_locations <- unique(eurasia_africa_data$location_id)
  cat("Processing coordinates for", length(unique_locations), "unique locations\n")
  
  for (location in unique_locations) {
    location_data <- eurasia_africa_data %>% 
      filter(location_id == location) %>% 
      slice(1)
    
    cat("Processing location:", location, "Country:", location_data$source_country, "\n")
    
    # For Europe/Asia, all should be country-level (no Brazil states)
    country_map <- world_map %>% 
      filter(name == location_data$source_country | 
               name_long == location_data$source_country |
               sovereignt == location_data$source_country |
               admin == location_data$source_country)
    
    if (nrow(country_map) > 0) {
      centroid <- country_map %>%
        st_centroid() %>%
        st_coordinates()
      
      location_coords <- bind_rows(location_coords, tibble(
        location_id = location,
        location_name = location_data$location_name,
        lon = as.numeric(centroid[1]),
        lat = as.numeric(centroid[2])
      ))
      cat("  Added coordinates\n")
    } else {
      cat("  Warning: Country not found, trying alternative matching\n")
      # Try alternative matching patterns
      alt_patterns <- c(
        str_remove_all(location_data$source_country, "Republic of |Kingdom of |State of |People's Republic of |Arab Republic of |Federal Democratic Republic of |Tunisian |Oriental "),
        str_extract(location_data$source_country, "\\b[A-Z][a-z]+\\b")
      )
      
      for (pattern in alt_patterns) {
        if (is.na(pattern) || pattern == "") next
        
        country_map_alt <- world_map %>% 
          filter(str_detect(name, regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE)) |
                   str_detect(name_long, regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE)) |
                   str_detect(admin, regex(paste0("\\b", pattern, "\\b"), ignore_case = TRUE)))
        
        if (nrow(country_map_alt) > 0) {
          centroid <- country_map_alt %>%
            slice(1) %>%
            st_centroid() %>%
            st_coordinates()
          
          location_coords <- bind_rows(location_coords, tibble(
            location_id = location,
            location_name = location_data$location_name,
            lon = as.numeric(centroid[1]),
            lat = as.numeric(centroid[2])
          ))
          cat("  Added coordinates using pattern:", pattern, "\n")
          break
        }
      }
    }
  }
  
  cat("Successfully processed", nrow(location_coords), "location coordinates\n")
  
  # Merge population data with coordinates
  plot_data <- eurasia_africa_data %>%
    left_join(location_coords, by = c("location_id", "location_name")) %>%
    filter(!is.na(lon), !is.na(lat))
  
  if (nrow(plot_data) == 0) {
    cat("No plot data after joining with coordinates.\n")
    return(NULL)
  }
  
  cat("Successfully mapped", length(unique(plot_data$location_id)), "locations for plotting\n")
  
  # Create COLORBLIND-FRIENDLY palette with consistent ordering
  all_populations <- sort(unique(plot_data$population))  # v1, v2, v3, etc. will sort correctly
  n_pops <- length(all_populations)
  
  colors <- get_colorblind_palette(n_pops)
  names(colors) <- all_populations
  
  cat("Using", n_pops, "colorblind-friendly colors for populations:", paste(all_populations, collapse = ", "), "\n")
  
  # Base world map focused on Europe/Asia/Africa with improved styling
  p_map <- ggplot() +
    geom_sf(data = world_map, fill = "lightgray", color = "white", linewidth = 0.2) +
    coord_sf(xlim = c(-20, 150), ylim = c(-40, 70)) +  # Focus on Europe/Asia/Africa
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold", margin = margin(t = 15, b = 8)),
      plot.subtitle = element_text(hjust = 0.5, size = 14, margin = margin(b = 15)),
      legend.position = "none",
      plot.margin = margin(t = 25, r = 15, b = 15, l = 15)  
    ) +
    labs(
      title = "Population Structure: Europe, Asia & Africa",
      subtitle = paste(dataset_name, "- K =", k_value, "| Standardized populations (v1, v2, v3...) | >80% assignment confidence")
    )
  
  # Add bar charts for each location
  unique_locations_plot <- unique(plot_data$location_id)
  cat("Adding", length(unique_locations_plot), "bar charts to map\n")
  
  # Bar chart size
  fixed_width <- 3.5   
  fixed_height <- 4.5  
  
  for (location in unique_locations_plot) {
    location_df <- plot_data %>% filter(location_id == location)
    
    if (nrow(location_df) == 0) next
    
    # Ensure consistent ordering by population name (v1, v2, v3, etc.)
    location_df <- location_df %>% arrange(population)
    total_samples <- unique(location_df$total_samples)
    
    # Create bar chart with colorblind-friendly colors
    location_bar_chart <- ggplot(location_df, aes(x = 1, y = proportion, fill = population)) +
      geom_col(width = 1, color = "black", linewidth = 0.3) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7)
      )
    
    # Get coordinates
    centroid_lon <- unique(location_df$lon)
    centroid_lat <- unique(location_df$lat)
    
    # Add bar chart to map
    p_map <- p_map +
      annotation_custom(
        grob = ggplotGrob(location_bar_chart),
        xmin = centroid_lon - fixed_width/2,
        xmax = centroid_lon + fixed_width/2,
        ymin = centroid_lat - fixed_height/2,
        ymax = centroid_lat + fixed_height/2
      )
  }
  
  # Create legend with improved styling
  dummy_data <- tibble(
    population = all_populations,
    y = 1
  )
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Population\n(Standardized)") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.7, "cm"),
      legend.text = element_text(size = 12, family = "sans"),
      legend.title = element_text(size = 13, face = "bold", family = "sans"),
      legend.margin = margin(15, 15, 15, 15)
    )
  
  # Extract legend
  legend <- get_legend(legend_plot)
  
  # Combine map and legend
  final_plot <- plot_grid(
    p_map, 
    legend, 
    ncol = 2, 
    rel_widths = c(0.8, 0.2),  
    align = "h"
  ) +
    theme(plot.margin = margin(t = 15, r = 15, b = 15, l = 15))
  
  # Save population mapping information
  mapping_info <- list(
    dataset_name = dataset_name,
    k_value = k_value,
    population_mapping = pop_mapping,
    total_populations = n_pops,
    colors_used = colors,
    map_type = "europe_asia_africa"
  )
  
  return(list(plot = final_plot, mapping_info = mapping_info))
}

# --- MAIN EXECUTION FOR UPGRADED WORLDWIDE PLOTS ---

cat("=== CREATING UPGRADED WORLDWIDE POPULATION STRUCTURE MAPS ===\n")

# Process each dataset for worldwide analysis
worldwide_plots_created <- 0
eurasia_plots_created <- 0
all_worldwide_mapping_info <- list()

for (i in 1:length(pop_data_list)) {
  file_name <- names(pop_data_list)[i]
  pop_data <- pop_data_list[[i]]
  
  cat("\n", rep("=", 60), "\n")
  cat("PROCESSING UPGRADED WORLDWIDE DATASET", i, ":", file_name, "\n")
  cat(rep("=", 60), "\n")
  
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("SKIPPING", file_name, "- no data loaded\n")
    next
  }
  
  # Extract dataset info (same as your original code)
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
  
  # Process data worldwide with standardization
  processed_result <- process_population_data_worldwide(pop_data, dataset_type, samples_data)
  
  if (is.null(processed_result)) {
    cat("FAILED: No processed worldwide data for", dataset_type, "\n")
    next
  }
  
  # Create UPGRADED worldwide map
  worldwide_map_result <- create_worldwide_map_with_bar_charts_UPGRADED(processed_result, dataset_type, k_value)
  
  if (!is.null(worldwide_map_result) && !is.null(worldwide_map_result$plot)) {
    # Display the plot
    print(worldwide_map_result$plot)
    
    # Save plot to UPGRADE directory
    safe_name <- str_replace_all(tolower(dataset_type), "[^a-z0-9]+", "_") 
    filename <- paste0("UPGRADE/worldwide_population_map_UPGRADED_", safe_name, "_k", k_value, ".png")
    ggsave(filename, worldwide_map_result$plot, width = 18, height = 14, dpi = 300)
    cat("SUCCESS: UPGRADED worldwide plot saved as", filename, "\n")
    
    # Store mapping information
    all_worldwide_mapping_info[[paste0(file_name, "_worldwide")]] <- worldwide_map_result$mapping_info
    
    worldwide_plots_created <- worldwide_plots_created + 1
  }
  
  # Create UPGRADED Europe/Asia/Africa focused map
  eurasia_map_result <- create_europe_asia_focused_map_UPGRADED(processed_result, dataset_type, k_value)
  
  if (!is.null(eurasia_map_result) && !is.null(eurasia_map_result$plot)) {
    # Display the plot
    print(eurasia_map_result$plot)
    
    # Save plot to UPGRADE directory
    filename_eurasia <- paste0("UPGRADE/eurasia_africa_population_map_UPGRADED_", safe_name, "_k", k_value, ".png")
    ggsave(filename_eurasia, eurasia_map_result$plot, width = 18, height = 14, dpi = 300)
    cat("SUCCESS: UPGRADED Europe/Asia/Africa plot saved as", filename_eurasia, "\n")
    
    # Store mapping information
    all_worldwide_mapping_info[[paste0(file_name, "_eurasia_africa")]] <- eurasia_map_result$mapping_info
    
    eurasia_plots_created <- eurasia_plots_created + 1
  }
  
  # Add a pause to see each plot
  if (interactive()) {
    readline(prompt = "Press [Enter] to continue to next dataset...")
  }
}

# --- SAVE WORLDWIDE MAPPING INFORMATION ---

# Save all worldwide mapping information to UPGRADE directory
saveRDS(all_worldwide_mapping_info, "UPGRADE/worldwide_population_mapping_info.rds")
cat("\nWorldwide population mapping information saved to UPGRADE/worldwide_population_mapping_info.rds\n")

# Create a worldwide summary report
worldwide_summary_report <- tibble(
  dataset = names(all_worldwide_mapping_info),
  k_value = map_chr(all_worldwide_mapping_info, ~ .x$k_value),
  total_populations = map_int(all_worldwide_mapping_info, ~ .x$total_populations),
  dataset_name = map_chr(all_worldwide_mapping_info, ~ .x$dataset_name),
  map_type = map_chr(all_worldwide_mapping_info, ~ .x$map_type)
)

write_csv(worldwide_summary_report, "UPGRADE/worldwide_population_analysis_summary.csv")
cat("Worldwide summary report saved to UPGRADE/worldwide_population_analysis_summary.csv\n")

cat("\n=== FINAL UPGRADED WORLDWIDE SUMMARY ===\n")
cat("Total UPGRADED worldwide plots successfully created:", worldwide_plots_created, "\n")
cat("Total UPGRADED Europe/Asia/Africa plots successfully created:", eurasia_plots_created, "\n")
cat("All files saved to UPGRADE/ directory with colorblind-friendly colors and standardized population names!\n")
cat("UPGRADED worldwide analysis completed!\n")

save.image("UPGRADE/worldwide_map_plotting_UPGRADED_data.Rda")


#------------------------------------------------------------------------------------------------------------------------------------------
# AUTHOR:Zeynep Sakaoglu
# DATE: 28.07.2025
# UPGRADED VERSION: Consistent population naming (v1, v2, v3...) and colorblind-friendly colors

# Load required libraries
library(tidyverse)    
library(ggplot2)      
library(RColorBrewer) 
library(gridExtra)    
library(sf)           
library(rnaturalearth) 
library(rnaturalearthdata) 
library(cowplot)      
library(readr)
library(viridis)      # For colorblind-friendly colors

# Create UPGRADE directory if it doesn't exist
if (!dir.exists("UPGRADE")) {
  dir.create("UPGRADE")
  cat("Created UPGRADE directory\n")
}

# --- TRUE COLORBLIND-FRIENDLY PALETTE FUNCTION (AVOIDS RED-GREEN) ---
get_colorblind_palette <- function(n) {
  if (n <= 8) {
    # Paul Tol's colorblind-safe palette - avoids red-green confusion
    colors <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", 
                "#66CCEE", "#AA3377", "#BBBBBB", "#000000")
    return(colors[1:n])
  } else if (n <= 12) {
    # Extended colorblind-safe palette
    colors <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", 
                "#66CCEE", "#AA3377", "#BBBBBB", "#332288",
                "#88CCEE", "#44AA99", "#117733", "#999933")
    return(colors[1:n])
  } else {
    # For many populations, use viridis plasma (purple-pink-yellow, no red-green)
    return(viridis::plasma(n))
  }
}

# --- POPULATION STANDARDIZATION FUNCTION ---
standardize_population_names <- function(pop_data) {
  # Get unique populations and sort them consistently
  unique_pops <- sort(unique(pop_data$population))
  
  # Create mapping from original names to v1, v2, v3, etc.
  pop_mapping <- tibble(
    original_population = unique_pops,
    standardized_population = paste0("v", 1:length(unique_pops))
  )
  
  cat("Population mapping:\n")
  for (i in 1:nrow(pop_mapping)) {
    cat("  ", pop_mapping$original_population[i], " -> ", pop_mapping$standardized_population[i], "\n")
  }
  
  # Apply mapping to data
  standardized_data <- pop_data %>%
    left_join(pop_mapping, by = c("population" = "original_population")) %>%
    select(-population) %>%
    rename(population = standardized_population)
  
  return(list(data = standardized_data, mapping = pop_mapping))
}

# --- 1. Data Loading and Preparation ---

# Country name standardization
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
  "Palestine", "West Bank",
  "Ethiopia", "Federal Democratic Republic of Ethiopia",
  "Sudan", "Republic of the Sudan"
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

# Load data files
samples_data <- read_tsv("data/vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv", show_col_types = FALSE) %>%
  rename(sample = sra_run_accession)

# Load population membership files
pop_membership_Americas_k13 <- read_tsv("data/population_membership_table.Americas_only.k13.tsv", show_col_types = FALSE)
pop_membership_Infantum_k13 <- read_tsv("data/population_membership_table.infantum_only.k13.tsv", show_col_types = FALSE)
pop_membership_All_k12 <- read_tsv("data/population_membership_table.allsamples.k12.tsv", show_col_types = FALSE)
pop_membership_Infantum_k8 <- read_tsv("data/population_membership_table.infantum_only.k8.tsv", show_col_types = FALSE)

# Create a named list with your loaded population dataframes
pop_data_list <- list(
  "Americas_only_k13" = pop_membership_Americas_k13,
  "infantum_only_k13" = pop_membership_Infantum_k13,
  "allsamples_k12" = pop_membership_All_k12,
  "infantum_only_k8" = pop_membership_Infantum_k8
)

#--------------------------PART 2-----------------------------------------------
# --- 2. Helper Functions ---

# Function to process population membership data and merge with sample locations
process_population_data <- function(pop_data, dataset_name, samples_data) {
  cat("\n=== Processing", dataset_name, "===\n")
  if (is.null(pop_data) || nrow(pop_data) == 0) {
    cat("No population data provided for processing (", dataset_name, ").\n")
    return(NULL)
  }
  
  # Filter for high confidence assignment (>0.8)
  pop_filtered <- pop_data %>% filter(population_proportion > 0.8)
  cat("Samples with >0.8 population proportion:", nrow(pop_filtered), "\n")
  
  if (nrow(pop_filtered) == 0) {
    cat("No samples with >0.8 population proportion found after filtering for", dataset_name, ".\n")
    return(NULL)
  }
  
  # STANDARDIZE POPULATION NAMES TO v1, v2, v3, etc.
  standardized_result <- standardize_population_names(pop_filtered)
  pop_filtered <- standardized_result$data
  pop_mapping <- standardized_result$mapping
  
  # Clean samples data
  samples_clean <- samples_data %>%
    left_join(country_fix, by = "source_country", relationship = "many-to-many") %>%
    mutate(source_country = coalesce(standard_country, source_country)) %>%
    select(-standard_country) %>%
    left_join(state_mapping, by = c("source_state_or_region" = "full_name"), relationship = "many-to-many") %>%
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
  
  return(list(data = state_populations, mapping = pop_mapping))
}

#------------------------------PART 3-------------------------------------------
# --- 3. UPGRADED Map Integration for Bar Charts ---

# Getting brazil map data
brazil_states <- ne_states(country = "Brazil", returnclass = "sf")

# Standardize state codes in the map data
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

# UPGRADED function to create map with colorblind-friendly colors and consistent population naming
create_map_with_bar_charts_UPGRADED <- function(processed_result, dataset_name, k_value) {
  if (is.null(processed_result) || is.null(processed_result$data) || nrow(processed_result$data) == 0) {
    cat("DEBUG: processed_result is NULL or empty.\n")
    return(NULL)
  }
  
  state_pop_data <- processed_result$data
  pop_mapping <- processed_result$mapping
  
  cat("Creating UPGRADED map for", dataset_name, "with", nrow(state_pop_data), "data points\n")
  
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
  
  # Create COLORBLIND-FRIENDLY palette with consistent ordering
  all_populations <- sort(unique(plot_data$population))  # v1, v2, v3, etc. will sort correctly
  n_pops <- length(all_populations)
  
  colors <- get_colorblind_palette(n_pops)
  names(colors) <- all_populations
  
  cat("Using", n_pops, "colorblind-friendly colors for populations:", paste(all_populations, collapse = ", "), "\n")
  
  # Base map with improved styling
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
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(t = 10, b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  
    ) +
    labs(
      title = "Population Structure in Brazilian States",
      subtitle = paste(dataset_name, "- K =", k_value, "| Standardized populations (v1, v2, v3...) | >80% assignment confidence")
    )
  
  # Uniform bar charts with consistent size
  unique_states <- unique(plot_data$source_state_or_region)
  cat("Processing", length(unique_states), "states for bar charts\n")
  
  # Fixed dimensions for consistency
  fixed_width <- 1.6   
  fixed_height <- 2.2  
  
  for (state_code in unique_states) {
    state_df <- plot_data %>% filter(source_state_or_region == state_code)
    
    if (nrow(state_df) == 0) next
    
    # Ensure consistent ordering by population name (v1, v2, v3, etc.)
    state_df <- state_df %>% arrange(population)
    total_samples <- unique(state_df$total_samples)
    
    cat("State:", state_code, "Samples:", total_samples, "Populations:", paste(state_df$population, collapse = ", "), "\n")
    
    # Create bar chart with colorblind-friendly colors
    state_bar_chart <- ggplot(state_df, aes(x = 1, y = proportion, fill = population)) +
      geom_col(width = 1, color = "black", linewidth = 0.2) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
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
  
  # Create legend with improved styling
  dummy_data <- tibble(
    population = all_populations,
    y = 1
  )
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Population\n(Standardized)") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 11, family = "sans"),
      legend.title = element_text(size = 12, face = "bold", family = "sans"),
      legend.margin = margin(10, 10, 10, 10)
    )
  
  # Extract legend
  legend <- get_legend(legend_plot)
  
  # Combine map and legend with proper spacing
  final_plot <- plot_grid(
    p_map, 
    legend, 
    ncol = 2, 
    rel_widths = c(0.75, 0.25),  
    align = "h"
  ) +
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  # Save population mapping information
  mapping_info <- list(
    dataset_name = dataset_name,
    k_value = k_value,
    population_mapping = pop_mapping,
    total_populations = n_pops,
    colors_used = colors
  )
  
  return(list(plot = final_plot, mapping_info = mapping_info))
}

#--------------------------PART 4-----------------------------------------------
# --- 4. MAIN EXECUTION with UPGRADED Output Control ---

cat("=== BRAZIL POPULATION STRUCTURE MAP ANALYSIS - UPGRADED VERSION ===\n")

# Checking if samples_data is loaded correctly
if (!exists("samples_data") || !("sample" %in% names(samples_data))) {
  cat("ERROR: samples_data not found or missing 'sample' column.\n")
  stop("Please check data loading.")
}

cat("Samples data loaded with", nrow(samples_data), "rows\n")

# Processing each dataset and create UPGRADED plots
plots_created <- 0
failed_datasets <- c()
all_mapping_info <- list()

for (i in 1:length(pop_data_list)) {
  file_name <- names(pop_data_list)[i]
  pop_data <- pop_data_list[[i]]
  
  cat("\n", rep("=", 60), "\n")
  cat("PROCESSING UPGRADED DATASET", i, "of", length(pop_data_list), ":", file_name, "\n")
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
  
  # Process data with population standardization
  processed_result <- process_population_data(pop_data, dataset_type, samples_data)
  
  if (is.null(processed_result)) {
    cat("FAILED: No processed data for", dataset_type, "\n")
    failed_datasets <- c(failed_datasets, file_name)
    next
  }
  
  # Create UPGRADED map
  map_result <- create_map_with_bar_charts_UPGRADED(processed_result, dataset_type, k_value)
  
  if (!is.null(map_result) && !is.null(map_result$plot)) {
    # Display the plot
    print(map_result$plot)
    
    # Save plot to UPGRADE directory
    safe_name <- str_replace_all(tolower(dataset_type), "[^a-z0-9]+", "_") 
    filename <- paste0("UPGRADE/brazil_population_map_UPGRADED_", safe_name, "_k", k_value, ".png")
    ggsave(filename, map_result$plot, width = 16, height = 12, dpi = 300)
    cat("SUCCESS: UPGRADED plot saved as", filename, "\n")
    
    # Store mapping information
    all_mapping_info[[file_name]] <- map_result$mapping_info
    
    plots_created <- plots_created + 1
    
    # Add a pause to see each plot
    if (interactive()) {
      readline(prompt = "Press [Enter] to continue to next plot...")
    }
  } else {
    cat("FAILED: Could not create UPGRADED plot for", dataset_type, "\n")
    failed_datasets <- c(failed_datasets, file_name)
  }
}

#----------------------------PART 5---------------------------------------------
# SAVE MAPPING INFORMATION AND FINAL SUMMARY

# Save all mapping information to UPGRADE directory
saveRDS(all_mapping_info, "UPGRADE/population_mapping_info.rds")
cat("\nPopulation mapping information saved to UPGRADE/population_mapping_info.rds\n")

# Create a summary report
summary_report <- tibble(
  dataset = names(all_mapping_info),
  k_value = map_chr(all_mapping_info, ~ .x$k_value),
  total_populations = map_int(all_mapping_info, ~ .x$total_populations),
  dataset_name = map_chr(all_mapping_info, ~ .x$dataset_name)
)

write_csv(summary_report, "UPGRADE/population_analysis_summary.csv")
cat("Summary report saved to UPGRADE/population_analysis_summary.csv\n")

cat("\n=== FINAL UPGRADED SUMMARY ===\n")
cat("Expected datasets: 4\n")
cat("Total UPGRADED plots successfully created:", plots_created, "\n")
cat("Failed datasets:", length(failed_datasets), "\n")
if (length(failed_datasets) > 0) {
  cat("Failed dataset names:\n")
  for (failed in failed_datasets) {
    cat("  -", failed, "\n")
  }
}
cat("All files saved to UPGRADE/ directory\n")
cat("UPGRADED script completed with colorblind-friendly colors and standardized population names!\n")

save.image("UPGRADE/brazil_map_plotting_UPGRADED_data.Rda")
