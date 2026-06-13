# WORLDWIDE POPULATION STRUCTURE MAP
# Modified from your original Brazil-focused code

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

# --- WORLDWIDE VERSION: Process population data globally ---

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
  
  return(location_populations)
}

# --- WORLDWIDE VERSION: Create map function ---

create_worldwide_map_with_bar_charts <- function(location_pop_data, dataset_name, k_value) {
  if (is.null(location_pop_data) || nrow(location_pop_data) == 0) {
    cat("DEBUG: location_pop_data is NULL or empty.\n")
    return(NULL)
  }
  
  cat("Creating worldwide map for", dataset_name, "with", nrow(location_pop_data), "data points\n")
  
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
  
  # Create color palette
  all_populations <- sort(unique(plot_data$population))
  n_pops <- length(all_populations)
  
  colors <- if (n_pops <= 12) {
    brewer.pal(max(3, n_pops), "Set3")[1:n_pops]
  } else {
    rainbow(n_pops)
  }
  names(colors) <- all_populations
  
  # Base world map
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
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(t = 10, b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  
    ) +
    labs(
      title = "Worldwide Population Structure",
      subtitle = paste(dataset_name, "- K =", k_value, "| Samples with >80% assignment confidence")
    )
  
  # Add bar charts for each location
  unique_locations <- unique(plot_data$location_id)
  cat("Processing", length(unique_locations), "locations for bar charts\n")
  
  # Adjust bar chart size for worldwide view
  fixed_width <- 2.5   # Slightly larger for worldwide view
  fixed_height <- 3.0  
  
  for (location in unique_locations) {
    location_df <- plot_data %>% filter(location_id == location)
    
    if (nrow(location_df) == 0) next
    
    location_df <- location_df %>% arrange(population)
    total_samples <- unique(location_df$total_samples)
    
    # Create bar chart
    location_bar_chart <- ggplot(location_df, aes(x = 1, y = proportion, fill = population)) +
      geom_col(width = 1, color = "black", linewidth = 0.2) +
      scale_fill_manual(values = colors, drop = FALSE) +
      scale_y_continuous(expand = c(0, 0)) +
      theme_void() +
      theme(
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
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
  
  # Create legend
  dummy_data <- tibble(
    population = all_populations,
    y = 1
  )
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Population") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.margin = margin(10, 10, 10, 10)
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
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  return(final_plot)
}

# --- EUROPE/ASIA FOCUSED VERSION ---

create_europe_asia_focused_map <- function(location_pop_data, dataset_name, k_value) {
  if (is.null(location_pop_data) || nrow(location_pop_data) == 0) {
    cat("DEBUG: location_pop_data is NULL or empty.\n")
    return(NULL)
  }
  
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
  
  # Create color palette
  all_populations <- sort(unique(plot_data$population))
  n_pops <- length(all_populations)
  
  colors <- if (n_pops <= 12) {
    brewer.pal(max(3, n_pops), "Set3")[1:n_pops]
  } else {
    rainbow(n_pops)
  }
  names(colors) <- all_populations
  
  # Base world map focused on Europe/Asia/Africa
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
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(t = 10, b = 5)),
      plot.subtitle = element_text(hjust = 0.5, size = 12, margin = margin(b = 10)),
      legend.position = "none",
      plot.margin = margin(t = 20, r = 10, b = 10, l = 10)  
    ) +
    labs(
      title = "Population Structure: Europe, Asia & Africa",
      subtitle = paste(dataset_name, "- K =", k_value, "| Samples with >80% assignment confidence")
    )
  
  # Add bar charts for each location
  unique_locations_plot <- unique(plot_data$location_id)
  cat("Adding", length(unique_locations_plot), "bar charts to map\n")
  
  # Bar chart size
  fixed_width <- 3.0   
  fixed_height <- 4.0  
  
  for (location in unique_locations_plot) {
    location_df <- plot_data %>% filter(location_id == location)
    
    if (nrow(location_df) == 0) next
    
    location_df <- location_df %>% arrange(population)
    total_samples <- unique(location_df$total_samples)
    
    # Create bar chart
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
  
  # Create legend
  dummy_data <- tibble(
    population = all_populations,
    y = 1
  )
  
  legend_plot <- ggplot(dummy_data, aes(x = 1, y = y, fill = population)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, name = "Population") +
    theme_void() +
    theme(
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),
      legend.text = element_text(size = 11),
      legend.title = element_text(size = 12, face = "bold"),
      legend.margin = margin(10, 10, 10, 10)
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
    theme(plot.margin = margin(t = 10, r = 10, b = 10, l = 10))
  
  return(final_plot)
}

# --- MAIN EXECUTION FOR EUROPE/ASIA FOCUSED PLOTS ---

cat("=== CREATING WORLDWIDE POPULATION STRUCTURE MAPS ===\n")

# Process each dataset for worldwide analysis
worldwide_plots_created <- 0
eurasia_plots_created <- 0

for (i in 1:length(pop_data_list)) {
  file_name <- names(pop_data_list)[i]
  pop_data <- pop_data_list[[i]]
  
  cat("\n", rep("=", 60), "\n")
  cat("PROCESSING WORLDWIDE DATASET", i, ":", file_name, "\n")
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
  
  # Process data worldwide
  location_pop_data <- process_population_data_worldwide(pop_data, dataset_type, samples_data)
  
  if (is.null(location_pop_data)) {
    cat("FAILED: No processed worldwide data for", dataset_type, "\n")
    next
  }
  
  # Create worldwide map (original)
  worldwide_map_plot <- create_worldwide_map_with_bar_charts(location_pop_data, dataset_type, k_value)
  
  if (!is.null(worldwide_map_plot)) {
    # Display the plot
    print(worldwide_map_plot)
    
    # Save plot 
    safe_name <- str_replace_all(tolower(dataset_type), "[^a-z0-9]+", "_") 
    filename <- paste0("worldwide_population_map_", safe_name, "_k", k_value, ".png")
    ggsave(filename, worldwide_map_plot, width = 16, height = 12, dpi = 300)
    cat("SUCCESS: Worldwide plot saved as", filename, "\n")
    worldwide_plots_created <- worldwide_plots_created + 1
  }
  
  # Create Europe/Asia/Africa focused map
  eurasia_map_plot <- create_europe_asia_focused_map(location_pop_data, dataset_type, k_value)
  
  if (!is.null(eurasia_map_plot)) {
    # Display the plot
    print(eurasia_map_plot)
    
    # Save plot 
    filename_eurasia <- paste0("eurasia_africa_population_map_", safe_name, "_k", k_value, ".png")
    ggsave(filename_eurasia, eurasia_map_plot, width = 16, height = 12, dpi = 300)
    cat("SUCCESS: Europe/Asia/Africa plot saved as", filename_eurasia, "\n")
    eurasia_plots_created <- eurasia_plots_created + 1
  }
  
  # Add a pause to see each plot
  if (interactive()) {
    readline(prompt = "Press [Enter] to continue to next dataset...")
  }
}

cat("\n=== WORLDWIDE PLOTS SUMMARY ===\n")
cat("Total worldwide plots successfully created:", worldwide_plots_created, "\n")
cat("Total Europe/Asia/Africa plots successfully created:", eurasia_plots_created, "\n")
cat("Analysis completed!\n")