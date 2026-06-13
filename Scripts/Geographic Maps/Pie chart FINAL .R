# MSL COVERAGE DISTRIBUTION MAPS
# L. infantum Genomics Study




# 1. SETUP AND CONFIGURATION ---------------------------------------------------

# 1.1. Load Required Libraries -------------------------------------------------

library(tidyverse)    
library(ggplot2)      
library(sf)         
library(rnaturalearth) 
library(rnaturalearthdata) 
library(ggforce)
library(cowplot)


# 1.2. Define Project Directories ----------------------------------------------

# Set base directory - ONLY CHANGE THIS LINE to run on different machines
BASE_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

METADATA_DIR <- file.path(BASE_DIR, "Metadata")
SCRIPTS_DIR  <- file.path(BASE_DIR, "Scripts")
FIGURES_DIR  <- file.path(BASE_DIR, "Figures/Geographical Work")

dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)


# 1.3. Load Metadata -----------------------------------------------------------

samples_data <- read_tsv(
  file.path(METADATA_DIR, "vcf_samples29_location_and_MSL_data_2025-06-19_corrected.tsv"),
  show_col_types = FALSE
) %>%
  rename(sample = sra_run_accession)



# Country name standardization 
country_fix <- tribble(
  ~source_country,   ~standard_country,
  "Brazil",          "Federative Republic of Brazil",
  "Paraguay",        "Republic of Paraguay",
  "Uruguay",         "Oriental Republic of Uruguay",
  "China",           "People's Republic of China",
  "Uzbekistan",      "Republic of Uzbekistan",
  "Albania",         "Republic of Albania",
  "France",          "French Republic",
  "Greece",          "Hellenic Republic",
  "Italy",           "Italian Republic",
  "Portugal",        "Portuguese Republic",
  "Spain",           "Kingdom of Spain",
  "Yugoslavia",      "Republic of Serbia",
  "Israel",          "State of Israel",
  "Turkiye",         "Republic of Turkey",
  "Yemen",           "Republic of Yemen",
  "Egypt",           "Arab Republic of Egypt",
  "Morocco",         "Kingdom of Morocco",
  "Tunisia",         "Tunisian Republic",
  "Honduras",        "Republic of Honduras",
  "Panama",          "Republic of Panama",
  "Palestine",       "West Bank",
  "Ethiopia",        "Federal Democratic Republic of Ethiopia",
  "Sudan",           "Republic of the Sudan",
  "Georgia",         "Georgia"
)

# Brazil state mapping
state_mapping <- tribble(
  ~full_name,             ~code,
  "Acre",                "AC",
  "Alagoas",             "AL",
  "Amapá",               "AP",
  "Amazonas",            "AM",
  "Bahia",               "BA",
  "Ceara",               "CE",
  "Distrito Federal",    "DF",
  "Espírito Santo",      "ES",
  "Goiás",               "GO",
  "Maranhão",            "MA",
  "Mato Grosso",         "MT",
  "Mato Grosso do Sul",  "MS",
  "Minas Gerais",        "MG",
  "Pará",                "PA",
  "Paraíba",             "PB",
  "Paraná",              "PR",
  "Pernambuco",          "PE",
  "Piauí",               "PI",
  "Rio de Janeiro",      "RJ",
  "Rio Grande do Norte", "RN",
  "Rio Grande do Sul",   "RS",
  "Rondônia",            "RO",
  "Roraima",             "RR",
  "Santa Catarina",      "SC",
  "São Paulo",           "SP",
  "Sergipe",             "SE",
  "Tocantins",           "TO"
)

# Apply fixes to data
cat("Processing samples data...\n")
samples_data_clean <- samples_data %>%
  left_join(country_fix, by = "source_country") %>%
  mutate(source_country = coalesce(standard_country, source_country)) %>%
  select(-standard_country) %>%
  left_join(state_mapping, by = c("source_state_or_region" = "full_name")) %>%
  mutate(source_state_or_region = if_else(!is.na(code), code, source_state_or_region)) %>%
  select(-code)

# Non-Americas partial deletions are artefactual — recode to 4 copies for plotting
non_americas_artefacts <- c(
  "SRR20748453", "SRR20748474", "SRR20748481", "SRR20748334",
  "ERR205736",   "ERR205760",   "ERR3550144",  "ERR271200",
  "ERR3550132",  "ERR205785",   "ERR304753",   "ERR205797",
  "ERR3550134",  "ERR3550142"
)

samples_data_clean <- samples_data_clean %>%
  mutate(
    # Override artefactual partial deletions in non-Americas samples
    coverage = if_else(sample %in% non_americas_artefacts, 4.0, coverage)
  )

# Creating coverage categories with ranges
samples_data_clean <- samples_data_clean %>%
  mutate(
    coverage_category = case_when(
      is.na(coverage) ~ "Missing",
      coverage >= 0 & coverage < 0.5 ~ "Cov0_Zero",      # Close to zero -> RED
      coverage >= 0.5 & coverage < 1.5 ~ "Cov1_One",     # Close to 1 -> ORANGE  
      coverage >= 1.5 & coverage < 2.5 ~ "Cov2_Two",     # Close to 2 -> YELLOW
      coverage >= 2.5 & coverage < 3.5 ~ "Cov3_Three",   # Close to 3 -> GREEN
      coverage >= 3.5 ~ "Cov4_Four",                      # Close to 4+ -> BLUE
      TRUE ~ "Other"
    )
  )

cat("Total samples processed:", nrow(samples_data_clean), "\n")

# Map data
cat("Loading map data...\n")
world_map <- ne_countries(scale = "medium", returnclass = "sf")
brazil_map <- ne_states(country = "Brazil", returnclass = "sf")

# Generate coordinate lookup tables from sf map objects
# Disable S2 geometry to avoid validity errors with some polygons
sf_use_s2(FALSE)

brazil_coords <- brazil_map %>%
  st_centroid() %>%
  mutate(
    centroid_longitude = st_coordinates(.)[, 1],
    centroid_latitude  = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(subdivision = postal, centroid_longitude, centroid_latitude)

world_coords <- world_map %>%
  st_centroid() %>%
  mutate(
    LONG = st_coordinates(.)[, 1],
    LAT  = st_coordinates(.)[, 2]
  ) %>%
  st_drop_geometry() %>%
  select(FULL_NAME = sovereignt, LAT, LONG)

# Build a bridge from formal country names (used in samples_data_clean) back to sovereignt names
# so world_coords join works correctly
# Deduplicate world_coords to mainland only — exclude overseas territories
# by keeping the single coordinate per country closest to its expected mainland location
world_coords_mainland <- world_coords %>%
  group_by(FULL_NAME) %>%
  slice(which.min(
    # For each country keep the entry with longitude closest to 0 and latitude closest to 30
    # This heuristic selects mainland Europe/Asia over distant overseas territories
    abs(LAT - 30) + abs(LONG - 20)
  )) %>%
  ungroup()

# Override specific known-problematic countries with correct mainland coordinates
mainland_overrides <- tribble(
  ~FULL_NAME,                        ~LAT,   ~LONG,
  "French Republic",                  46.0,    2.0,
  "France",                           46.0,    2.0,
  "People's Republic of China",       35.0,  105.0,
  "China",                            35.0,  105.0,
  "State of Israel",                  31.5,   34.9,
  "Israel",                           31.5,   34.9,
  "Republic of Serbia",               44.2,   20.8,
  "Georgia",                          42.2,   43.5
)

world_coords_mainland <- world_coords_mainland %>%
  rows_update(mainland_overrides, by = "FULL_NAME", unmatched = "ignore")

country_to_sovereignt <- country_fix %>%
  mutate(sovereignt = source_country) %>%
  select(FULL_NAME = standard_country, sovereignt) %>%
  left_join(world_coords_mainland %>% select(FULL_NAME, LAT, LONG),
            by = c("sovereignt" = "FULL_NAME")) %>%
  select(FULL_NAME, LAT, LONG) %>%
  bind_rows(world_coords_mainland) %>%
  group_by(FULL_NAME) %>%
  slice(1) %>%
  ungroup()

# Define colors 
coverage_colors <- c(
  "Cov0_Zero" = "#A50026",     
  "Cov1_One" = "#F99858",      
  "Cov2_Two" = "#EAECCC",     
  "Cov3_Three" = "#98CAE1",    
  "Cov4_Four" = "#364B9A"        
)




# 2. DIAGNOSTIC INFORMATION ----------------------------------------------------

cat("COVERAGE DISTRIBUTION CHECK")
coverage_table <- table(samples_data_clean$coverage_category, useNA = "ifany")
cat("Coverage distribution by ranges:\n")
print(coverage_table)
cat("Percentage distribution:\n")
print(round(prop.table(coverage_table) * 100, 1))

cat("\n=== DIAGNOSTIC INFO ===\n")
cat("Unique countries in samples:", length(unique(samples_data_clean$source_country)), "\n")
cat("Brazil samples:", sum(samples_data_clean$source_country == "Federative Republic of Brazil", na.rm = TRUE), "\n")
cat("Brazil coordinate columns:", paste(colnames(brazil_coords), collapse = ", "), "\n")
cat("World coordinate columns:", paste(colnames(country_to_sovereignt), collapse = ", "), "\n")

# Show some sample data
cat("\nSample countries in data:\n")
print(head(unique(samples_data_clean$source_country), 10))

cat("\nSample Brazil states in data:\n")
brazil_states <- unique(samples_data_clean$source_state_or_region[
  samples_data_clean$source_country == "Federative Republic of Brazil"
])
print(head(brazil_states, 10))


# 2.1. BRAZIL MAP
cat("PROCESSING BRAZIL DATA")
# Filter to AM subset samples only, then to Brazil
am_samples <- read_tsv(
  file.path("/home/johansson/Desktop/Leishmania Manuscript/Data/population_structure/population_membership_table_Americas_only_k13_AM_notation.tsv"),
  show_col_types = FALSE
) %>%
  select(sample = 1)

brazil_samples <- samples_data_clean %>%
  inner_join(am_samples, by = "sample") %>%
  filter(source_country == "Federative Republic of Brazil")

cat("Brazil samples found:", nrow(brazil_samples), "\n")

if(nrow(brazil_samples) > 0) {
  cat("Available states in brazil_coords:\n")
  print(head(unique(brazil_coords$subdivision), 10))
  
  brazil_data <- brazil_samples %>%
    left_join(brazil_coords, by = c("source_state_or_region" = "subdivision"))
  
  cat("After joining coordinates - rows with coordinates:", 
      sum(!is.na(brazil_data$centroid_longitude) & !is.na(brazil_data$centroid_latitude)), "\n")
  
  brazil_data <- brazil_data %>%
    filter(!is.na(centroid_longitude) & !is.na(centroid_latitude))
  
  if(nrow(brazil_data) > 0) {
    cat("Creating Brazil pie chart data...\n")
    
    canonical_order <- c("Cov0_Zero", "Cov1_One", "Cov2_Two", "Cov3_Three", "Cov4_Four")
    
    # Keep data in long format for ggforce - no pivot_wider needed
    # Radius formula: pmax(0.3, pmin(total * 0.2, 4))
    brazil_pie <- brazil_data %>%
      count(source_state_or_region, centroid_longitude, centroid_latitude, coverage_category) %>%
      mutate(coverage_category = factor(coverage_category, levels = canonical_order)) %>%
      group_by(source_state_or_region, centroid_longitude, centroid_latitude) %>%
      arrange(source_state_or_region, coverage_category) %>%
      mutate(
        total = sum(n),
        prop  = n / total,
        pie_r = pmax(0.3, pmin(total * 0.2, 2)),
        end_angle   = cumsum(prop) * 2 * pi,
        start_angle = lag(end_angle, default = 0)
      ) %>%
      ungroup() %>%
      mutate(coverage_category = as.character(coverage_category))
    
    cat("Brazil pie data created with", length(unique(brazil_pie$source_state_or_region)), "locations\n")
    
    cat("Creating Brazil plot...\n")
    
    brazil_pie$coverage_category <- factor(
      brazil_pie$coverage_category,
      levels = canonical_order
    )
    
    p_brazil <- ggplot() +
      geom_sf(data = brazil_map, fill = "lightgray", color = "white", size = 0.3) +
      geom_arc_bar(
        data = brazil_pie,
        aes(
          x0 = centroid_longitude, y0 = centroid_latitude,
          r0 = 0, r = pie_r,
          start = start_angle, end = end_angle,
          fill = coverage_category
        ),
        color = "black", linewidth = 0.3, alpha = 0.8
      ) +
      scale_fill_manual(
        values = coverage_colors,
        breaks = canonical_order,
        limits = canonical_order,
        labels = c(
          "Cov0_Zero"  = "0",
          "Cov1_One"   = "1",
          "Cov2_Two"   = "2",
          "Cov3_Three" = "3",
          "Cov4_Four"  = "4"
        ),
        name = "Coverage:"
      ) +
      guides(fill = guide_legend(ncol = 5)) + 
      coord_sf() +
      theme_void() +
      theme(
        plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
        legend.position = "bottom",
        legend.text  = element_text(size = 10)
      ) +
      labs(title = "Brazil: Coverage Distribution by State")
    
    cat("Displaying Brazil plot...\n")
    print(p_brazil)
    
    # Size legend inset for Brazil
    # Radius formula: pmax(0.3, pmin(total * 0.2, 4))
    # Representative n values spanning full data range (n=1 to n=130)
    legend_ns_br <- c(10, 5, 2, 1)
    legend_rs_br <- pmax(0.3, pmin(legend_ns_br * 0.2, 2))
    
    gap_br          <- 0.3
    y_centres_br    <- numeric(length(legend_rs_br))
    y_centres_br[1] <- legend_rs_br[1]
    for (i in 2:length(legend_rs_br)) {
      y_centres_br[i] <- y_centres_br[i - 1] - legend_rs_br[i - 1] - gap_br - legend_rs_br[i]
    }
    
    legend_df_br <- data.frame(
      x     = rep(0, length(legend_ns_br)),
      y     = y_centres_br,
      r     = legend_rs_br,
      label = as.character(legend_ns_br)
    )
    
    p_size_legend_br <- ggplot(legend_df_br) +
      geom_circle(aes(x0 = x, y0 = y, r = r), fill = "grey70", colour = "black",
                  linewidth = 0.3) +
      geom_text(aes(x = max(legend_rs_br) + 0.2, y = y, label = label),
                size = 2.8, hjust = 0, vjust = 0.5) +
      annotate("text", x = 0, y = y_centres_br[1] + legend_rs_br[1] + 0.5,
               label = "Sample sizes", size = 2.8, hjust = 0.5) +
      coord_fixed(clip = "off",
                  xlim = c(-max(legend_rs_br), max(legend_rs_br) + 2),
                  ylim = c(min(y_centres_br) - max(legend_rs_br),
                           y_centres_br[1] + legend_rs_br[1] + 1)) +
      theme_void()
    
    p_brazil_final <- ggdraw(p_brazil) +
      draw_plot(p_size_legend_br,
                x = 0.01, y = 0.01,
                width = 0.10, height = 0.40)
    
    ggsave(file.path(FIGURES_DIR, "brazil_coverage_map.jpg"),
           p_brazil_final, width = 12, height = 8, dpi = 300)
    cat("Brazil map saved\n")
  } else {
    cat("No Brazil data with coordinates found\n")
  }
} else {
  cat("No Brazil samples found in dataset\n")
}



# 2.2. WORLD MAP  

cat("PROCESSING WORLD DATA")
world_data <- samples_data_clean %>%
  left_join(country_to_sovereignt, by = c("source_country" = "FULL_NAME"))

cat("After joining world coordinates - rows with coordinates:", 
    sum(!is.na(world_data$LAT) & !is.na(world_data$LONG)), "\n")

cat("Sample country matches:\n")
sample_matches <- samples_data_clean %>%
  left_join(country_to_sovereignt, by = c("source_country" = "FULL_NAME")) %>%
  filter(!is.na(LAT)) %>%
  select(source_country, LAT, LONG) %>%
  distinct() %>%
  head(5)
print(sample_matches)

world_data <- world_data %>%
  filter(!is.na(LAT) & !is.na(LONG))

if(nrow(world_data) > 0) {
  cat("Creating world pie chart data...\n")
  
  canonical_order <- c("Cov0_Zero", "Cov1_One", "Cov2_Two", "Cov3_Three", "Cov4_Four")
  
  # Keep data in long format for ggforce - no pivot_wider needed
  # Radius formula: pmax(3, pmin(total^0.4 * 2.0, 14))
  world_pie <- world_data %>%
    count(source_country, LONG, LAT, coverage_category) %>%
    mutate(coverage_category = factor(coverage_category, levels = canonical_order)) %>%
    group_by(source_country, LONG, LAT) %>%
    arrange(source_country, coverage_category) %>%
    mutate(
      total = sum(n),
      prop  = n / total,
      pie_r = pmax(3, pmin(total^0.4 * 2.0, 14)),
      end_angle   = cumsum(prop) * 2 * pi,
      start_angle = lag(end_angle, default = 0)
    ) %>%
    ungroup() %>%
    rename(longitude = LONG, latitude = LAT) %>%
    mutate(coverage_category = as.character(coverage_category))
  
  cat("World pie data created with", length(unique(world_pie$source_country)), "locations\n")
  
  cat("Creating world plot...\n")
  
  world_pie$coverage_category <- factor(
    world_pie$coverage_category,
    levels = canonical_order
  )
  
  p_world <- ggplot() +
    geom_sf(data = world_map, fill = "lightgray", color = "white", size = 0.2) +
    geom_arc_bar(
      data = world_pie,
      aes(
        x0 = longitude, y0 = latitude,
        r0 = 0, r = pie_r,
        start = start_angle, end = end_angle,
        fill = coverage_category
      ),
      color = "black", linewidth = 0.3, alpha = 0.8
    ) +
    scale_fill_manual(
      values = coverage_colors,
      breaks = canonical_order,
      limits = canonical_order,
      labels = c(
        "Cov0_Zero"  = "0",
        "Cov1_One"   = "1",
        "Cov2_Two"   = "2",
        "Cov3_Three" = "3",
        "Cov4_Four"  = "4"
      ),
      name = "Coverage:"
    ) +
    guides(fill = guide_legend(ncol = 5)) +
    coord_sf(xlim = c(-180, 180), ylim = c(-60, 80)) +
    theme_void() +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 14, face = "bold"),
      legend.position = "bottom",
      legend.text  = element_text(size = 10)
    ) +
    labs(title = "Global Coverage Distribution by Country")
  
  # Size legend inset for world map
  # Radius formula: pmax(3, pmin(total^0.4 * 2.0, 14))
  # Representative n values spanning full data range (n=1 to n=378)
  legend_ns <- c(378, 50, 10, 1)   # largest first (top to bottom)
  legend_rs <- pmax(3, pmin(legend_ns^0.4 * 2.0, 14))
  
  gap       <- 2
  y_centres <- numeric(length(legend_rs))
  y_centres[1] <- legend_rs[1]
  for (i in 2:length(legend_rs)) {
    y_centres[i] <- y_centres[i - 1] - legend_rs[i - 1] - gap - legend_rs[i]
  }
  
  legend_df <- data.frame(
    x     = rep(0, length(legend_ns)),
    y     = y_centres,
    r     = legend_rs,
    label = as.character(legend_ns)
  )
  
  p_size_legend <- ggplot(legend_df) +
    geom_circle(aes(x0 = x, y0 = y, r = r), fill = "grey70", colour = "black",
                linewidth = 0.3) +
    geom_text(aes(x = max(legend_rs) + 1.5, y = y, label = label),
              size = 2.8, hjust = 0, vjust = 0.5) +
    annotate("text", x = 0, y = y_centres[1] + legend_rs[1] + 3,
             label = "Sample sizes", size = 2.8, hjust = 0.5) +
    coord_fixed(clip = "off",
                xlim = c(-max(legend_rs), max(legend_rs) + 8),
                ylim = c(min(y_centres) - max(legend_rs),
                         y_centres[1] + legend_rs[1] + 6)) +
    theme_void()
  
  p_world_final <- ggdraw(p_world) +
    draw_plot(p_size_legend,
              x = 0.00, y = 0.00,
              width = 0.09, height = 0.45)
  
  print(p_world_final)
  
  ggsave(file.path(FIGURES_DIR, "world_coverage_map.jpg"),
         p_world_final, width = 14, height = 8, dpi = 300)
  cat("World map saved\n")
} else {
  cat("No world data with coordinates found\n")
}




# 3. FINAL SUMMARY -------------------------------------------------------------

cat("FINAL SUMMARY")
coverage_table_final <- table(samples_data_clean$coverage_category)
cat("Coverage distribution by ranges:\n")
print(coverage_table_final)
cat("Percentage distribution:\n")
print(round(prop.table(coverage_table_final) * 100, 1))

cat("\nActual coverage value examples by category:\n")
coverage_examples <- samples_data_clean %>%
  group_by(coverage_category) %>%
  summarise(
    min_coverage = min(coverage, na.rm = TRUE),
    max_coverage = max(coverage, na.rm = TRUE),
    mean_coverage = mean(coverage, na.rm = TRUE),
    count = n()
  ) %>%
  arrange(coverage_category)
print(coverage_examples)

cat("\nScript completed!\n")