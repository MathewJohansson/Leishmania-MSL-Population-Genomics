# CALCULATING NEI GENETIC DISTANCE — DIPLOID WHOLE-GENOME SNP CALLS
# Subset: Americas-only (AM, n=386; excludes Old World + L. donovani + outliers)
# Output: Figures/Phylogenetics Work/NJ Trees/Whole-Genome/Americas-only/
# LOCAL MODE: loads nj_tree from .nwk file only. No RData, no VCF.
# MSL0-only plots (D3.4, D3.5) commented out — require dosage_filtered from Viking RData.




## 1. CLEAN UP AND LOAD LIBRARIES ----------------------------------------------

rm(list = ls())
library(ape)
library(tidyverse)
library(ggtree)
library(gtools)
library(ggpubr)




## 2. PATHS --------------------------------------------------------------------

LOCAL_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

OUT_DIR  <- file.path(LOCAL_DIR, "Figures/Phylogenetics Work/NJ Trees/Whole-Genome/Americas-only")
DATA_OUT <- file.path(LOCAL_DIR, "Data/Phylogenetic_Data/NJ_Trees/Whole-Genome/Americas-only")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)




## 3. LOAD TREE ----------------------------------------------------------------

nj_tree <- read.tree(file.path(DATA_OUT, "diploid_whole_genome_Americas_only_Neis_tree.nwk"))




## 4. METADATA AND PALETTE -----------------------------------------------------

samples <- read_tsv(file.path(LOCAL_DIR, "Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv"))
names(samples)[1] <- "sra_run_accession"

pop_membership <- read_tsv(file.path(LOCAL_DIR, "Data/population_structure/population_membership_table_Americas_only_k13_AM_notation.tsv"))

# Named Paul Tol palette — AM notation (K=13)
paul_tol_palette <- setNames(
  c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
    "#AA3377", "#BBBBBB", "#EE7733", "#009988", "#994F00",
    "#332288", "#9955BB", "#FFAABB"),
  paste0("AM", 1:13)
)

# Region recode: Brazil → state; elsewhere → country name
region_recode <- c(
  "Le Vigan (30)"                                = "Gard",
  "Sommières (30)"                               = "Gard",
  "chodak namangan region in the pap district"   = "Namangan",
  "oltinkan namangan region in the pap district" = "Namangan",
  "Miyun (faubourg de Pékin)"                    = "Beijing",
  "Meshkin-Shahr (Est de l'Azerbaijan)"          = "East Azerbaijan",
  "Azerbaidjan (Est)"                            = "East Azerbaijan",
  "Hanovil"                                      = "Hanover",
  "Kassab"                                       = "Latakia",
  "Lisbonne"                                     = "Lisbon",
  "Marakech"                                     = "Marrakech",
  "Istambul"                                     = "Istanbul",
  "AYDIN"                                        = "Aydın"
)

# Coastal vs inland populations (for geographic plot)
coastal_pops <- c("AM2", "AM6", "AM9")
inland_pops  <- c("AM1", "AM10", "AM11")

# Build meta_combined from tree tip labels
meta_combined <- data.frame(sra_run_accession = nj_tree$tip.label) |>
  left_join(
    samples |> select(sra_run_accession, MSL_coverage_rounded,
                      source_country, source_state_or_region),
    by = "sra_run_accession"
  ) |>
  left_join(
    pop_membership |> select(sample, population, population_proportion),
    by = c("sra_run_accession" = "sample")
  ) |>
  mutate(
    display_region = if_else(source_country == "Brazil",
                             source_state_or_region, source_country),
    display_region = recode(display_region, !!!region_recode)
  ) |>
  filter(population_proportion >= 0.9 | is.na(population_proportion)) |>
  distinct(sra_run_accession, .keep_all = TRUE)

pops_present        <- mixedsort(unique(na.omit(meta_combined$population)))
pop_colours_present <- setNames(paul_tol_palette[pops_present], pops_present)

print(paste("Samples in meta_combined after filter:", nrow(meta_combined)))
print(paste("Populations present:", paste(pops_present, collapse = ", ")))




## 5. PLOT: MSL STATUS (panel A of D3.3) ---------------------------------------
# MSL plot uses a temporary tree with partial deletions (MSL1/2/3) dropped —
# only MSL0 and MSL4 samples shown here.

partial_deletion_ids <- samples |>
  filter(!is.na(MSL_coverage_rounded), !MSL_coverage_rounded %in% c(0, 4)) |>
  pull(sra_run_accession)

nj_tree_msl <- drop.tip(nj_tree,
                        intersect(partial_deletion_ids, nj_tree$tip.label))

tree_meta_msl <- data.frame(sra_run_accession = nj_tree_msl$tip.label) |>
  left_join(
    samples |> select(sra_run_accession, MSL_coverage_rounded),
    by = "sra_run_accession"
  ) |>
  mutate(MSL_coverage = factor(as.character(MSL_coverage_rounded), levels = c("0", "4")))

nj_plot.msl <- ggtree(nj_tree_msl, layout = "daylight") %<+% tree_meta_msl +
  geom_tippoint(aes(color = MSL_coverage), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#BB5566", "4" = "#4477AA"), name = "MSL Copies",
                     na.translate = FALSE) +
  theme(legend.position = "right")

ggsave(file.path(OUT_DIR, "diploid_whole_genome_Americas_only_NJ_tree_msl.jpeg"),
       plot = nj_plot.msl, width = 10, height = 10, dpi = 300)




## 6. PLOT: POPULATION (panel B of D3.3) ---------------------------------------

tree_meta_pop <- meta_combined |>
  select(sra_run_accession, population) |>
  filter(!is.na(population)) |>
  mutate(population = factor(population, levels = mixedsort(unique(na.omit(population)))))

unassigned_ids <- setdiff(nj_tree$tip.label, tree_meta_pop$sra_run_accession)
nj_tree_pop    <- drop.tip(nj_tree, unassigned_ids)

nj_plot.populations <- suppressWarnings(
  ggtree(nj_tree_pop, layout = "daylight") %<+% tree_meta_pop +
    geom_tippoint(aes(color = population), size = 2.5, alpha = 0.8) +
    scale_color_manual(values = pop_colours_present, name = "Population",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 8),
          legend.title = element_text(size = 10)) +
    guides(color = guide_legend(override.aes = list(size = 3)))
)

ggsave(file.path(OUT_DIR, "diploid_whole_genome_Americas_only_NJ_tree_population.jpeg"),
       plot = nj_plot.populations, width = 10, height = 10, dpi = 300)

# D3.3: Combined MSL + Population panel
nj_plot.msl.populations <- ggarrange(
  nj_plot.populations, nj_plot.msl,
  ncol = 2, nrow = 1, common.legend = FALSE, legend = "bottom",
  labels = c("A", "B"), font.label = list(size = 16, face = "bold")
)

ggsave(file.path(OUT_DIR, "diploid_whole_genome_Americas_only_NJ_tree_msl_populations.jpeg"),
       plot = nj_plot.msl.populations, width = 20, height = 10, dpi = 300)




## 7. PLOT: MSL0-ONLY BY POPULATION (D3.4) + GEOGRAPHIC (D3.5) -----------------
## COMMENTED OUT — requires dosage_filtered from full RData (too large for local)

# meta_msl0_only <- meta_combined |> filter(MSL_coverage_rounded == 0)
# ...




## 8. PLOT: BY COUNTRY (D3.1) --------------------------------------------------
# Uses ALL tree tips — 90% filter not applied

tree_meta_country <- data.frame(sra_run_accession = nj_tree$tip.label) |>
  left_join(
    samples |> select(sra_run_accession, source_country),
    by = "sra_run_accession"
  ) |>
  mutate(source_country = as.factor(source_country))

countries_present <- mixedsort(unique(na.omit(as.character(tree_meta_country$source_country))))
n_countries       <- length(countries_present)
country_palette   <- setNames(
  colorRampPalette(c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE",
                     "#AA3377","#BBBBBB","#EE7733","#009988","#994F00",
                     "#332288","#9955BB","#FFAABB","#BB5566","#44AA99"))(n_countries),
  countries_present
)

nj_plot_country <- suppressWarnings(
  ggtree(nj_tree, layout = "daylight") %<+% tree_meta_country +
    geom_tippoint(aes(color = source_country), size = 4, alpha = 0.85) +
    scale_color_manual(values = country_palette, name = "Country",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 18),
          legend.title = element_text(size = 20),
          legend.key.height = unit(1, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
)

ggsave(file.path(OUT_DIR, "diploid_whole_genome_Americas_only_NJ_tree_country.jpeg"),
       plot = nj_plot_country, width = 20, height = 20, dpi = 300)




## 9. PLOT: BY REGION (D3.2) ---------------------------------------------------
# Uses ALL tree tips — 90% filter not applied
# Region: Brazil → macro-region; all other countries → continent/region

# Non-Brazil countries shown as individual countries (Americas-only subset)
country_to_region <- c(
  "Honduras"  = "Honduras",
  "Paraguay"  = "Paraguay",
  "Uruguay"   = "Uruguay"
)

brazil_state_to_region <- c(
  "Pará"                = "Brazil North",
  "Bahia"               = "Brazil Northeast",
  "Ceara"               = "Brazil Northeast",
  "Maranhão"            = "Brazil Northeast",
  "Pernambuco"          = "Brazil Northeast",
  "Piauí"               = "Brazil Northeast",
  "Rio Grande do Norte" = "Brazil Northeast",
  "Sergipe"             = "Brazil Northeast",
  "Distrito Federal"    = "Brazil Central-West",
  "Mato Grosso"         = "Brazil Central-West",
  "Mato Grosso do Sul"  = "Brazil Central-West",
  "Espírito Santo"      = "Brazil Southeast",
  "Minas Gerais"        = "Brazil Southeast",
  "Rio de Janeiro"      = "Brazil Southeast",
  "São Paulo"           = "Brazil Southeast",
  "Rio Grande do Sul"   = "Brazil South",
  "Santa Catarina"      = "Brazil South"
)

tree_meta_region <- data.frame(sra_run_accession = nj_tree$tip.label) |>
  left_join(
    samples |> select(sra_run_accession, source_country, source_state_or_region),
    by = "sra_run_accession"
  ) |>
  mutate(
    display_region = if_else(
      source_country == "Brazil",
      brazil_state_to_region[source_state_or_region],
      country_to_region[source_country]
    ),
    display_region = as.factor(display_region)
  )

regions_present <- mixedsort(unique(na.omit(as.character(tree_meta_region$display_region))))
n_regions       <- length(regions_present)
region_palette  <- setNames(
  colorRampPalette(c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE",
                     "#AA3377","#BBBBBB","#EE7733","#009988","#994F00",
                     "#332288","#9955BB","#FFAABB","#BB5566","#44AA99"))(n_regions),
  regions_present
)

nj_plot_region <- suppressWarnings(
  ggtree(nj_tree, layout = "daylight") %<+% tree_meta_region +
    geom_tippoint(aes(color = display_region), size = 4, alpha = 0.85) +
    scale_color_manual(values = region_palette, name = "Region",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 20),
          legend.title = element_text(size = 22),
          legend.key.height = unit(1, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
)

ggsave(file.path(OUT_DIR, "diploid_whole_genome_Americas_only_NJ_tree_region.jpeg"),
       plot = nj_plot_region, width = 20, height = 20, dpi = 300)



