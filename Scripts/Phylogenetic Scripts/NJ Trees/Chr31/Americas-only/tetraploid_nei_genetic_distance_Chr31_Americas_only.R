# CALCULATING NEI GENETIC DISTANCE — TETRAPLOID CHR31 SNP CALLS
# Subset: Americas-only (AM, n=332 MSL0/MSL4 strains)
# Input VCF: tetraploid.snps.filter.Americas_only_MSL0and4_samples.vcf.gz
# Output: Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only
# UPDATED: AM notation, named Paul Tol palette, region recode, separate country/region plots




## 1. CLEAN UP AND LOAD LIBRARIES ----------------------------------------------
rm(list = ls())
library(vcfR)
library(adegenet)
library(ape)
library(tidyverse)
library(StAMPP)
library(ggtree)
library(ggnewscale)
library(gtools)
library(ggpubr)




## 2. PATHS --------------------------------------------------------------------

LOCAL_DIR <- "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript"

# Define and create the output directory explicitly
OUT_DIR  <- file.path(LOCAL_DIR, "Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only")
DATA_OUT <- file.path(LOCAL_DIR, "Data/Phylogenetic_Data/NJ_Trees/Chr31/Americas-only")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DATA_OUT, recursive = TRUE, showWarnings = FALSE)




## 3. LOAD DATA ----------------------------------------------------------------

samples <- read_tsv("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
names(samples)[1] <- "sra_run_accession"

outliers_to_remove <- c("ERR205787", "SRR8608748")

# Named Paul Tol palette — AM notation
paul_tol_palette <- setNames(
  c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE","#AA3377","#BBBBBB",
    "#EE7733","#009988","#994F00","#332288","#9955BB","#FFAABB"),
  paste0("AM", 1:13)
)

vcf <- read.vcfR("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/tetraploid.snps.filter.Americas_only_MSL0and4_samples.vcf.gz")




## 4. PROCESS VCF FOR DOSAGE DATA ----------------------------------------------

chrom      <- getCHROM(vcf)
pos        <- as.numeric(getPOS(vcf))
keep_snps  <- !(chrom == "LinJ.31" & pos >= 1181282 & pos <= 1192406)
vcf_no_msl <- vcf[keep_snps, ]

gt_strings <- extract.gt(vcf_no_msl, element = "GT")

count_alt_alleles <- function(x) lengths(regmatches(x, gregexpr("1", x)))

dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings), ncol = ncol(gt_strings)
)
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)

dosage_matrix[is.na(gt_strings) | gt_strings == "./." | gt_strings == "./././."] <- NA
dosage_matrix <- t(dosage_matrix)

dosage_imputed <- apply(dosage_matrix, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE); x
})

dosage_imputed  <- dosage_imputed[!rownames(dosage_imputed) %in% outliers_to_remove, ]
print(paste("Removed outliers. New sample count:", nrow(dosage_imputed)))

column_vars     <- apply(dosage_imputed, 2, var, na.rm = TRUE)
keep_cols       <- !is.na(column_vars) & column_vars > 0
dosage_filtered <- dosage_imputed[, keep_cols]

print(paste("Total SNPs initially:", ncol(dosage_imputed)))
print(paste("Remaining SNPs for StAMPP:", ncol(dosage_filtered)))




## 5. PREPARE METADATA ---------------------------------------------------------

# V-notation → AM-notation lookup
notation_lookup <- c(
  V1="AM1", V2="AM2", V3="AM3", V4="AM4", V5="AM5", V6="AM6", V7="AM7",
  V8="AM8", V9="AM9", V10="AM10", V11="AM11", V12="AM12", V13="AM13"
)

# Region recode: Brazil → state; all other countries → country name
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

meta_combined <- data.frame(sra_run_accession = rownames(dosage_filtered)) |>
  left_join(
    samples |> select(sra_run_accession, population, MSL_coverage_rounded,
                      source_country, source_geographic_region,
                      source_state_or_region, host_species),
    by = "sra_run_accession"
  ) |>
  mutate(
    population     = notation_lookup[population],
    display_region = if_else(source_country == "Brazil",
                             source_state_or_region, source_country),
    display_region = recode(display_region, !!!region_recode)
  )

pops_present        <- mixedsort(unique(na.omit(meta_combined$population)))
pop_colours_present <- setNames(paul_tol_palette[pops_present], pops_present)

# Coastal vs inland populations (AM notation)
coastal_pops <- c("AM2", "AM6", "AM9")
inland_pops  <- c("AM1", "AM10", "AM11")




## 6. CALCULATE NEI'S GENETIC DISTANCE (StAMPP) --------------------------------

allele_freqs <- dosage_filtered / 4

stampp_df <- as.data.frame(allele_freqs)
stampp_df$Sample <- rownames(stampp_df)
stampp_df <- stampp_df |>
  left_join(meta_combined |> select(sra_run_accession, population),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = population) |>
  mutate(Pop = ifelse(is.na(Pop), "Unassigned", Pop),
         Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs)))

stampp_obj      <- stamppConvert(stampp_df, type = "r")
nei_dist_matrix <- stamppNeisD(stampp_obj, pop = FALSE)

stamppPhylip(nei_dist_matrix,
             file = "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/Chr31/Americas-only/tetraploid_chr31_Americas_only_Neis_distance.txt")

nj_tree <- nj(as.dist(nei_dist_matrix))
write.tree(nj_tree,
           file = "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/Chr31/Americas-only/tetraploid_chr31_Americas_only_Neis_tree.nwk")

print("StAMPP Nei's Genetic Distance calculated and exported successfully!")




## 6. PLOT: MSL STATUS (panel A of E2.3) ---------------------------------------

tree_meta <- meta_combined |>
  select(sra_run_accession, population, MSL_coverage_rounded) |>
  distinct() |>
  mutate(MSL_coverage = ifelse(is.na(MSL_coverage_rounded), "Unknown",
                               as.character(MSL_coverage_rounded)),
         MSL_coverage = as.factor(MSL_coverage))

nj_plot.msl <- ggtree(nj_tree, layout = "daylight") %<+% tree_meta +
  geom_tippoint(aes(color = MSL_coverage), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#BB5566", "4" = "#4477AA"), name = "MSL Copies") +
  theme(legend.position = "right",
        legend.text  = element_text(size = 16),
        legend.title = element_text(size = 18))

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only/tetraploid_chr31_NJ_tree_msl.jpeg",
  plot = nj_plot.msl, width = 10, height = 10, dpi = 300)




## 7. PLOT: POPULATION (panel B of E2.3) ---------------------------------------
# Drop unassigned samples (below 90% membership threshold) from population panel only.
# MSL, country, and region plots retain the full tree.

tree_meta_pop <- tree_meta |>
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
          legend.text  = element_text(size = 16),
          legend.title = element_text(size = 18)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only/tetraploid_chr31_NJ_tree_population.jpeg",
  plot = nj_plot.populations, width = 10, height = 10, dpi = 300)

nj_plot.msl.populations <- ggarrange(
  nj_plot.populations, nj_plot.msl,
  ncol = 2, nrow = 1, common.legend = FALSE, legend = "bottom",
  labels = c("A", "B"), font.label = list(size = 16, face = "bold")
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only/tetraploid_chr31_NJ_tree_msl_populations.jpeg",
  plot = nj_plot.msl.populations, width = 20, height = 10, dpi = 300)




## 8. PLOT: MSL0-ONLY BY POPULATION (E2.4) -------------------------------------

meta_msl0_only    <- meta_combined |> filter(MSL_coverage_rounded == 0)
dosage_msl0       <- dosage_filtered[rownames(dosage_filtered) %in% meta_msl0_only$sra_run_accession, ]
col_vars_msl0     <- apply(dosage_msl0, 2, var, na.rm = TRUE)
dosage_msl0_clean <- dosage_msl0[, !is.na(col_vars_msl0) & col_vars_msl0 > 0]

print(paste("MSL0 Samples:", nrow(dosage_msl0_clean)))
print(paste("Informative SNPs within MSL0 clade:", ncol(dosage_msl0_clean)))

allele_freqs_msl0 <- dosage_msl0_clean / 4
stampp_df_msl0    <- as.data.frame(allele_freqs_msl0)
stampp_df_msl0$Sample <- rownames(stampp_df_msl0)
stampp_df_msl0 <- stampp_df_msl0 |>
  left_join(meta_msl0_only |> select(sra_run_accession, population),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = population) |>
  mutate(Ploidy = 4, Format = "freq") |>
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs_msl0)))

stampp_obj_msl0 <- stamppConvert(stampp_df_msl0, type = "r")
nei_dist_msl0   <- stamppNeisD(stampp_obj_msl0, pop = FALSE)
nj_tree_msl0    <- nj(as.dist(nei_dist_msl0))

tree_meta_msl0 <- meta_msl0_only |>
  select(sra_run_accession, population) |>
  distinct() |>
  mutate(population = factor(population, levels = mixedsort(unique(population))))

nj_plot_msl0_pops <- suppressWarnings(
  ggtree(nj_tree_msl0, layout = "daylight") %<+% tree_meta_msl0 +
    geom_tippoint(aes(color = population), 
                  size = 3, 
                  alpha = 0.9) +
    scale_color_manual(values = pop_colours_present, name = "Population",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 12),
          legend.title = element_text(size = 14)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only/tetraploid_chr31_NJ_tree_MSL0_Only.jpeg",
  plot = nj_plot_msl0_pops, width = 12, height = 10, dpi = 300)




## 9. PLOT: MSL0 GEOGRAPHIC — COASTAL VS INLAND (E2.5) -------------------------

# Drop MSL0 samples not in either geographic zone (AM3, AM5, AM7, AM12)
unzoned_ids <- meta_msl0_only |>
  filter(!population %in% c(coastal_pops, inland_pops)) |>
  pull(sra_run_accession)

nj_tree_msl0_geo <- drop.tip(nj_tree_msl0,
                             intersect(unzoned_ids, nj_tree_msl0$tip.label))

tree_meta_msl0_geo <- data.frame(sra_run_accession = nj_tree_msl0_geo$tip.label) |>
  left_join(
    meta_msl0_only |> select(sra_run_accession, population),
    by = "sra_run_accession"
  ) |>
  distinct() |>
  mutate(
    geographic_zone = case_when(
      population %in% coastal_pops ~ "East Coast/South (Established Zone)",
      population %in% inland_pops  ~ "Inland/North (Expansion Zone)",
      TRUE ~ NA_character_
    ),
    geographic_zone = factor(geographic_zone,
                             levels = c("East Coast/South (Established Zone)", "Inland/North (Expansion Zone)"))
  )

my_cols_no_grey <- c(
  "East Coast/South (Established Zone)" = "#4477AA",
  "Inland/North (Expansion Zone)"       = "#EE7733"
)

nj_plot_msl0_geo <- suppressWarnings(
  ggtree(nj_tree_msl0_geo, layout = "daylight") %<+% tree_meta_msl0_geo +
    geom_tippoint(aes(color = geographic_zone, shape = geographic_zone),
                  size = 3.5, alpha = 0.9) +
    scale_color_manual(values = my_cols_no_grey, name = "", na.translate = FALSE) +
    scale_shape_manual(values = c(
      "East Coast/South (Established Zone)" = 16,
      "Inland/North (Expansion Zone)"       = 17), name = "", na.translate = FALSE) +
    theme(legend.position = "bottom",
          legend.text  = element_text(size = 21),
          legend.title = element_text(size = 23)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/Chr31/Americas-only/tetraploid_chr31_NJ_tree_MSL0_Geographic.jpeg",
  plot = nj_plot_msl0_geo, width = 16, height = 16, dpi = 300)




## 10. PLOT: BY COUNTRY (E2.1) -------------------------------------------------

tree_meta_country <- data.frame(sra_run_accession = nj_tree$tip.label) |>
  left_join(
    samples |> select(sra_run_accession, source_country),
    by = "sra_run_accession"
  ) |>
  mutate(source_country = as.factor(source_country))

# Consistent palette logic (Matches Whole-Genome aesthetic)
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
          legend.text  = element_text(size = 20),
          legend.title = element_text(size = 22),
          legend.key.height = unit(0.8, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
)

ggsave(
  file.path(OUT_DIR, "tetraploid_chr31_Americas_only_NJ_tree_country.jpeg"),
  plot = nj_plot_country, width = 20, height = 20, dpi = 300)




## 11. PLOT: BY REGION (E2.2) --------------------------------------------------

brazil_state_to_region <- c(
  "Pará" = "Brazil North", "Bahia" = "Brazil Northeast", "Ceara" = "Brazil Northeast",
  "Maranhão" = "Brazil Northeast", "Pernambuco" = "Brazil Northeast", 
  "Piauí" = "Brazil Northeast", "Rio Grande do Norte" = "Brazil Northeast", 
  "Sergipe" = "Brazil Northeast", "Distrito Federal" = "Brazil Central-West", 
  "Mato Grosso" = "Brazil Central-West", "Mato Grosso do Sul" = "Brazil Central-West",
  "Espírito Santo" = "Brazil Southeast", "Minas Gerais" = "Brazil Southeast", 
  "Rio de Janeiro" = "Brazil Southeast", "São Paulo" = "Brazil Southeast",
  "Rio Grande do Sul" = "Brazil South", "Santa Catarina" = "Brazil South"
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
      source_country 
    ),
    display_region = as.factor(display_region)
  )

# Hardcoded distinct palette for AM regions — avoids blue-heavy interpolation
region_palette <- c(
  "Brazil Central-West"      = "#CCBB44",
  "Brazil North"             = "#4477AA",
  "Brazil Northeast"         = "#EE6677",
  "Brazil South"             = "#AA3377",
  "Brazil Southeast"         = "#228833",
  "Honduras"                 = "#994F00",
  "Paraguay"                 = "#EE7733",
  "Uruguay"                  = "#66CCEE"
)

nj_plot_region <- suppressWarnings(
  ggtree(nj_tree, layout = "daylight") %<+% tree_meta_region +
    geom_tippoint(aes(color = display_region), size = 4, alpha = 0.85) +
    scale_color_manual(values = region_palette, name = "Region",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 20),
          legend.title = element_text(size = 22),
          legend.key.height = unit(0.8, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 5), ncol = 1))
)

ggsave(
  file.path(OUT_DIR, "tetraploid_chr31_Americas_only_NJ_tree_region.jpeg"),
  plot = nj_plot_region, width = 20, height = 20, dpi = 300)




## 12. SAVE WORKSPACE ----------------------------------------------------------

save.image("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/Chr31/Americas-only/tetraploid_chr31_Americas_only_Neis.RData")



