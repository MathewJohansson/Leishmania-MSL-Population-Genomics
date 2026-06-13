# NEI'S GENETIC DISTANCE — 100kb WINDOWED SNP CALLS
# Subset: L. infantum-only (LI, n=461; excludes 2 L. donovani + 2 outliers)
# Input VCF: tetraploid_LI_LinJ.31.MSL-centred-window-100kb.vcf.gz
# Ploidy: TETRAPLOID (4)
# NOTE: MSL locus sits within window but is excluded from distance calculations by StAMPP (no filter needed — MSL region has no informative SNPs after dosage filtering in MSL0 samples)
# Window: LinJ.31:1,136,844–1,236,844 (100kb centred on MSL midpoint 1,186,844)
# Output: Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/




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




## 2. LOAD DATA ----------------------------------------------------------------

samples <- read_tsv("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Metadata/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv")
names(samples)[1] <- "sra_run_accession"

ldonovani_to_exclude <- c("ERR3550127", "ERR205758")
outliers_to_remove   <- c("ERR205787", "SRR8608748", ldonovani_to_exclude)

# Named Paul Tol palette — LI notation
paul_tol_palette <- setNames(
  c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE","#AA3377","#BBBBBB",
    "#EE7733","#009988","#994F00","#332288","#9955BB","#FFAABB"),
  paste0("LI", 1:13)
)

vcf <- read.vcfR("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/100kb_windows/L.infantum-only/tetraploid_LI_LinJ.31.MSL-centred-window-100kb.vcf.gz")




## 3. PROCESS VCF FOR DOSAGE DATA ----------------------------------------------

gt_strings <- extract.gt(vcf, element = "GT")

count_alt_alleles <- function(x) lengths(regmatches(x, gregexpr("1", x)))

dosage_matrix <- matrix(
  count_alt_alleles(gt_strings),
  nrow = nrow(gt_strings), ncol = ncol(gt_strings)
)
rownames(dosage_matrix) <- rownames(gt_strings)
colnames(dosage_matrix) <- colnames(gt_strings)

# Tetraploid missing data format
dosage_matrix[is.na(gt_strings) | gt_strings == "././."] <- NA
dosage_matrix <- t(dosage_matrix)

dosage_imputed <- apply(dosage_matrix, 2, function(x) {
  x[is.na(x)] <- mean(x, na.rm = TRUE); x
})

dosage_imputed  <- dosage_imputed[!rownames(dosage_imputed) %in% outliers_to_remove, ]
print(paste("Removed L. donovani + outliers. New sample count:", nrow(dosage_imputed)))

column_vars     <- apply(dosage_imputed, 2, var, na.rm = TRUE)
keep_cols       <- !is.na(column_vars) & column_vars > 0
dosage_filtered <- dosage_imputed[, keep_cols]

print(paste("Total SNPs initially:", ncol(dosage_imputed)))
print(paste("Remaining SNPs for StAMPP:", ncol(dosage_filtered)))




## 4. PREPARE METADATA ---------------------------------------------------------

notation_lookup <- c(
  V1="LI1", V2="LI2", V3="LI3", V4="LI4", V5="LI5", V6="LI6", V7="LI7",
  V8="LI8", V9="LI9", V10="LI10", V11="LI11", V12="LI12", V13="LI13"
)

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

country_to_continent <- c(
  "Albania"                  = "Europe",
  "China"                    = "Asia",
  "Egypt"                    = "Africa",
  "France"                   = "Europe",
  "Georgia"                  = "Asia",
  "Greece"                   = "Europe",
  "Honduras"                 = "Central America",
  "Islamic Republic of Iran" = "Asia",
  "Israel"                   = "Middle East",
  "Italy"                    = "Europe",
  "Morocco"                  = "Africa",
  "Palestine"                = "Middle East",
  "Paraguay"                 = "South America (non-Brazil)",
  "Portugal"                 = "Europe",
  "Spain"                    = "Europe",
  "Syrian Arab Republic"     = "Middle East",
  "Tunisia"                  = "Africa",
  "Turkiye"                  = "Europe",
  "Uruguay"                  = "South America (non-Brazil)",
  "Uzbekistan"               = "Asia",
  "Yemen"                    = "Middle East",
  "Yugoslavia"               = "Europe"
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




## 5. CALCULATE NEI'S GENETIC DISTANCE (StAMPP) --------------------------------

allele_freqs <- dosage_filtered / 4   # TETRAPLOID

stampp_df <- as.data.frame(allele_freqs)
stampp_df$Sample <- rownames(stampp_df)
stampp_df <- stampp_df |>
  left_join(meta_combined |> select(sra_run_accession, population),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = population) |>
  mutate(Pop = ifelse(is.na(Pop), "Unassigned", Pop),
         Ploidy = 4, Format = "freq") |>   # TETRAPLOID
  select(Sample, Pop, Ploidy, Format, all_of(colnames(allele_freqs)))

stampp_obj      <- stamppConvert(stampp_df, type = "r")
nei_dist_matrix <- stamppNeisD(stampp_obj, pop = FALSE)

stamppPhylip(nei_dist_matrix,
             file = "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_Neis_distance.txt")

nj_tree <- nj(as.dist(nei_dist_matrix))
write.tree(nj_tree,
           file = "/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_Neis_tree.nwk")

print("StAMPP Nei's Genetic Distance calculated and exported successfully!")




## 6. PLOT: MSL STATUS (panel A of F1.3) ---------------------------------------

tree_meta <- meta_combined |>
  select(sra_run_accession, population, MSL_coverage_rounded) |>
  distinct() |>
  mutate(MSL_coverage = ifelse(is.na(MSL_coverage_rounded), "Unknown",
                               as.character(MSL_coverage_rounded)),
         MSL_coverage = as.factor(MSL_coverage))

# Prune tree to MSL0 and MSL4 only for MSL plot
msl_keep <- tree_meta |> filter(MSL_coverage %in% c("0", "4")) |> pull(sra_run_accession)
msl_drop  <- nj_tree$tip.label[!nj_tree$tip.label %in% msl_keep]
nj_tree_msl_only <- drop.tip(nj_tree, msl_drop)

nj_plot.msl <- ggtree(nj_tree_msl_only, layout = "daylight") %<+% tree_meta +
  geom_tippoint(aes(color = MSL_coverage), size = 2.5, alpha = 0.8) +
  scale_color_manual(values = c("0" = "#BB5566", "4" = "#4477AA"), name = "MSL Copies") +
  theme(legend.position = "right",
        legend.text  = element_text(size = 14),
        legend.title = element_text(size = 16))

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/tetraploid_100kb_Linfantum_only_NJ_tree_msl.jpeg",
  plot = nj_plot.msl, width = 10, height = 10, dpi = 300)




## 7. PLOT: POPULATION (panel B of F1.3) ---------------------------------------
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
          legend.text  = element_text(size = 14),
          legend.title = element_text(size = 16)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/tetraploid_100kb_Linfantum_only_NJ_tree_population.jpeg",
  plot = nj_plot.populations, width = 10, height = 10, dpi = 300)

nj_plot.msl.populations <- ggarrange(
  nj_plot.populations, nj_plot.msl,
  ncol = 2, nrow = 1, common.legend = FALSE, legend = "bottom",
  labels = c("A", "B"), font.label = list(size = 16, face = "bold")
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/tetraploid_100kb_Linfantum_only_NJ_tree_msl_populations.jpeg",
  plot = nj_plot.msl.populations, width = 20, height = 10, dpi = 300)




## 8. PLOT: MSL0-ONLY BY POPULATION (F1.4) -------------------------------------

meta_msl0_only    <- meta_combined |> filter(MSL_coverage_rounded == 0)
dosage_msl0       <- dosage_filtered[rownames(dosage_filtered) %in% meta_msl0_only$sra_run_accession, ]
col_vars_msl0     <- apply(dosage_msl0, 2, var, na.rm = TRUE)
dosage_msl0_clean <- dosage_msl0[, !is.na(col_vars_msl0) & col_vars_msl0 > 0]

print(paste("MSL0 Samples:", nrow(dosage_msl0_clean)))
print(paste("Informative SNPs within MSL0 clade:", ncol(dosage_msl0_clean)))

allele_freqs_msl0 <- dosage_msl0_clean / 4   # TETRAPLOID
stampp_df_msl0    <- as.data.frame(allele_freqs_msl0)
stampp_df_msl0$Sample <- rownames(stampp_df_msl0)
stampp_df_msl0 <- stampp_df_msl0 |>
  left_join(meta_msl0_only |> select(sra_run_accession, population),
            by = c("Sample" = "sra_run_accession")) |>
  rename(Pop = population) |>
  mutate(Ploidy = 4, Format = "freq") |>   # TETRAPLOID
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
    geom_tippoint(aes(color = population), size = 3, alpha = 0.9) +
    scale_color_manual(values = pop_colours_present, name = "Population",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 10),
          legend.title = element_text(size = 12)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/tetraploid_100kb_Linfantum_only_NJ_tree_MSL0_Only.jpeg",
  plot = nj_plot_msl0_pops, width = 12, height = 10, dpi = 300)




## 9. PLOT: BY COUNTRY (F1.1) --------------------------------------------------

tree_meta_country <- meta_combined |>
  select(sra_run_accession, population, source_country) |>
  distinct() |>
  mutate(population = factor(population, levels = mixedsort(unique(na.omit(population)))))

countries_present <- sort(unique(na.omit(tree_meta_country$source_country)))
country_colours   <- setNames(scales::hue_pal()(length(countries_present)), countries_present)

nj_plot_country <- suppressWarnings(
  ggtree(nj_tree, layout = "daylight") %<+% tree_meta_country +
    geom_tippoint(aes(color = source_country), size = 3.5, alpha = 0.85, shape = 16) +
    scale_color_manual(values = country_colours, name = "Country", na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 10),
          legend.title = element_text(size = 12)) +
    guides(color = guide_legend(override.aes = list(size = 4)))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/tetraploid_100kb_Linfantum_only_NJ_tree_country.jpeg",
  plot = nj_plot_country, width = 12, height = 12, dpi = 300)
print("Country-annotated NJ tree saved.")




## 10. PLOT: BY REGION (F1.2) --------------------------------------------------

tree_meta_region <- data.frame(sra_run_accession = nj_tree$tip.label) |>
  left_join(
    samples |> select(sra_run_accession, source_country, source_state_or_region),
    by = "sra_run_accession"
  ) |>
  mutate(
    display_region = if_else(
      source_country == "Brazil",
      brazil_state_to_region[source_state_or_region],
      country_to_continent[source_country]
    ),
    display_region = as.factor(display_region)
  )

regions_present <- mixedsort(unique(na.omit(as.character(tree_meta_region$display_region))))
n_regions       <- length(regions_present)

paul_tol_15 <- c("#4477AA","#EE6677","#228833","#CCBB44","#66CCEE",
                 "#AA3377","#BBBBBB","#EE7733","#009988","#994F00",
                 "#332288","#9955BB","#FFAABB","#BB5566","#44AA99")
region_palette <- setNames(paul_tol_15[seq_len(n_regions)], regions_present)

nj_plot_region <- suppressWarnings(
  ggtree(nj_tree, layout = "daylight") %<+% tree_meta_region +
    geom_tippoint(aes(color = display_region), size = 4, alpha = 0.85) +
    scale_color_manual(values = region_palette, name = "Region",
                       na.translate = FALSE) +
    theme(legend.position = "right",
          legend.text  = element_text(size = 16),
          legend.title = element_text(size = 18),
          legend.key.height = unit(0.8, "cm")) +
    guides(color = guide_legend(override.aes = list(size = 4), ncol = 1))
)

ggsave(
  "~/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Figures/Phylogenetics Work/NJ Trees/100kb Windows/L.infantum-only/tetraploid_100kb_Linfantum_only_NJ_tree_region.jpeg",
  plot = nj_plot_region, width = 20, height = 20, dpi = 300)




## 11. SAVE WORKSPACE ----------------------------------------------------------

save.image("/home/johansson/Documents/Bioinformatics/Leishmania/Leishmania Manuscript/Data/Phylogenetic_Data/NJ_Trees/100kb_windows/L.infantum-only/tetraploid_100kb_Linfantum_only_Neis.RData")



