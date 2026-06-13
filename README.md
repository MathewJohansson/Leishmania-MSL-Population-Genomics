# Population Genomics of Miltefosine Resistance in *Leishmania infantum*

## Overview

*Leishmania infantum* causes visceral leishmaniasis, affecting thousands globally. This project investigates the evolutionary dynamics of miltefosine resistance through population genomic analysis of the Miltefosine Sensitivity Locus (MSL) on chromosome 31.

**Central Research Question:** Is the spread of MSL deletions (conferring drug resistance) driven by positive selection or neutral demographic processes (e.g., drift)?

## Key Findings (Preliminary)

- Clear Old World/New World phylogenetic split at whole-genome and locus-specific scales.
- 13 distinct populations identified within *L. infantum* via ADMIXTURE (K=13).
- MSL deletions predominantly in Brazilian populations; absent or rare in Old World samples.
- Ploidy variation (diploid vs tetraploid) correlates with geography and MSL status.
- MSL locus shows divergent phylogenetic topology vs. genome-wide baseline, suggesting localized evolutionary processes.

---

## Analytical Approach

### 1. Hierarchical Population Structure (ADMIXTURE)

ADMIXTURE analysis identified cryptic population structure across three complementary datasets:

- **All samples + outgroup (K=12):** Species-level structure separating *L. infantum* from *L. donovani*.
- **L. infantum only (K=13):** Fine-scale cryptic populations within *L. infantum*.
- **Americas only (K=13):** High-resolution New World structure where MSL deletions predominate.

Cross-validation error minimisation determined optimal K for each dataset. All three datasets represented in the repository; focus figures below use K=13 *L. infantum* only as it directly addresses MSL variation within the species of interest.

**Status:** ✅ Complete  
**Figures:** `figures/admixture/`  
**Results:** `results/population_structure/`

---

### 2. Multi-Scale Phylogenetic Analysis

Neighbour-Joining phylogenies computed at four genomic scales to detect selection signatures on the MSL region:

- **Whole-genome:** Neutral evolutionary baseline across all chromosomes (diploid, multiple datasets).
- **Chromosome 31 (MSL locus):** Locus-specific evolution under potential selection (tetraploid).
- **100kb windows around MSL:** Regional selection signal (tetraploid, chr31 centred).
- **20kb windows around MSL:** Fine-scale locus-specific evolution (tetraploid, centred on MSL genes).

Trees stratified by ploidy level (diploid vs tetraploid) and annotated with geographic region, host species, MSL copy number, and population ancestry.

**Rationale:** If the MSL region is under selection, phylogenetic topology and branch lengths on chromosome 31 will deviate from genome-wide evolutionary patterns. Comparative phylogenetics at multiple scales isolates the geographic and genomic regions where selection is acting.

**Methods:** Nei genetic distances computed from filtered SNP matrices; Neighbour-Joining trees constructed with `ape::nj()` in R and visualised with `ggtree`. Distance matrices retained for NeighborNet analysis (SplitsTree). All phylogenetic figures annotated by region, country, population, host species, MSL copy number, and ploidy.

**Status:** ✅ Complete  
**Figures:** `figures/phylogenetics/nj_trees/`  
**Scripts:** `scripts/Phylogenetic Scripts/NJ Trees/`

---

### 3. Hardy-Weinberg Equilibrium Testing

Tests whether MSL allele frequencies deviate from neutral expectations within populations, distinguishing selection from demographic history.

**Methods:** Comprehensive HWE testing across all populations (≥90% membership threshold) and windowed HWE analysis on chromosome 31 (50kb windows, 5kb step). Two ploidy models: pseudodiploidised (for SNP-based tests) and raw tetraploid calls.

**Status:** ✅ Complete  
**Figures:** `figures/HWE/`  
**Results:** `tables/HWE/`  
**Scripts:** `scripts/HWE/`

---

### 4. Population Diversity & Differentiation

Comprehensive population genetic statistics quantifying within- and between-population genetic variation:

**Nucleotide Diversity (π):** Per-population nucleotide diversity on chr31 and genome-wide (50kb windows, 5kb step).

**Tajima's D:** Scaled measure of nucleotide diversity deviations; sensitive to demographic processes and selection. Per-population and genome-wide.

**Pairwise Fixation Index (Fst):** Measures population differentiation; identifies chr31 as divergent locus in some population pairs. Pairwise Fst computed across all ≥90% confident samples; sliding window Fst on chr31 (50kb windows, 5kb step).

**Nucleotide Divergence (Dxy):** Absolute divergence between MSL0 and MSL4 populations; quantifies evolutionary distance independent of within-population variation.

**Status:** ✅ Complete (individual analyses), ⏳ Unified framework pending final VCF validation  
**Figures:** `figures/population_diversity/`, `figures/unified_popgen/`  
**Results:** `tables/`, `data/Unified_Popgen/`  
**Scripts:** `scripts/Population Diversity and Differentiation/`

---

### 5. Geographic Distribution & Coverage Analysis

Spatial analysis of MSL copy number variation and sample distribution across global sampling locations.

**MSL Coverage:** Read depth-based inference of MSL copy number (0, 2, 3, or 4 copies). Coverage normalised to genome-wide median; stratified by geographic region, population, and host species.

**Geographic Maps:** Pie chart and bar chart overlays showing population composition and MSL copy number by country/region.

**Status:** ✅ Complete  
**Figures:** `figures/geographical_work/`, `figures/MSL Coverage Normalisation/`, `figures/MSL Verification Non-Americas/`  
**Scripts:** `scripts/Coverage/`, `scripts/Geographic Maps/`

---

## Clinical Relevance

MSL deletions confer resistance to **miltefosine**, the only oral treatment for leishmaniasis. Understanding whether resistant strains spread via:

- **Selection:** Faster than expected, requires urgent treatment policy intervention.
- **Drift:** Random spread driven by demographic history and founder effects.

...has direct implications for treatment policy, surveillance strategies, and prediction of resistance spread in new geographic regions.

---

## Methods Summary

**Population Structure (ADMIXTURE):**
- ADMIXTURE v1.3 with K=1-20 tested on pruned SNP matrices (LD-pruned, minor allele frequency ≥5%).
- Cross-validation error minimisation for optimal K selection.
- High-confidence assignments: ≥90% ancestry threshold for population membership.

**Phylogenetics:**
- Nei genetic distances computed from biallelic SNP genotypes (ape package, R).
- Neighbour-Joining trees constructed with `ape::nj()` and rooted where applicable.
- Visualisation with `ggtree` (R) including multiple annotation layers: MSL status, population, geographic region, host species, ploidy.
- Distance matrices and Newick-format tree files retained for reproducibility.

**Population Genetics:**
- π: windowed nucleotide diversity per population (50kb windows, 5kb step) via VCFtools.
- Tajima's D: standardised measure of nucleotide diversity, computed in 50kb windows with 5kb step.
- Fst: pairwise population differentiation (stampp package, R) and sliding-window Fst on chr31 via VCFtools.
- Dxy: nucleotide divergence between MSL0 and MSL4 populations (pegas package, R) and unified framework (custom R scripts).

**Hardy-Weinberg Equilibrium:**
- Comprehensive HWE: Chi-square test on diploidised genotypes for each population.
- Windowed HWE: Chi-square test on 50kb windows of chromosome 31 with 5kb step.
- Threshold for ≥90% population membership; separate ploidy models for SNP and tetraploid data.

**Data:**
- 465 *L. infantum* isolates, global sampling.
- Whole-genome sequencing (VCF format, biallelic SNPs).
- Geographic metadata, host species, MSL copy number (coverage-based).
- Ploidy stratification: diploid (whole-genome) and tetraploid (chr31) calls.

---

## Repository Structure

```
.
├── ADMIXTURE/                          # ADMIXTURE input/output files
│   └── Analysis Outputs/               # Population assignments, summaries, labels
├── Data/                               # Processed data and summary statistics
│   ├── Scripts/                        # Utility scripts (vcf2phylip)
│   ├── Unified_Popgen/                 # π, Tajima's D, Fst, Dxy results & verification
│   ├── sample_population_MSL_status.tsv
│   └── samples_outside_Americas_*.tsv
├── Figures/                            # High-quality publication-ready figures
│   ├── ADMIXTURE/                      # Population structure figures (3 datasets)
│   ├── Phylogenetics Work/             # IQ-TREE circular/rectangular trees
│   ├── Phylogenetics Work/NJ Trees/    # Neighbour-Joining trees at 4 scales
│   ├── HWE/                            # Hardy-Weinberg figures (comprehensive + windowed)
│   ├── Geographical Work/              # Maps, pie charts, coverage overlays
│   ├── Population Diversity and Differentiation/  # Fst, Tajima's D, π, Dxy
│   ├── MSL Coverage Normalisation/     # Coverage plots
│   ├── MSL Verification Non-Americas/  # Coverage verification for non-Americas strains
│   ├── Read Coverage/                  # Chromosome-wide coverage heatmaps
│   ├── Turkish Strain Coverage Verification/
│   └── Unified_Popgen/                 # Windowed population statistics figures
├── Metadata/                           # Sample metadata and project notes
│   ├── vcf_samples29_location_and_MSL_data_*.tsv
│   └── MSL Project Data.md
├── Scripts/                            # All analysis scripts (R, Perl, bash)
│   ├── ADMIXTURE/
│   ├── Coverage/
│   ├── Geographic Maps/
│   ├── HWE/
│   ├── Phylogenetic Scripts/
│   │   ├── IQ-TREE analyses (archived)
│   │   └── NJ Trees/                   # Current phylogenetic pipeline
│   ├── Population Diversity and Differentiation/
│   │   ├── Fst/
│   │   ├── PCA/
│   │   ├── Tajima's D/
│   │   ├── Pegas_Dxy/
│   │   └── Unified Population Statistics/
│   ├── diploidise_tetraploid_vcf_v2.pl
│   └── Script Inventory
├── Tables/                             # Summary tables and supplementary data
│   ├── HWE/                            # HWE results by population
│   ├── Chr31_TajD_Per_Population_Summary.csv
│   ├── population_mapping_table.tsv
│   └── *summary*.xlsx
├── Turkish Strain Coverage Verification/
├── README.md                           # This file
├── paul_tol_colorblind_palette_*.rds   # Colour palettes for figures
└── neis_txt_to_nex.py                  # Utility for distance matrix conversion
```

---

## Software & Dependencies

**R (v4.4.2):**
- Data manipulation: tidyverse, data.table
- Phylogenetics: ape, ggtree, treeio, phytools
- Population genetics: pegas, stampp, adegenet
- Visualisation: ggplot2, RColorBrewer, ggtext

**External Tools:**
- ADMIXTURE v1.3 (population structure)
- VCFtools (VCF processing, population genetics statistics)
- SplitsTree (NeighborNet, distance tree visualisation)

**Utilities:**
- vcf2phylip (VCF to phylip conversion)
- diploidise_tetraploid_vcf_v2.pl (tetraploid to pseudodiploid conversion)

---

## Data Sources

Genomic data: 465 *L. infantum* isolates, global sampling (Dr. Daniel Jeffares, University of York, unpublished).  
Geographic metadata and host species compiled from sampling records and literature.  
MSL copy number inferred from read depth mapping to reference genome.

---

## Project Status

**Completed:**
- ✅ ADMIXTURE population structure (three complementary datasets: species-level, species-only, Americas-only)
- ✅ Multi-scale phylogenetic analysis (Neighbour-Joining trees at whole-genome, chr31, 100kb, 20kb scales)
- ✅ Ploidy stratification (diploid and tetraploid analyses)
- ✅ Hardy-Weinberg equilibrium testing (comprehensive + windowed on chr31)
- ✅ Population diversity statistics (π, Tajima's D, Fst, Dxy)
- ✅ Geographic distribution and MSL coverage analysis
- ✅ MSL copy number verification in non-Americas strains

**Pending Final Validation:**
- ⏳ Unified population statistics framework (awaiting final VCF QC from supervisor)

**Broader Context:**
- 🔄 Manuscript preparation (collaborative, multi-author, in active development)
- 🔄 Integration of HWE and geographic collaborator results (Zeynep Sakaoglu)

---

## Author

**Mathew Johansson**  
MSc Bioinformatics, University of York  
Voluntary Research Associate, Jeffares Lab  
GitHub: [@MathewJohansson](https://github.com/MathewJohansson)

---

## Acknowledgements

- Dr. Daniel Jeffares (University of York) — supervision, genomic data generation, ongoing feedback
- Zeynep Sakaoglu — Hardy-Weinberg equilibrium and geographic distribution analyses collaboration

---

## Licence

**Code:** MIT Licence (upon manuscript publication)  
**Data:** Available upon reasonable request pending publication

---

## Key References

**ADMIXTURE:**  
Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655–1664.

**Phylogenetics & Distance-Based Methods:**  
Saitou, N. & Nei, M. (1987). The Neighbour-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution*, 4(4), 406–425.

**Population Genetics (VCFtools):**  
Danecek, P., Auton, A., Abecasis, G., et al. (2011). The variant call format and VCFtools. *Bioinformatics*, 27(15), 2156–2158.

**Visualisation (ggplot2, ggtree):**  
Wickham, H. (2016). *ggplot2: Elegant Graphics for Data Analysis* (2nd ed.). Springer.  
Yu, G., Smith, D.K., Zhu, H., Zhao, Y. & Guan, Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. *Methods in Ecology and Evolution*, 8(1), 28–36.

---

**Note:** This is an active research project demonstrating comprehensive bioinformatics workflow, reproducible research practices, and evolutionary/clinical genomics expertise. All scripts and figures are under development for multi-author manuscript publication.
