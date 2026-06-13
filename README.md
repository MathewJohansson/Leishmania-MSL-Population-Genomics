# Population Genomics of Miltefosine Resistance in *Leishmania infantum*

## Overview

*Leishmania infantum* causes visceral leishmaniasis (kala-azar), a severe and often fatal parasitic disease affecting \~300,000 people globally per year. For decades, miltefosine has been the only oral treatment available, however, miltefosine-resistant strains of *L. infantum* have emerged in South America, threatening treatment efficacy in endemic regions.

**The Problem:** Resistance to miltefosine is associated with deletions of the Miltefosine Sensitivity Locus (MSL), a 4-gene region on chromosome 31. If this resistance is spreading due to positive selection, intervention is urgent. If it spreads through neutral demographic processes (founder effects, genetic drift), the trajectory is more predictable but still concerning.

**This Project:** We investigate the evolutionary dynamics of MSL deletions across global *L. infantum* populations using population genomic analysis. By comparing phylogenetic and population genetic patterns at the MSL locus versus genome-wide baselines, we distinguish selection from drift with direct implications for predicting resistance spread and informing treatment policy.

**Central Research Question:** Is the spread of MSL deletions (conferring drug resistance) driven by positive selection or neutral demographic processes (e.g., drift)?

## Key Findings (Preliminary)

-   Clear Old World/New World phylogenetic split at whole-genome and locus-specific scales.
-   13 distinct populations identified within *L. infantum* via ADMIXTURE (K=13).
-   MSL deletions predominantly in Brazilian populations; absent or rare in Old World samples.
-   Ploidy variation (diploid vs tetraploid) correlates with geography and MSL status.
-   MSL locus shows divergent phylogenetic topology vs. genome-wide baseline, suggesting localised evolutionary processes.

------------------------------------------------------------------------

## Analytical Approach

### 1. Hierarchical Population Structure (ADMIXTURE)

ADMIXTURE analysis identified cryptic population structure across three complementary datasets:

-   **All samples + outgroup (K=12):** Species-level structure separating *L. infantum* from an *L. donovani* outgroup.
-   **L. infantum only (K=13):** Fine-scale cryptic populations within *L. infantum*.
-   **Americas only (K=13):** High-resolution New World structure where MSL deletions predominate.

Cross-validation error minimisation determined optimal K for each dataset. All three datasets are represented in the repository; focus figures below use K=13 *L. infantum* only as it directly addresses MSL variation within the species of interest.

**Status:** ✅ Complete\
**Figures:** `figures/admixture/`\
**Results:** `results/population_structure/`

------------------------------------------------------------------------

### 2. Multi-Scale Phylogenetic Analysis

Neighbour-Joining phylogenies were computed at four genomic scales to detect selection signatures on the MSL region at varying levels:

-   **Whole-genome:** Neutral evolutionary baseline across all chromosomes (diploid, multiple datasets).
-   **Chromosome 31 (MSL locus):** Locus-specific evolution under potential selection (tetraploid).
-   **100kb windows around MSL:** Regional selection signal (tetraploid, chr31 centred).
-   **20kb windows around MSL:** Fine-scale locus-specific evolution (tetraploid, centred on the MSL genes).

Trees were stratified by ploidy level (diploid vs tetraploid) and annotated with geographic region, host species, MSL copy number, and population ancestry.

**Rationale:** If the MSL region is under selection, phylogenetic topology and branch lengths on chromosome 31 will deviate from genome-wide evolutionary patterns. Comparative phylogenetics at multiple scales isolates the geographic and genomic regions where selection is acting.

**Methods:** Nei genetic distances computed from filtered SNP matrices; Neighbour-Joining trees constructed with `ape::nj()` in R and visualised with `ggtree`. Distance matrices retained for NeighborNet analysis (SplitsTree). All phylogenetic figures annotated by region, country, population, host species, MSL copy number, and ploidy.

**Status:** ✅ Complete\
**Figures:** `figures/phylogenetics/nj_trees/`\
**Scripts:** `scripts/Phylogenetic Scripts/NJ Trees/`

------------------------------------------------------------------------

### 3. Hardy-Weinberg Equilibrium Testing

Tests whether MSL allele frequencies deviate from neutral expectations within populations, distinguishing selection from demographic history.

**Methods:** Comprehensive HWE testing across all populations (≥90% membership threshold) and windowed HWE analysis on chromosome 31 (50kb windows, 5kb step). Two ploidy models: pseudodiploidised (for SNP-based tests) and raw tetraploid calls.

**Status:** ✅ Complete\
**Figures:** `figures/HWE/`\
**Results:** `tables/HWE/`\
**Scripts:** `scripts/HWE/`

------------------------------------------------------------------------

### 4. Population Diversity & Differentiation

Four comprehensive population genetic statistics quantifying within- and between-population genetic variation:

1.  **Nucleotide Diversity (π):** Per-population nucleotide diversity on chr31 and genome-wide (50kb windows, 5kb step).
2.  **Tajima's D:** Scaled measure of nucleotide diversity deviations; sensitive to demographic processes and selection. Per-population and genome-wide.
3.  **Pairwise Fixation Index (Fst):** Measures population differentiation; identifies chr31 as divergent locus in some population pairs. Pairwise Fst computed across all ≥90% confident samples; sliding window Fst on chr31 (50kb windows, 5kb step).
4.  **Nucleotide Divergence (Dxy):** Absolute divergence between MSL0 and MSL4 populations; quantifies evolutionary distance independent of within-population variation.

**Status:** ✅ Complete (individual analyses)

⏳ Unified framework pending final VCF validation. All four statistics will be combined into a single approach, where only two scripts will exist. This section will be updated once the unified work is complete.\
**Figures:** `figures/population_diversity/`, `figures/unified_popgen/`\
**Results:** `tables/`, `data/Unified_Popgen/`\
**Scripts:** `scripts/Population Diversity and Differentiation/`

------------------------------------------------------------------------

### 5. Geographic Distribution & Coverage Analysis

Spatial analysis of MSL copy number variation and sample distribution across global sampling locations.

**MSL Coverage:** Read depth-based inference of MSL copy number (0, 1, 2, 3, or 4 copies). Coverage normalised to genome-wide median; stratified by geographic region and population.

**Geographic Maps:** Pie chart and bar chart overlays showing population composition and MSL copy number by country/region.

![MSL Coverage by State](Figures/Geographical%20Work/brazil_coverage_map.jpg)

**Figure 1:** *Geographic distribution of MSL copy number variation across Brazilian states. Pie charts show the proportion of samples with 0 (red), 1 (orange), 2 (light yellow), 3 (light teal), or 4 (dark blue) MSL gene copies per location. Circle size indicates sample size (legend: 1, 2, 5, 10 samples). MSL deletions (0 copies, red) are concentrated in southeastern and northeastern Brazilian states, particularly São Paulo and Bahia, whilst northern and southern regions show higher frequency of intact MSL (blue, 4 copies). This geographic clustering suggests regional epidemiological patterns and potential selection or founder effects driving MSL resistance spread in Brazil.*

**Status:** ✅ Complete

**Figures:** `figures/geographical_work/`, `figures/MSL Coverage Normalisation/`, `figures/MSL Verification Non-Americas/`

**Scripts:** `scripts/Coverage/`, `scripts/Geographic Maps/`

------------------------------------------------------------------------

## Clinical Relevance

MSL deletions confer resistance to **miltefosine**, the only oral treatment for leishmaniasis. Understanding whether resistant strains spread via:

-   **Selection:** Faster than expected, requires urgent treatment policy intervention.
-   **Drift:** Random spread driven by demographic history and founder effects.

...has direct implications for treatment policy, surveillance strategies, and prediction of resistance spread in new geographic regions.

------------------------------------------------------------------------

## Methods Summary

**Population Structure (ADMIXTURE):**

-   ADMIXTURE v1.3 with K=1-20 tested on pruned SNP matrices.

-   Cross-validation error minimisation for optimal K selection.

-   High-confidence assignments: ≥90% ancestry threshold for population membership.

**Phylogenetics:**

-   Nei genetic distances computed from biallelic SNP genotypes (ape package, R).

-   Neighbour-Joining trees constructed with `ape::nj()` and rooted where applicable.

-   Visualisation with `ggtree` (R) including multiple annotation layers: MSL status, population, geographic region, host species, ploidy.

-   Distance matrices and Newick-format tree files retained for reproducibility.

**Population Genetics:**

-   π: windowed nucleotide diversity per population (50kb windows, 5kb step) via VCFtools.

-   Tajima's D: standardised measure of nucleotide diversity, computed in 50kb windows with 5kb step.

-   Fst: pairwise population differentiation (stampp package, R) and sliding-window Fst on chr31 via VCFtools.

-   Dxy: nucleotide divergence between MSL0 and MSL4 populations (pegas package, R) and unified framework (custom R scripts).

**Hardy-Weinberg Equilibrium:**

-   Comprehensive HWE: Chi-square test on diploidised genotypes for each population.

-   Windowed HWE: Chi-square test on 50kb windows of chromosome 31 with 5kb step.

-   Threshold for ≥90% population membership; separate ploidy models for SNP and tetraploid data.

**Data:**

-   463 *L. infantum* isolates, global sampling.

-   Whole-genome sequencing (VCF format, biallelic SNPs).

-   Geographic metadata, host species, MSL copy number (coverage-based).

-   Ploidy stratification: diploid (whole-genome) and tetraploid (chr31) calls.

------------------------------------------------------------------------

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
│   │   └── NJ Trees/
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
├── Turkish Strain Coverage Verification/
├── README.md                           # This file
├── paul_tol_colorblind_palette_*.rds   # Colour palettes for figures
└── neis_txt_to_nex.py                  # Utility for distance matrix conversion
```

------------------------------------------------------------------------

## Software & Computational Tools

**R (v4.4.2) & Packages:**

-   Data manipulation: tidyverse (Wickham et al. 2019), data.table (Dowle & Srinivasan 2021)

-   Phylogenetics: ape (Paradis & Schliep 2019), ggtree (Yu et al. 2017), treeio (Wang et al. 2020), phytools (Revell 2012)

-   Population genetics: pegas (Paradis 2010), stampp (Pembleton et al. 2013), adegenet (Jombart 2008)

-   Visualisation: ggplot2 (Wickham 2016), RColorBrewer (Neuwirth 2022), ggtext (Wilke & Wiernik 2022)

**External Tools:**

-   ADMIXTURE v1.3 (Alexander et al. 2009) — Bayesian population structure inference

-   VCFtools (Danecek et al. 2011) — VCF processing, filtering, population genetics statistics

-   SplitsTree (Huson & Bryant 2006) — NeighborNet and phylogenetic network visualisation

**Custom Utilities:**

-   `vcf2phylip` — VCF to PHYLIP/NEXUS format conversion

-   `diploidise_tetraploid_vcf_v2.pl` — Conversion of tetraploid to pseudodiploid genotypes for SNP-based analyses

-   A complete list of R package dependencies is available in `requirements.R`.

------------------------------------------------------------------------

## Data Sources

Genomic data: 465 *L. infantum* isolates, global sampling (Dr. Daniel Jeffares, University of York, unpublished).\
Geographic metadata and host species compiled from sampling records and literature.\
MSL copy number inferred from read depth mapping to reference genome.

------------------------------------------------------------------------

## Project Status

**Completed:**

-   ✅ ADMIXTURE population structure (three complementary datasets: species-level, species-only, Americas-only)

-   ✅ Multi-scale phylogenetic analysis (Neighbour-Joining trees at whole-genome, chr31, 100kb, 20kb scales)

-   ✅ Ploidy stratification (diploid and tetraploid analyses)

-   ✅ Hardy-Weinberg equilibrium testing (comprehensive + windowed on chr31)

-   ✅ Population diversity statistics (π, Tajima's D, Fst, Dxy)

-   ✅ Geographic distribution and MSL coverage analysis

-   ✅ MSL copy number verification in non-Americas strains

**Pending Final Validation:**

-   ⏳ Unified population statistics framework (awaiting final VCF QC from supervisor)

**Broader Context:**

-   🔄 Manuscript preparation (this is in active development, being prepared by myself, with assistance from Dr Daniel Jeffares

-   🔄 Integration of HWE and geographic collaborator results (Zeynep Sakaoglu)

------------------------------------------------------------------------

## Author

**Mathew Johansson**\
MSc Bioinformatics, University of York\
Voluntary Research Associate, Jeffares Lab\
GitHub: [\@MathewJohansson](https://github.com/MathewJohansson)

------------------------------------------------------------------------

## Acknowledgements

-   [Dr. Daniel Jeffares](https://www.york.ac.uk/biology/research/molecular-microbiology/daniel-jeffares/) (University of York) — supervision, ongoing feedback
-   [Zeynep Sakaoglu](https://tr.linkedin.com/in/zeynepsakaoglu) — Hardy-Weinberg equilibrium and geographic distribution analyses collaboration
-   [Dr. Eliza Cupolillo](https://scholar.google.com/citations?user=S5jebv6GW0IC&hl=pt-BR) (Instituto Oswaldo Cruz) — genomic data generation

------------------------------------------------------------------------

## Licence

**Code:** MIT Licence (upon manuscript publication)\
**Data:** Available upon reasonable request pending publication

------------------------------------------------------------------------

## References

Alexander, D.H., Novembre, J. & Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. Genome Research, 19(9), 1655–1664.

Danecek, P., Auton, A., Abecasis, G., et al. (2011). The variant call format and VCFtools. Bioinformatics, 27(15), 2156–2158.

Dowle, M. & Srinivasan, A. (2021). data.table: Extension of data.frame. R package version 1.14.2.

Huson, D.H. & Bryant, D. (2006). Application of phylogenetic networks in evolutionary studies. Molecular Biology and Evolution, 23(2), 254–267.

Jombart, T. (2008). adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics, 24(11), 1403–1405.

Neuwirth, E. (2022). RColorBrewer: ColorBrewer Palettes. R package version 1.1-3.

Paradis, E. (2010). pegas: an R package for population genetics with an integrated-modular environment. Bioinformatics, 26(3), 419–420.

Paradis, E. & Schliep, K. (2019). ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R. Bioinformatics, 35(3), 526–528.

Pembleton, L.W., Cogan, N.O.I. & Forster, J.W. (2013). StAMPP: an R package for calculation of genetic differentiation and structure of mixed-ploidy level populations. Molecular Ecology Resources, 13(5), 946–952.

Revell, L.J. (2012). phytools: an R package for phylogenetic comparative biology (and other things). Methods in Ecology and Evolution, 3(2), 217–223.

Saitou, N. & Nei, M. (1987). The Neighbour-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution, 4(4), 406–425.

Wang, L.G., Lam, T.T.Y., Xu, S., Dai, Z., Zhou, L., Feng, T., et al. (2020). treeio: an R package for phylogenetic tree input and output with richly annotated and associated data. Molecular Biology and Evolution, 37(2), 599–603.

Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis (2nd ed.). Springer.

Wickham, H., Averick, M., Bryan, J., Chang, W., McGowan, L.D., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T.L., Miller, E., Bache, S.M., Müller, K., Ooms, J., Robinson, D., Seidel, D., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. and Yutani, H. (2019). Welcome to the Tidyverse. *Journal of Open Source Software*, [online] 4(43), p.1686. doi:<https://doi.org/10.21105/joss.01686.>

Wilke, C.O. & Wiernik, B.M. (2022). ggtext: Improved Text Rendering Support for "ggplot2". R package version 0.1.1.

Yu, G., Smith, D.K., Zhu, H., Zhao, Y. & Guan, Y. (2017). ggtree: an R package for visualization and annotation of phylogenetic trees with their covariates and other associated data. Methods in Ecology and Evolution, 8(1), 28–36.

------------------------------------------------------------------------

**Note:** This is an active research project demonstrating comprehensive bioinformatics workflow, reproducible research practices, and evolutionary/clinical genomics expertise. All scripts and figures are under development for multi-author manuscript publication.
