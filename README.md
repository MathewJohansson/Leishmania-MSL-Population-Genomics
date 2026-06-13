# Population Genomics of Miltefosine Resistance in *Leishmania infantum*

## Overview

*Leishmania infantum* causes visceral leishmaniasis, affecting thousands globally. This project investigates the evolutionary dynamics of miltefosine resistance through population genomic analysis of the Miltefosine Sensitivity Locus (MSL) on chromosome 31.

**Central Research Question:** Is the spread of MSL deletions (conferring drug resistance) driven by positive selection or neutral demographic processes (e.g., drift)?

## Key Findings (Preliminary)

- Clear Old World/New World phylogenetic split.
- 13 distinct populations identified within *L. infantum* via ADMIXTURE.
- MSL deletions predominantly in Brazilian populations.
- Ploidy variation (diploid vs tetraploid) correlates with geography.

---

## Analytical Approach

### 1. Hierarchical Population Structure (ADMIXTURE)

ADMIXTURE analysis identified cryptic population structure across three complementary datasets:

- **All samples + outgroup (K=12):** Species-level structure separating *L. infantum* from *L. donovani*.
- **L. infantum only (K=13):** Fine-scale cryptic populations within *L. infantum*.
- **Americas only (K=13):** High-resolution New World structure where MSL deletions predominate.

Cross-validation error minimisation determined optimal K for each dataset.

*The figures below focus on the K=13 L. infantum only analysis, as this dataset directly addresses MSL variation within the species of interest. Complete results for all three datasets are available in the `figures/admixture/` directory.*

#### Cross-Validation Error Justification (K=13)

![Cross-Validation Error K=13](readme_images/linfantum_only_cv_error_plot_k1_to_k20_for_k13_only.jpeg)

*Figure 1: Cross-validation error minimisation confirms K=13 as optimal for L. infantum population structure. The lowest CV error at K=13 indicates this number of populations best explains the genetic structure within the species where MSL variation occurs.*

#### Population Assignments (K=13, L. infantum Only)

![ADMIXTURE K=13 Population Assignments](readme_images/linfantum_only_population_assignments_k13.jpeg)

*Figure 2: ADMIXTURE-inferred population assignments (K=13) showing 13 distinct populations within L. infantum. Brazilian populations dominate the dataset, with clear geographic and genetic structure. This fine-scale resolution reveals cryptic subdivisions relevant to MSL resistance patterns.*

#### Ancestry Proportions (K=13)

![ADMIXTURE K=13 Ancestry Barplot](readme_images/linfantum_only_stacked_ancestry_barplot_k13.jpeg)

*Figure 3: Stacked ancestry barplot showing admixture proportions for each sample (K=13). Samples are ordered by dominant population assignment. Bar heights represent the proportion of ancestry from each of the 13 populations, revealing admixture patterns and population boundaries.*

---

### 2. Multi-Scale Phylogenetic Analysis

Comparing phylogenetic signal at three scales to detect selection on the MSL region:

- **Whole-genome phylogenies:** Neutral evolutionary baseline across all chromosomes.
- **100kb windows around MSL:** Regional selection signal on chromosome 31.
- **20kb windows around MSL:** Fine-scale locus-specific evolution.

Phylogenies stratified by ploidy level (diploid vs tetraploid strains) and annotated with geographic region, host species, and MSL copy number variation.

**Rationale:** If the MSL region is under selection, phylogenetic patterns on chromosome 31 will differ from genome-wide evolutionary history.

#### Whole-Genome Phylogeny (Neutral Baseline)

![Whole-Genome Phylogeny](readme_images/whole_genome_phylogeny_circular.jpeg)

*Figure 4: Circular maximum likelihood phylogenetic tree based on whole-genome SNPs (L. infantum only, n=463 samples). Tips are coloured by geographic region and annotated with host species (circles = human, triangles = dog) and MSL copy number (symbol size indicates ploidy). Clear Old World/New World split evident with strong bootstrap support.*

#### 100kb Window Around MSL (Regional Selection Signal)

![100kb Window Phylogeny](readme_images/100kb_window_phylogeny_circular.jpeg)

*Figure 5: Phylogenetic tree based on 100kb genomic window surrounding the MSL locus on chromosome 31. Comparison with whole-genome tree reveals topological differences in New World lineages, suggesting localized selection pressure on the MSL region in American populations where miltefosine resistance is prevalent.*

#### 20kb Window Around MSL (Fine-Scale Evolution)

![20kb Window Phylogeny](readme_images/20kb_window_phylogeny_circular.jpeg)

*Figure 6: High-resolution phylogeny from 20kb window directly flanking the MSL genes. Clustering patterns differ from genome-wide phylogeny, with MSL deletions (small symbols indicating low copy number) concentrated in specific Brazilian clades, providing evidence for selection acting on the MSL locus.*

---

### 3. Geographic Distribution Analysis ⏳

Spatial analysis of MSL copy number variation and population distribution across global sampling locations.

**Status:** Collaborative analysis complete, integration in progress.

---

### 4. Hardy-Weinberg Equilibrium Testing ⏳

Tests whether MSL allele frequencies deviate from neutral expectations within populations, distinguishing selection from demographic history.

**Status:** Analysis complete, integration in progress.

---

## Clinical Relevance

MSL deletions confer resistance to **miltefosine**, the only oral treatment for leishmaniasis. Understanding whether resistant strains spread via:

- **Selection:** Faster than expected, requires urgent intervention.
- **Drift:** Random spread, reflects demographic history.

...has direct implications for treatment policy and surveillance strategies.

---

## Methods Summary

**Population Structure:**
- ADMIXTURE v1.3 with K=1-20 tested.
- Cross-validation error minimisation for optimal K selection.
- High-confidence assignments: ≥90% ancestry threshold.

**Phylogenetics:**
- IQ-TREE2 maximum likelihood phylogenies with 1000 ultrafast bootstrap replicates.
- Model selection via ModelFinder, phylogenies visualized with ggtree (R).
- ggtree (R) for visualisation with geographic/ploidy annotation.

**Data:**
- 465 *L. infantum* isolates, global sampling.
- Whole-genome sequencing data (VCF format).
- Geographic metadata, MSL copy number (coverage-based).

---

## Repository Structure
```
├── data/
│   ├── metadata/              # Sample information, MSL coverage data
│   ├── population_structure/  # ADMIXTURE Q-matrices and CV results
│   └── phylogenetic/          # IQ-TREE output files (.treefile)
├── scripts/
│   ├── 01_population_structure_admixture.R
│   ├── Whole_Genome_Linfantum_With_Outgroup_Phylogenetics.R
│   ├── Whole_Genome_Linfantum_Only_Phylogenetics.R
│   ├── Whole_Genome_Americas_Only_Phylogenetics.R
│   ├── 100kb_Window_Linfantum_Only_Phylogenetics.R
│   ├── 100kb_Window_Americas_Only_Phylogenetics.R
│   ├── 20kb_Window_Linfantum_Only_Phylogenetics.R
│   └── 20kb_Window_Americas_Only_Phylogenetics.R
├── results/
│   ├── population_structure/  # Population assignments and summaries
│   └── tables/                # Summary statistics
├── figures/
│   ├── admixture/             # ADMIXTURE plots (all datasets)
│   └── phylogenetics/         # Phylogenetic tree figures (all scales)
└── readme_images/             # Optimized figures for README display
```

---

## Software & Dependencies

**R (v4.4.2):**
- tidyverse, ggtree, treeio, ape, phytools.
- RColorBrewer, ggtext.

**External:**
- ADMIXTURE v1.3.
- IQ-TREE2.
- vcftools (preprocessing).

---

## Data Sources

Genomic data: 465 *L. infantum* isolates from global sampling (Jeffares lab, University of York, unpublished).

---

## Project Status

**Completed:**
- ✅ ADMIXTURE population structure analysis (three datasets).
- ✅ Multi-scale phylogenetic analyses (whole-genome, 100kb, 20kb windows).
- ✅ Ploidy stratification and MSL-region specific analyses.

**In Progress:**
- ⏳ Integration of Hardy-Weinberg equilibrium testing.
- ⏳ Integration of geographic distribution analyses.
- ⏳ Manuscript preparation (collaborative, multi-author).

---

## Author

**Mathew Johansson**  
MSc Bioinformatics, University of York  
Voluntary Research Associate, Jeffares Lab  
GitHub: [@MathewJohansson](https://github.com/MathewJohansson)

---

## Acknowledgements

- Dr. Daniel Jeffares (University of York) - Supervision, data generation.
- Zeynep Sakaoglu - Geographical and HWE analyses collaboration.

---

## Licence

Code: MIT Licence (upon manuscript publication).  
Data: Available upon reasonable request pending publication.

---

## Key References

**ADMIXTURE:**  
Alexander, D.H., Novembre, J. and Lange, K. (2009). Fast model-based estimation of ancestry in unrelated individuals. *Genome Research*, 19(9), 1655-1664.

**IQ-TREE:**  
Nguyen, L.T., Schmidt, H.A., von Haeseler, A. and Minh, B.Q. (2015). IQ-TREE: a fast and effective stochastic algorithm for estimating maximum-likelihood phylogenies. *Molecular Biology and Evolution*, 32(1), 268-274.

---

**Note:** This is an active research project demonstrating bioinformatics workflow, reproducible research practices, and clinical/evolutionary genomics expertise. Analyses are being refined for manuscript preparation as a collaborative, multi-author publication.
