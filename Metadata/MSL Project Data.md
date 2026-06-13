

# Data for the L. infantum MSL project

**Author:** Daniel Jeffares  
**Date:** 2025-06-25  
**Summary:** This summarizes the data files for students

This file is best viewed in a plain text editor can display markdown format, like VSCode or even RStudio.



## Locations of data

- **Google Drive folder:** <https://drive.google.com/drive/folders/1KZEDF40uesLBgudGusjyI5Dnk_v20bN9?usp=drive_link>
- **Viking directory:** `/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN`
- **Local directory (Dan's Mac):** `/Users/dj757/gd/modules/BIO71M/BIO71M_2025_Linfantum_MSL_study`



## Summary of data

Within the Google Drive folder we have:

- **README_students.md:** this file
- **markdown_basics.md:** a simple description for how to write markdown format files (like this one)
- **phylogeny:** phylogeny data
- **population_structure:** population structure files
- **other_data:** sample metadata, country and Brazil state centroids, gene coverage data (all genes, all strains)
  - The file `vcf_samples29_location_and_MSL_data_2025-06-19.tsv` contains metadata and also MSL copy numbers.
- **plots:** plots
- **scripts:** scripts



## Sample information: strain names, source_country, species

The file `vcf_samples29_location_and_MSL_data_2025-06-19.tsv` describes all the samples.
It is a tab-separated file that can be read into R using `read_tsv`.

The most important fields are described below. This information can be used to annotate maps and phylogenies.

**Fields of the samples file:** `vcf_samples29_location_and_MSL_data_2025-06-19.tsv`

- **sra_run_accession:** the unique strain name, e.g. SRR20748484, SRR20748483
  - These actually refer to the unique ID of the genome sequence in the short read archive: <https://www.ncbi.nlm.nih.gov/sra/>
- **who_id:** the international code for each strain, e.g. MCAN/BR/1980/CR3, MHOM/ES/1992/LLM373
- **source_geographic_region:** which region of the world the strain came from: Americas, Central.East.Asia, East.Africa, Europe, Middle.East, North.Africa
- **source_country:** which country the strain came from: e.g. Brazil for MCAN/BR/1980/CR3, Spain for MHOM/ES/1992/LLM373
- **country_code:** A two letter code for the country, e.g. BR Brazil, ES Spain
- **source_state_or_region:** mainly for samples from Brazil, which state of Brazil the strain came from
- **host_species:** which species the strain came from, usually Human or dog
- **host_code:** a code for species MCAN dog, MHOM human
- **coverage:** the median coverage for the 4 genes in the MSL region (where coverage 2.0 = diploid)
- **coverage_rounded:** the closest integer of the coverage



## Phylogenies

Phylogenetic trees were run with various subsets of strains. Sequence files were generated from the VCF SNP files using vcf2phylyp. Phylogenies were generated using IQ-TREE2, which uses a maximum likelihood method.

The three subsets were:

- **Linfantum_samples_with_outgroup:** all samples. Best set to use so you can set the East African samples as the outgroup.
- **Linfantum_samples_only:** Excludes these two East African samples. Contains L. infantum from Europe, Brazil etc.
- **Americas_samples:** only samples from South and Central America. Mostly Brazil.

The most useful file within each directory is the `.treefile`. The other files are described at the end of the log file, e.g.
`Linfantum_samples_with_outgroup_2025-03-12_Fast_bootstrap_tree.log`

The `.treefile` can be read and plotted with ggtree like so:

```r
# load libraries
library(ggtree)
library(treeio)

# Read the tree file
tree <- read.tree("Linfantum_samples_with_outgroup_2025-03-12_Fast_bootstrap_tree.treefile")

# Or alternatively, if you want to preserve bootstrap values as node labels
tree <- read.newick("Linfantum_samples_with_outgroup_2025-03-12_Fast_bootstrap_tree.treefile")

# Create a basic tree plot
ggtree(tree)

# and so on
```



## Population structure information

To avoid all our analyses being subject to the Wahlund effect (reduction of heterozygotes due to population structure), we need to do our analysis within populations, not for all samples at once.

To enable this, I have used [Admixture](https://dalexander.github.io/admixture/index.html) to divide the samples into populations.

Admixture divides the samples into K populations (I ran K=1 to K=20 in our case). To determine what the most appropriate K value is we use the K with the lowest cross validation (CV) error.

I ran Admixture with K 1..20 with multiple subsets of the data.
The CV values for each run with different K values can be found here.
The K value(s) with the lowest K values are considered to explain the populations best.

**CV files:**

- Americas_samples_n388.CV.txt
- Linfantum_samples_with_outgroup_2025-03-12.CV.txt
- infantum_only_n463.CV.txt

The main population that each sample belongs to can be found in these files.
They also contain the proportion of each strain's genome that belongs to this population.
When working with within-population analysis, I suggest filtering for samples that have ≥ 90% of their genome belonging to a population. You will note that for the infantum_only set, I have given you two outputs with different K values.

**Population membership files:**

- population_membership_table.allsamples.k12.tsv
- population_membership_table.infantum_only.k8.tsv
- population_membership_table.infantum_only.k13.tsv
- population_membership_table.Americas_only.k13.tsv



## VCF files

These will be useful only when Viking is available.

- **Original VCF from the dating paper:** Ldon_Linf_qual30_MAC1_MQM20.SNPs.filt2.recode.vcf
  - From: `/mnt/scratch/projects/biol-leish-2019/PROJECTS/LDSC-DATING-2024/vcf-files/Ldon_Linf_qual30_MAC1_MQM20.SNPs.filt2.recode.vcf`

- **VCF with only those samples that we have information for** in `sampldating_paper_table_1_samples_2025-03-12.tsv`
  - `Ldon_Linf_qual30_MAC1_MQM20.SNPs.filt2.samples.2025-03-13.vcf.gz`
  - Contains the samples listed in `Ldon_Linf_samples.2025-03-13.txt`

- **VCF with only Linfantum and two outgroups:** `Linfantum_samples_with_outgroup_2025-03-12.vcf`
  - Samples are listed in `Linf_sample_names_with_outgroup_2025-03-13.txt`
  - Full information about these strains is in: `Linf_samples_with_outgroup_2025-03-13.tsv`
  - `Linfantum_samples_with_outgroup_2025-03-12.vcf.gz` (only Linfantum and two outgroups, for phylogeny)

**Other important files:**

- `Ldon_Linf_qual30_MAC1_MQM20.SNPs.filt2.samples.2025-03-13.vcf`
- Sample list: `dating_paper_table_1_samples_2025-03-12.tsv` (all sample information)
- Sample list: `Linfantum_samples_2025-03-12.tsv` (all samples in the VCF)
- `Linfantum_samples_with_outgroups_2025-03-12.tsv` (all Linfantum samples in the VCF)
- Reference genome and annotation: `reference_genome/TriTrypDB-66_LinfantumJPCM5_Genome_Plus_maxicircle.fasta`
  - This includes Leishmania infantum strain JPCM5 kinetoplast, complete sequence GenBank: BK010877.1
- Sample information for phylogeny: `Linf_samples_with_outgroup_2025-03-13.tsv`
- Reference genome chromosomes: `TriTrypDB-66_LinfantumJPCM5_Genome_Plus_maxicircle.chromosomes.bed`
- BAM files: `/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/bams`



## Summary of available data

1. Strain list, including sample names, SRA accessions
2. Genomic data in the form of VCF files (should correspond exactly to the strain list)
3. Phylogeny raw data, in the form of newick files
4. Locations: source country, and for Brazil the source state of all strains
5. BAM files with coverage information
6. MSL and chromosome depth for all strains (from the BAMs)



## Reference genome

**Reference genome:** `TriTrypDB-66_LinfantumJPCM5_Genome_Plus_maxicircle.fasta`

Location: `reference_genome/TriTrypDB-66_LinfantumJPCM5_Genome_Plus_maxicircle.fasta` and related files.



## Gene coverage data

The file `LdonLinf_all_samples_gene_coverages.tsv` has the coverage for each gene in 829 samples. Coverage data is also available in the Viking directories under:

- `/mnt/scratch/projects/biol-leish-2019/PROJECTS/LDSC-DATING-2024/Mario_masters_MSL_results`

Which contains coverage files:

- `LdonLinf_all_samples_gene_coverages.tsv`
- `LdonLinf_all_samples_MSL_and_31_coverage.txt`

BAM files are located at:

- `/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/bams`

We have all 899 BAM files for samples in the dating paper VCF and all 465 BAM files for samples in the Linfantum+outgroups VCF.




