#!/bin/bash
#SBATCH --job-name=unified_nonchr31
#SBATCH --output=/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/log_files/unified_nonchr31_%j.out
#SBATCH --error=/mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/log_files/unified_nonchr31_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=bgv512@york.ac.uk

# Load R
module load R/4.4.1-gfbf-2023b

# Create output directory if it doesn't exist
mkdir -p /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/popgen_unified

# Run script
Rscript /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/scripting/Unified_Population_Statistics_NonChr31.R \
--vcf_path  /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/vcf_files/Linfantum_samples_no_outgroup_biallelic.2025-06-10.noLinJ.31.SNPsonly.vcf.gz \
--meta_path /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/scripting/vcf_samples29_location_and_MSL_data_2025-06-19_corrected_with_population_info.tsv \
--mem_path  /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/scripting/population_membership_table_Americas_only_k13_AM_notation.tsv \
--out_dir   /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN/popgen_unified