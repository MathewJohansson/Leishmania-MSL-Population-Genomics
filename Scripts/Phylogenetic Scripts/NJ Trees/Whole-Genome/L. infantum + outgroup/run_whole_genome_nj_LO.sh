#!/usr/bin/env bash
#SBATCH --time=12:00:00
#SBATCH --error=whole_genome_nj_LO_%A.e
#SBATCH --output=whole_genome_nj_LO_%A.log
#SBATCH --account=biol-leish-2019
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G

module load R/4.4.1-gfbf-2023b

cd /mnt/scratch/projects/biol-leish-2019/PROJECTS/2025-LINFANTUM-POPGEN

Rscript scripting/tetraploid_nei_genetic_distance_whole_genome_Linfantum_plus_outgroup.R
