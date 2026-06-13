# R Package Requirements for Leishmania MSL Population Genomics
# Complete list of dependencies used in this analysis pipeline

# Core Analysis & Population Genetics
library(adegenet)              # Multivariate analysis of genetic markers
library(ape)                   # Phylogenetic analysis and Neighbour-Joining trees
library(pegas)                 # Population genetics
library(StAMPP)                # Population differentiation (Fst, etc.)
library(phytools)              # Phylogenetic comparative methods
library(treeio)                # Tree file I/O
library(vcfR)                  # VCF file handling and manipulation

# Data Manipulation & I/O
library(tidyverse)             # Data wrangling workflow (includes: dplyr, ggplot2, readr, tidyr, stringr, forcats, tibble)
library(dplyr)                 # Data manipulation
library(readr)                 # Efficient file reading
library(stringr)               # String manipulation
library(tibble)                # Enhanced data frames
library(tidyr)                 # Data reshaping
library(forcats)               # Factor manipulation

# Visualisation & Graphics
library(ggplot2)               # Grammar of graphics
library(ggtree)                # Phylogenetic tree visualization
library(ggtext)                # Enhanced text rendering for ggplot2
library(ggpubr)                # Publication-ready plots
library(ggforce)               # Additional geometric layers
library(ggnewscale)            # Multiple aesthetic scales
library(cowplot)               # Plot composition and alignment
library(gridExtra)             # Grid-based plot layouts
library(patchwork)             # Intuitive plot composition
library(RColorBrewer)          # Colour palettes
library(viridis)               # Colourblind-friendly colour scales
library(colorspace)            # Colour space manipulation
library(Polychrome)            # Accessible colour palettes for many categories

# Geographic & Spatial
library(rnaturalearth)         # Natural Earth geographic data
library(rnaturalearthdata)     # Data for rnaturalearth
library(sf)                    # Spatial features (vector data)

# Utilities
library(gtools)                # Miscellaneous tools (mixedsort for natural sorting)
library(argparse)              # Command-line argument parsing
