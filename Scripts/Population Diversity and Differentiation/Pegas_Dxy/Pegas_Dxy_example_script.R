# SUMMARY -----------------------------------------------------------------------
# This R script screens a set of FASTA files prior to temporal signal analysis.
# It calculates population genetic metrics and writes results in batches of 100.
#
# Metrics computed per region:
#   1. n_seg       : number of segregating sites
#   2. pi          : nucleotide diversity across all samples
#   3. dxy         : mean absolute divergence between two populations
#   4. dxy_pi      : dxy / pi ratio
#   5. mean_r2     : mean pairwise r^2 LD across all SNP pairs (pegas::LDscan)
#   6. max_r2      : max pairwise r^2 LD across all SNP pairs
#
# Usage (command line):
#   Rscript temporal_signal_prelim_metrics3.R \
#     --fasta_dir   /path/to/fastas \
#     --meta_path   /path/to/metadata.tsv \
#     --output_dir  /path/to/output \
#     --basal_pop1  "PopulationA" \
#     --basal_pop2  "PopulationB" \
#     --min_seg     20

# LOAD LIBRARIES ----------------------------------------------------------------
suppressPackageStartupMessages({
	library(argparse)
	library(ape)
	library(pegas)
	library(tidyverse)
})

# PARSE COMMAND-LINE ARGUMENTS --------------------------------------------------
parse_cli_args <- function() {
	parser <- ArgumentParser(
		description = "Screen genomic regions by population genetic metrics."
	)
	
	parser$add_argument("--fasta_dir",  type="character", required=TRUE)
	parser$add_argument("--meta_path",  type="character", required=TRUE)
	parser$add_argument("--output_dir", type="character", required=TRUE)
	parser$add_argument("--basal_pop1", type="character", required=TRUE)
	parser$add_argument("--basal_pop2", type="character", required=TRUE)
	parser$add_argument("--pop_col",    type="character", default="Population")
	parser$add_argument("--id_col",     type="character", default="SRA.Run.Accession")
	parser$add_argument("--min_seg",    type="integer",   default=20L)
	parser$add_argument("--dxy_model",  type="character", default="raw")
	
	args <- parser$parse_args()
	return(args)
}

args <- parse_cli_args()

# IDENTIFY FILES ----------------------------------------------------------------
fasta_files <- list.files(
	args$fasta_dir,
	pattern     = "\\.(fa|fasta|fna|fas)$",
	full.names  = TRUE,
	ignore.case = TRUE
)

meta     <- read_tsv(args$meta_path, show_col_types = FALSE)
pop1_ids <- meta[[args$id_col]][meta[[args$pop_col]] == args$basal_pop1]
pop2_ids <- meta[[args$id_col]][meta[[args$pop_col]] == args$basal_pop2]

output_dir <- args$output_dir
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# LOGGING SETUP -----------------------------------------------------------------
log_path <- file.path(output_dir, paste0("screener_", format(Sys.time(), "%Y-%m-%dT%H%M%S"), ".log"))
log_con  <- file(log_path, open = "wt")
sink(log_con, split = TRUE)

log_msg <- function(...) {
	cat(paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", ..., "\n"), sep = "")
}

# STARTUP LOG -------------------------------------------------------------------
log_msg("========================================")
log_msg(paste0("Script:      temporal_signal_prelim_metrics3.R"))
log_msg(paste0("Started:     ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
log_msg(paste0("fasta_dir:   ", args$fasta_dir))
log_msg(paste0("meta_path:   ", args$meta_path))
log_msg(paste0("output_dir:  ", args$output_dir))
log_msg(paste0("basal_pop1:  ", args$basal_pop1, "  (n=", length(pop1_ids), " samples)"))
log_msg(paste0("basal_pop2:  ", args$basal_pop2, "  (n=", length(pop2_ids), " samples)"))
log_msg(paste0("min_seg:     ", args$min_seg))
log_msg(paste0("dxy_model:   ", args$dxy_model))
log_msg(paste0("n_regions:   ", length(fasta_files)))
log_msg("========================================")

if (length(pop1_ids)  == 0) stop(sprintf("No samples found for basal_pop1: '%s'", args$basal_pop1), call. = FALSE)
if (length(pop2_ids)  == 0) stop(sprintf("No samples found for basal_pop2: '%s'", args$basal_pop2), call. = FALSE)
if (length(fasta_files) == 0) stop(sprintf("No FASTA files found in: %s", args$fasta_dir), call. = FALSE)

# PREPARE OUTPUT FILE -----------------------------------------------------------
tsv_path <- file.path(
	output_dir,
	paste0("region_screen_", format(Sys.time(), "%Y-%m-%dT%H%M%S"), ".tsv")
)

header <- tibble(
	region    = character(),
	n_samples = integer(),
	n_sites   = integer(),
	n_seg     = integer(),
	pi        = numeric(),
	dxy       = numeric(),
	dxy_pi    = numeric(),
	mean_r2   = numeric(),   # mean pairwise r^2 across all SNP pairs
	max_r2    = numeric()    # max pairwise r^2 across all SNP pairs
)
write_tsv(header, tsv_path)

# BESPOKE FUNCTIONS -------------------------------------------------------------

#' Load FASTA as DNAbin, subsetting to metadata samples only.
#' Returns NULL if unreadable or fewer than 2 matched samples.
load_fasta_as_dnabin <- function(fasta_path, target_ids) {
	aln  <- tryCatch(read.dna(fasta_path, format = "fasta"), error = function(e) NULL)
	if (is.null(aln)) return(NULL)
	keep <- intersect(rownames(aln), target_ids)
	if (length(keep) < 2) return(NULL)
	return(aln[keep, , drop = FALSE])
}

#' Compute Dxy — mean absolute pairwise divergence between two populations.
compute_dxy <- function(aln, pop1_ids, pop2_ids, model = "raw") {
	p1 <- intersect(rownames(aln), pop1_ids)
	p2 <- intersect(rownames(aln), pop2_ids)
	if (length(p1) == 0 || length(p2) == 0) return(NA_real_)
	combined <- aln[c(p1, p2), , drop = FALSE]
	d  <- as.matrix(dist.dna(combined, model = model, pairwise.deletion = TRUE))
	n1 <- length(p1)
	n2 <- length(p2)
	between <- d[seq_len(n1), seq_len(n2) + n1, drop = FALSE]
	mean(between, na.rm = TRUE)
}

#' Compute mean and max r^2 LD across all SNP pairs in the alignment.
#' Uses pegas::LDscan() which returns pairwise r values (not r^2).
#' Requires at least 2 segregating sites; returns NA if LDscan fails.
#'
#' Note: high mean_r2 indicates strong linkage across the region — useful as
#' a flag that recombination has NOT broken down LD, consistent with a single
#' non-recombining locus suitable for phylogenetic/clock analysis.
#' Low mean_r2 may indicate residual recombination or mutational saturation.
compute_ld <- function(aln, n_seg) {
	if (n_seg < 2) return(list(mean_r2 = NA_real_, max_r2 = NA_real_))
	r_vals <- tryCatch({
		capture.output(
			r <- suppressMessages(LDscan(aln, what = "r")),
			type = "output"
		)
		as.numeric(r)
	},
	error = function(e) NULL
	)
	if (is.null(r_vals) || length(r_vals) == 0) {
		return(list(mean_r2 = NA_real_, max_r2 = NA_real_))
	}
	r2_vals <- r_vals^2
	list(
		mean_r2 = mean(r2_vals, na.rm = TRUE),
		max_r2  = max(r2_vals,  na.rm = TRUE)
	)
}


# MAIN LOOP ---------------------------------------------------------------------
run_start   <- proc.time()
buffer_list <- list()
batch_size  <- 100L
n_succeeded <- 0L
n_skipped   <- 0L
all_ids     <- meta[[args$id_col]]

for (i in seq_along(fasta_files)) {
	
	input_fasta     <- fasta_files[[i]]
	alignment_label <- tools::file_path_sans_ext(basename(input_fasta))
	
	# Load alignment
	aln <- load_fasta_as_dnabin(input_fasta, target_ids = all_ids)
	if (is.null(aln)) { n_skipped <- n_skipped + 1L; next }
	
	# Segregating sites pre-filter
	n_seg <- length(seg.sites(aln))
	if (n_seg < args$min_seg) { n_skipped <- n_skipped + 1L; next }
	
	# Compute metrics
	pi_val       <- nuc.div(aln)
	dxy_val      <- compute_dxy(aln, pop1_ids, pop2_ids, model = args$dxy_model)
	dxy_pi_ratio <- if (!is.na(dxy_val) && !is.na(pi_val) && pi_val > 0) dxy_val / pi_val else NA_real_
	ld           <- compute_ld(aln, n_seg)
	
	# Add to buffer
	buffer_list[[length(buffer_list) + 1]] <- tibble(
		region    = alignment_label,
		n_samples = nrow(aln),
		n_sites   = ncol(aln),
		n_seg     = n_seg,
		pi        = pi_val,
		dxy       = dxy_val,
		dxy_pi    = dxy_pi_ratio,
		mean_r2   = ld$mean_r2,
		max_r2    = ld$max_r2
	)
	
	# Batch write
	if (length(buffer_list) >= batch_size || i == length(fasta_files)) {
		if (length(buffer_list) > 0) {
			batch_df    <- bind_rows(buffer_list)
			write_tsv(batch_df, tsv_path, append = TRUE)
			n_succeeded <- n_succeeded + nrow(batch_df)
			buffer_list <- list()
			log_msg(sprintf("Processed %d/%d: Batch sync complete.", i, length(fasta_files)))
		}
	}
}

# SORT OUTPUT -------------------------------------------------------------------
final_df <- read_tsv(tsv_path, show_col_types = FALSE) |>
	arrange(desc(dxy), desc(dxy_pi))
write_tsv(final_df, tsv_path)
log_msg("Output sorted by dxy desc, dxy_pi desc.")

# CLEANUP -----------------------------------------------------------------------
elapsed <- round((proc.time() - run_start)["elapsed"])
log_msg("========================================")
log_msg(sprintf("Run complete. Total succeeded: %d, Skipped: %d", n_succeeded, n_skipped))
log_msg(sprintf("Results saved to: %s", tsv_path))
log_msg(sprintf("Elapsed: %d seconds", elapsed))
log_msg("========================================")

sink()
close(log_con)