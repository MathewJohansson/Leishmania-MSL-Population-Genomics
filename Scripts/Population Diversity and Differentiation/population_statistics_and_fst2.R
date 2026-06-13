#!/usr/bin/env Rscript

# population_statistics_and_fst2.R
# Fast summary pipeline from a split run-table:
# - mean pairwise FST for listed population pairs only
# - mean nucleotide diversity (pi) and Tajima's D for non-redundant populations
#
# Input table: whitespace-delimited, uses only first 3 columns:
#   col1 = pop1 sample file
#   col2 = pop2 sample file
#   col3 = input bgzipped VCF
#
# Outputs:
#   PREFIX.summary.tsv
#   PREFIX.run.log

usage <- function() {
  cat(
"Usage:
  Rscript population_statistics_and_fst2.R \\
    --input-table toy6.split.smcruntable.txt \\
    --output-prefix toy6.stats

Optional arguments:
  --stat-window-size INT   Window size for pi and Tajima's D (default: 10000)
  --fst-window-size INT    Window size for FST (default: 1000)
  --fst-window-step INT    Step size for FST (default: 1000)
  --work-dir DIR           Directory for intermediate vcftools files (default: PREFIX.tmp)
  --force                  Re-run vcftools even if output files already exist
  --help                   Show this message
\n"
  )
}

parse_args <- function(args) {
  opt <- list(
    input_table = NULL,
    output_prefix = NULL,
    stat_window_size = 10000L,
    fst_window_size = 1000L,
    fst_window_step = 1000L,
    work_dir = NULL,
    force = FALSE,
    help = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    key <- args[[i]]
    if (key %in% c("--help", "-h")) {
      opt$help <- TRUE
      i <- i + 1L
      next
    }
    if (key == "--force") {
      opt$force <- TRUE
      i <- i + 1L
      next
    }

    if (i == length(args)) {
      stop("Missing value for argument: ", key)
    }
    value <- args[[i + 1L]]

    if (key == "--input-table") {
      opt$input_table <- value
    } else if (key == "--output-prefix") {
      opt$output_prefix <- value
    } else if (key == "--stat-window-size") {
      opt$stat_window_size <- as.integer(value)
    } else if (key == "--fst-window-size") {
      opt$fst_window_size <- as.integer(value)
    } else if (key == "--fst-window-step") {
      opt$fst_window_step <- as.integer(value)
    } else if (key == "--work-dir") {
      opt$work_dir <- value
    } else {
      stop("Unknown argument: ", key)
    }
    i <- i + 2L
  }

  opt
}

safe_id <- function(x) {
  y <- gsub("[^A-Za-z0-9._-]", "_", x)
  y <- gsub("_+", "_", y)
  y
}

read_numeric_col_mean <- function(path, colname) {
  if (!file.exists(path)) {
    return(NA_real_)
  }
  dat <- tryCatch(
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = ""),
    error = function(e) NULL
  )
  if (is.null(dat) || !(colname %in% names(dat)) || nrow(dat) == 0) {
    return(NA_real_)
  }
  vals <- suppressWarnings(as.numeric(dat[[colname]]))
  if (all(is.na(vals))) {
    return(NA_real_)
  }
  mean(vals, na.rm = TRUE)
}

make_union_keep <- function(pop1_file, pop2_file, out_file) {
  s1 <- readLines(pop1_file, warn = FALSE)
  s2 <- readLines(pop2_file, warn = FALSE)
  s <- unique(trimws(c(s1, s2)))
  s <- s[nzchar(s)]
  writeLines(s, out_file)
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  opt <- parse_args(args)

  if (isTRUE(opt$help)) {
    usage()
    quit(status = 0)
  }

  if (is.null(opt$input_table) || is.null(opt$output_prefix)) {
    usage()
    stop("--input-table and --output-prefix are required")
  }
  if (!file.exists(opt$input_table)) {
    stop("Input table not found: ", opt$input_table)
  }
  if (!nzchar(Sys.which("vcftools"))) {
    stop("vcftools not found in PATH")
  }

  if (is.null(opt$work_dir)) {
    opt$work_dir <- paste0(opt$output_prefix, ".tmp")
  }
  dir.create(opt$work_dir, showWarnings = FALSE, recursive = TRUE)

  log_file <- paste0(opt$output_prefix, ".run.log")
  summary_file <- paste0(opt$output_prefix, ".summary.tsv")

  log_con <- file(log_file, open = "wt")
  on.exit(close(log_con), add = TRUE)

  log_msg <- function(...) {
    msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste(..., collapse = ""))
    writeLines(msg, con = log_con)
    message(msg)
  }

  run_vcftools <- function(args_vec, desc) {
    cmd_txt <- paste(c("vcftools", shQuote(args_vec)), collapse = " ")
    log_msg("START ", desc, " | ", cmd_txt)
    out <- tryCatch(
      system2("vcftools", args = args_vec, stdout = TRUE, stderr = TRUE),
      error = function(e) structure(paste("SYSTEM2_ERROR:", conditionMessage(e)), status = 1L)
    )
    status <- attr(out, "status")
    if (is.null(status)) {
      status <- 0L
    }
    if (length(out) > 0) {
      writeLines(paste0("    ", out), con = log_con)
    }
    log_msg("END   ", desc, " | exit_status=", status)
    status
  }

  raw_lines <- readLines(opt$input_table, warn = FALSE)
  raw_lines <- trimws(raw_lines)
  raw_lines <- raw_lines[nzchar(raw_lines)]
  raw_lines <- raw_lines[!grepl("^#", raw_lines)]
  if (length(raw_lines) == 0) {
    stop("No non-empty, non-comment lines found in input table")
  }

  parsed <- lapply(raw_lines, function(line) {
    fields <- strsplit(line, "\\s+")[[1]]
    if (length(fields) < 3) {
      return(NULL)
    }
    data.frame(pop1 = fields[1], pop2 = fields[2], vcf = fields[3], stringsAsFactors = FALSE)
  })
  parsed <- parsed[!vapply(parsed, is.null, logical(1))]
  if (length(parsed) == 0) {
    stop("No valid rows with at least 3 columns found in input table")
  }

  input_df <- do.call(rbind, parsed)
  input_df$row_id <- seq_len(nrow(input_df))

  missing_pop1 <- input_df$pop1[!file.exists(input_df$pop1)]
  missing_pop2 <- input_df$pop2[!file.exists(input_df$pop2)]
  missing_vcf <- input_df$vcf[!file.exists(input_df$vcf)]
  if (length(missing_pop1) + length(missing_pop2) + length(missing_vcf) > 0) {
    stop(
      "Missing files detected. pop1 missing: ", paste(unique(missing_pop1), collapse = ", "),
      " | pop2 missing: ", paste(unique(missing_pop2), collapse = ", "),
      " | vcf missing: ", paste(unique(missing_vcf), collapse = ", ")
    )
  }

  # Deduplicate triplets for running, while preserving original rows for final output.
  unique_triplets <- unique(input_df[, c("pop1", "pop2", "vcf")])
  unique_triplets$triplet_id <- seq_len(nrow(unique_triplets))

  # Build unique population list per VCF for pi/TajimaD.
  pop_pairs <- unique(rbind(
    data.frame(vcf = unique_triplets$vcf, pop = unique_triplets$pop1, stringsAsFactors = FALSE),
    data.frame(vcf = unique_triplets$vcf, pop = unique_triplets$pop2, stringsAsFactors = FALSE)
  ))
  pop_pairs$pop_id <- seq_len(nrow(pop_pairs))

  log_msg("Input rows: ", nrow(input_df), "; unique pair+vcf triplets: ", nrow(unique_triplets), "; unique population+vcf items: ", nrow(pop_pairs))

  # Cache for population stats (mean pi, mean TajimaD) keyed by vcf||pop
  pop_cache <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(pop_pairs))) {
    vcf <- pop_pairs$vcf[i]
    pop <- pop_pairs$pop[i]
    pop_key <- paste(vcf, pop, sep = "||")

    file_stub <- paste0(
      "pop_", pop_pairs$pop_id[i], "__",
      safe_id(basename(pop)), "__",
      safe_id(basename(vcf))
    )
    out_prefix <- file.path(opt$work_dir, file_stub)
    pi_file <- paste0(out_prefix, ".windowed.pi")
    taj_file <- paste0(out_prefix, ".Tajima.D")

    if (isTRUE(opt$force) || !file.exists(pi_file)) {
      args_pi <- c("--gzvcf", vcf, "--keep", pop, "--window-pi", as.character(opt$stat_window_size), "--out", out_prefix)
      status_pi <- run_vcftools(args_pi, paste0("pi for ", basename(pop), " on ", basename(vcf)))
      if (status_pi != 0L) {
        log_msg("WARN  pi failed for ", pop_key)
      }
    }

    if (isTRUE(opt$force) || !file.exists(taj_file)) {
      args_taj <- c("--gzvcf", vcf, "--keep", pop, "--TajimaD", as.character(opt$stat_window_size), "--out", out_prefix)
      status_taj <- run_vcftools(args_taj, paste0("TajimaD for ", basename(pop), " on ", basename(vcf)))
      if (status_taj != 0L) {
        log_msg("WARN  TajimaD failed for ", pop_key)
      }
    }

    pop_cache[[pop_key]] <- list(
      mean_pi = read_numeric_col_mean(pi_file, "PI"),
      mean_tajimaD = read_numeric_col_mean(taj_file, "TajimaD")
    )
  }

  # FST for listed pairs only.
  fst_cache <- new.env(parent = emptyenv())

  for (i in seq_len(nrow(unique_triplets))) {
    pop1 <- unique_triplets$pop1[i]
    pop2 <- unique_triplets$pop2[i]
    vcf <- unique_triplets$vcf[i]
    triplet_key <- paste(pop1, pop2, vcf, sep = "||")

    file_stub <- paste0(
      "fst_", unique_triplets$triplet_id[i], "__",
      safe_id(basename(pop1)), "__",
      safe_id(basename(pop2)), "__",
      safe_id(basename(vcf))
    )
    out_prefix <- file.path(opt$work_dir, file_stub)
    fst_file <- paste0(out_prefix, ".windowed.weir.fst")

    # Make a pair-specific keep file (union of pop1 and pop2 samples), per user request.
    pair_keep_file <- paste0(out_prefix, ".pair.keep.txt")
    if (isTRUE(opt$force) || !file.exists(pair_keep_file)) {
      make_union_keep(pop1, pop2, pair_keep_file)
    }

    if (isTRUE(opt$force) || !file.exists(fst_file)) {
      args_fst <- c(
        "--gzvcf", vcf,
        "--keep", pair_keep_file,
        "--weir-fst-pop", pop1,
        "--weir-fst-pop", pop2,
        "--fst-window-size", as.character(opt$fst_window_size),
        "--fst-window-step", as.character(opt$fst_window_step),
        "--out", out_prefix
      )
      status_fst <- run_vcftools(args_fst, paste0("FST for ", basename(pop1), " vs ", basename(pop2), " on ", basename(vcf)))
      if (status_fst != 0L) {
        log_msg("WARN  FST failed for ", triplet_key)
      }
    }

    fst_cache[[triplet_key]] <- read_numeric_col_mean(fst_file, "WEIGHTED_FST")
  }

  # Build summary table in original input order.
  summary_df <- data.frame(
    row_id = input_df$row_id,
    pop1_file = input_df$pop1,
    pop2_file = input_df$pop2,
    vcf_file = input_df$vcf,
    mean_weighted_fst = NA_real_,
    mean_pi_pop1 = NA_real_,
    mean_pi_pop2 = NA_real_,
    mean_tajimaD_pop1 = NA_real_,
    mean_tajimaD_pop2 = NA_real_,
    status = "OK",
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(summary_df))) {
    pop1 <- summary_df$pop1_file[i]
    pop2 <- summary_df$pop2_file[i]
    vcf <- summary_df$vcf_file[i]

    triplet_key <- paste(pop1, pop2, vcf, sep = "||")
    key1 <- paste(vcf, pop1, sep = "||")
    key2 <- paste(vcf, pop2, sep = "||")

    summary_df$mean_weighted_fst[i] <- fst_cache[[triplet_key]]

    if (!is.null(pop_cache[[key1]])) {
      summary_df$mean_pi_pop1[i] <- pop_cache[[key1]]$mean_pi
      summary_df$mean_tajimaD_pop1[i] <- pop_cache[[key1]]$mean_tajimaD
    }
    if (!is.null(pop_cache[[key2]])) {
      summary_df$mean_pi_pop2[i] <- pop_cache[[key2]]$mean_pi
      summary_df$mean_tajimaD_pop2[i] <- pop_cache[[key2]]$mean_tajimaD
    }

    if (is.na(summary_df$mean_weighted_fst[i]) ||
        is.na(summary_df$mean_pi_pop1[i]) ||
        is.na(summary_df$mean_pi_pop2[i]) ||
        is.na(summary_df$mean_tajimaD_pop1[i]) ||
        is.na(summary_df$mean_tajimaD_pop2[i])) {
      summary_df$status[i] <- "WARN_NA"
    }
  }

  write.table(summary_df, file = summary_file, sep = "\t", row.names = FALSE, quote = FALSE)
  log_msg("Wrote summary table: ", summary_file)
  log_msg("Run complete")

  cat("Done\n")
  cat("Summary:", summary_file, "\n")
  cat("Log:", log_file, "\n")
}

main()


