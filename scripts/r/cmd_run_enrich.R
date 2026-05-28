#!/usr/bin/env Rscript

# 役割:
#   DESeq2 stats から ORA/GSEA enrichment を独立実行する。
# 入力:
#   stats CSV、genome、出力先、各種 enrichment オプション。
# 出力:
#   method/contrast ごとの enrichment table と plot。

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
analysis_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)
source(file.path(script_dir, "analysis_plot_spec.R"))
source(file.path(script_dir, "analysis_runtime.R"))
paths <- pic_runtime_paths(analysis_root)
cfg <- pic_runtime_defaults()
genome_map_path <- paths$genome_map_file
enrichment_lib_dir <- paths$enrichment_lib_dir
biomart_lookup_file <- paths$biomart_lookup_file
biomart_cache_dir <- paths$biomart_cache_dir
biomart_ortholog_dir <- paths$biomart_ortholog_dir

source(file.path(script_dir, "deseq2_data.R"))
source(file.path(script_dir, "plot_common.R"))
source(file.path(script_dir, "enrichment.R"))
source(file.path(script_dir, "prepare_enrichment_libs.R"))

parse_enrichment_args <- function(args) {
  defaults <- list(
    stats = NULL,
    genome = NULL,
    out_dir = NULL,
    project_name = NULL,
    threads = cfg$threads,
    methods = cfg$enrich_methods_csv,
    deg_clusters = NULL,
    prepare_libs = FALSE,
    categories = NULL
  )
  provided_keys <- character()

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }

    if (identical(arg, "--help")) {
      help_file <- file.path(paths$help_dir, "help_enrich.txt")
      cat(readChar(help_file, file.info(help_file)$size, useBytes = TRUE))
      cat("\n")
      quit(save = "no", status = 0)
    }

    if (identical(arg, "--prepare-libs")) {
      defaults$prepare_libs <- TRUE
      i <- i + 1
      next
    }

    if (grepl("=", arg, fixed = TRUE)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=(.*)$", "\\1", arg)
    } else {
      key <- sub("^--", "", arg)
      i <- i + 1
      if (i > length(args)) {
        stop(sprintf("Option requires a value: --%s", key), call. = FALSE)
      }
      value <- args[[i]]
    }

    key <- gsub("-", "_", key)
    if (!key %in% names(defaults)) {
      stop(sprintf("Unknown option: --%s", gsub("_", "-", key)), call. = FALSE)
    }

    provided_keys <- c(provided_keys, key)
    defaults[[key]] <- value
    i <- i + 1
  }

  if (isTRUE(defaults$prepare_libs)) {
    invalid <- setdiff(unique(provided_keys), c("genome", "categories"))
    if (length(invalid) > 0) {
      stop("--prepare-libs accepts only --genome and --categories.", call. = FALSE)
    }
    if (is.null(defaults$genome) || defaults$genome == "") {
      stop("--prepare-libs requires --genome <name>.", call. = FALSE)
    }
    if (!is.null(defaults$categories) && defaults$categories != "") {
      req <- strsplit(as.character(defaults$categories), ",", fixed = TRUE)[[1]]
      req <- toupper(trimws(req))
      req <- req[req != ""]
      if (length(req) == 0) {
        stop("--categories must include at least one category.", call. = FALSE)
      }
      known <- toupper(as.character(enrichment_resources))
      unknown <- setdiff(req, known)
      if (length(unknown) > 0) {
        stop(
          sprintf(
            "Unknown --categories: %s (known: %s)",
            paste(unknown, collapse = ","),
            paste(known, collapse = ",")
          ),
          call. = FALSE
        )
      }
      defaults$categories <- req
    } else {
      defaults$categories <- NULL
    }
    defaults$genome <- as.character(defaults$genome)
    return(defaults)
  }

  if (is.null(defaults$stats) || is.null(defaults$genome) || is.null(defaults$out_dir)) {
    stop("Missing required options: --stats, --genome, --out-dir", call. = FALSE)
  }

  defaults$stats <- normalizePath(defaults$stats, winslash = "/", mustWork = TRUE)
  defaults$out_dir <- normalizePath(defaults$out_dir, winslash = "/", mustWork = FALSE)
  if (!is.null(defaults$deg_clusters) && defaults$deg_clusters != "") {
    defaults$deg_clusters <- normalizePath(defaults$deg_clusters, winslash = "/", mustWork = FALSE)
  } else {
    defaults$deg_clusters <- NULL
  }
  defaults$genome <- as.character(defaults$genome)

  defaults$threads <- as.integer(defaults$threads)
  if (is.na(defaults$threads) || defaults$threads < 1) {
    stop("--threads must be an integer >= 1.", call. = FALSE)
  }

  if (is.null(defaults$project_name) || defaults$project_name == "") {
    defaults$project_name <- basename(defaults$stats)
    defaults$project_name <- sub("\\.csv$", "", defaults$project_name)
    defaults$project_name <- sub("^stats_", "", defaults$project_name)
  }

  defaults$methods <- strsplit(as.character(defaults$methods), ",", fixed = TRUE)[[1]]
  defaults
}

main <- function() {
  args <- parse_enrichment_args(commandArgs(trailingOnly = TRUE))

  if (isTRUE(args$prepare_libs)) {
    message(sprintf("[INFO] Preparing local enrichment gson libraries for genome: %s", args$genome))
    if (is.null(args$categories)) {
      message("[INFO] Updating categories: all")
    } else {
      message(sprintf("[INFO] Updating categories: %s", paste(args$categories, collapse = ",")))
    }
    prepare_enrichment_libs(target_genome = args$genome, target_resources = args$categories)
    return(invisible(NULL))
  }

  message(sprintf("[INFO] Reading stats table: %s", args$stats))
  stats <- suppressMessages(readr::read_csv(args$stats, show_col_types = FALSE, progress = FALSE))
  dir.create(args$out_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("[INFO] Running enrichment for genome: %s", args$genome))
  run_enrichment(
    stats = stats,
    genome = args$genome,
    out_dir = args$out_dir,
    project_name = args$project_name,
    methods = args$methods,
    enrichment_lib_dir = enrichment_lib_dir,
    genome_map_file = genome_map_path,
    biomart_lookup_file = biomart_lookup_file,
    biomart_cache_dir = biomart_cache_dir,
    biomart_ortholog_dir = biomart_ortholog_dir,
    deg_cluster_file = args$deg_clusters,
    parallel_threads = args$threads
  )
  message("[INFO] Enrichment workflow completed")
}

main()
