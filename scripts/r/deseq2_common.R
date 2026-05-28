required_packages <- c(
  "readr", "dplyr", "tidyr", "purrr", "tibble",
  "DESeq2", "ggplot2", "ggrepel", "ggsci", "ggbeeswarm", "scales", "pheatmap",
  "clusterProfiler"
)

load_required_packages <- function() {
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Missing required R packages: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(lapply(required_packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

parse_deseq2_args <- function(args) {
  cfg <- pic_runtime_defaults()
  defaults <- list(
    deftable = NULL,
    count_dir = NULL,
    genome = NULL,
    out_dir = NULL,
    fdr = cfg$fdr
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }

    if (identical(arg, "--help")) {
      help_file <- file.path(analysis_root, "help", "help_deseq2.txt")
      cat(readChar(help_file, file.info(help_file)$size, useBytes = TRUE))
      cat("\n")

      registered <- "  (none)"
      if (exists("biomart_lookup_file", inherits = TRUE) &&
          file.exists(get("biomart_lookup_file", inherits = TRUE))) {
        lookup <- suppressMessages(readr::read_tsv(
          get("biomart_lookup_file", inherits = TRUE),
          show_col_types = FALSE,
          progress = FALSE,
          na = c("", "NA")
        ))
        for (nm in colnames(lookup)) {
          lookup[[nm]] <- gsub("\r", "", as.character(lookup[[nm]]), fixed = TRUE)
        }
        if (nrow(lookup) > 0 && "genome" %in% colnames(lookup) && "biomart_dataset" %in% colnames(lookup)) {
          registered_tbl <- lookup |>
            dplyr::filter(!is.na(.data$genome), .data$genome != "", !is.na(.data$biomart_dataset), .data$biomart_dataset != "") |>
            dplyr::distinct(.data$genome, .data$biomart_dataset) |>
            dplyr::arrange(.data$genome)
          if (nrow(registered_tbl) > 0) {
            registered <- paste(sprintf("  - %s (%s)", registered_tbl$genome, registered_tbl$biomart_dataset), collapse = "\n")
          }
        }
      }

      cat(
        paste(
          "Registered biomart genomes:",
          registered,
          sep = "\n"
        ),
        "\n"
      )
      quit(save = "no", status = 0)
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

    defaults[[key]] <- value
    i <- i + 1
  }

  if (is.null(defaults$deftable) || is.null(defaults$count_dir) || is.null(defaults$genome) || is.null(defaults$out_dir)) {
    stop("Missing required options: --deftable, --count-dir, --genome, --out-dir", call. = FALSE)
  }

  defaults$deftable <- normalizePath(defaults$deftable, winslash = "/", mustWork = TRUE)
  defaults$count_dir <- normalizePath(defaults$count_dir, winslash = "/", mustWork = TRUE)
  defaults$out_dir <- normalizePath(defaults$out_dir, winslash = "/", mustWork = FALSE)
  defaults$genome <- as.character(defaults$genome)
  defaults$fdr <- as.numeric(defaults$fdr)
  defaults$ntop <- as.integer(cfg$ntop)

  if (is.na(defaults$fdr) || defaults$fdr <= 0 || defaults$fdr >= 1) {
    stop("--fdr must be a number between 0 and 1.", call. = FALSE)
  }

  if (is.na(defaults$ntop) || defaults$ntop < 2) {
    stop("Internal default ntop must be an integer >= 2.", call. = FALSE)
  }

  defaults$project_name <- basename(defaults$deftable)
  defaults$project_name <- sub("\\.tsv$", "", defaults$project_name)
  defaults$project_name <- sub("^deftable_", "", defaults$project_name)

  defaults
}

read_deftable <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("deftable was not found: %s", path), call. = FALSE)
  }

  if (dir.exists(path)) {
    stop(sprintf("deftable must be a TSV file, but a directory was given: %s", path), call. = FALSE)
  }

  def <- tryCatch(
    suppressMessages(
      readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        locale = readr::locale(encoding = "UTF-8")
      )
    ),
    error = function(err) {
      stop(
        paste(
          sprintf("Failed to read deftable: %s", path),
          "Make sure --deftable points to a TSV file created by pic.",
          sprintf("Original error: %s", conditionMessage(err))
        ),
        call. = FALSE
      )
    }
  )

  required_columns <- c("count_prefix", "barcode", "sample", "group")
  missing_columns <- setdiff(required_columns, colnames(def))
  if (length(missing_columns) > 0) {
    stop(
      paste(
        sprintf("deftable is missing columns: %s", paste(missing_columns, collapse = ", ")),
        "Make sure --deftable points to a genome-specific deftable TSV created by pic."
      ),
      call. = FALSE
    )
  }

  dplyr::filter(def, grepl("[^.]", .data$group))
}

build_contrasts_from_deftable <- function(def) {
  groups <- unique(def$group)
  if (length(groups) < 2) {
    stop("At least two groups are required to build contrasts.", call. = FALSE)
  }

  contrast_pairs <- combn(groups, 2, simplify = FALSE)
  names(contrast_pairs) <- vapply(
    contrast_pairs,
    function(pair) paste0(gsub("[^[:alnum:]]+", "_", pair[[1]]), "_vs_", gsub("[^[:alnum:]]+", "_", pair[[2]])),
    character(1)
  )

  lapply(contrast_pairs, function(pair) c("group", pair[[1]], pair[[2]]))
}

build_num_umi_gene <- function(mat) {
  sample_columns <- colnames(mat)[5:ncol(mat)]
  summary <- tibble::tibble(
    sample = sample_columns,
    `#UMIs` = colSums(mat[, 5:ncol(mat), drop = FALSE]),
    `#Genes` = colSums(mat[, 5:ncol(mat), drop = FALSE] > 0)
  )

  dplyr::mutate(summary, `UMIs/Genes` = `#UMIs` / `#Genes`)
}
