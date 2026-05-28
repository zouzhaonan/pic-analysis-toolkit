#!/usr/bin/env Rscript

# 役割:
#   normalized count table と遺伝子リストから group ごとの発現図を PDF で作る。
# 入力:
#   normalized count table、deftable、遺伝子リスト、出力先。
# 出力:
#   gene ごとの normalized count 図を複数ページに並べた PDF。

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
analysis_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

source(file.path(script_dir, "analysis_plot_spec.R"))
source(file.path(script_dir, "deseq2_common.R"))

load_plot_packages <- function() {
  suppressPackageStartupMessages(library(readr))
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(tidyr))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(ggbeeswarm))
  suppressPackageStartupMessages(library(ggsci))
}

print_usage <- function() {
  help_file <- file.path(analysis_root, "lib", "help", "help_plot-expression.txt")
  cat(readChar(help_file, file.info(help_file)$size, useBytes = TRUE))
  cat("\n")
}

parse_args <- function(args) {
  parsed <- list(
    normalized_count = NULL,
    deftable = NULL,
    gene_list = NULL,
    out = NULL
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (identical(arg, "--help")) {
      print_usage()
      quit(save = "no", status = 0)
    }

    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }

    key <- sub("^--", "", arg)
    i <- i + 1
    if (i > length(args)) {
      stop(sprintf("Option requires a value: --%s", key), call. = FALSE)
    }
    value <- args[[i]]

    if (key == "count") {
      parsed$normalized_count <- value
    } else if (key == "deftable") {
      parsed$deftable <- value
    } else if (key == "gene-list") {
      parsed$gene_list <- value
    } else if (key == "out") {
      parsed$out <- value
    } else {
      stop(sprintf("Unknown option: --%s", key), call. = FALSE)
    }

    i <- i + 1
  }

  if (is.null(parsed$normalized_count) || is.null(parsed$deftable) || is.null(parsed$gene_list) || is.null(parsed$out)) {
    stop("Missing required options: --count, --deftable, --gene-list, --out", call. = FALSE)
  }

  parsed$normalized_count <- normalizePath(parsed$normalized_count, winslash = "/", mustWork = TRUE)
  parsed$deftable <- normalizePath(parsed$deftable, winslash = "/", mustWork = TRUE)
  parsed$gene_list <- normalizePath(parsed$gene_list, winslash = "/", mustWork = TRUE)
  parsed$out <- normalizePath(parsed$out, winslash = "/", mustWork = FALSE)

  parsed
}

read_gene_list <- function(path) {
  genes <- readLines(path, warn = FALSE, encoding = "UTF-8")
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  unique(genes)
}

read_normalized_count <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

extract_sample_columns <- function(normcount) {
  known_annotation <- c("ens_gene", "ext_gene", "biotype", "chr")
  colnames(normcount)[!colnames(normcount) %in% known_annotation]
}

build_sample_group_table <- function(deftable_path, sample_columns) {
  deftable <- read_deftable(deftable_path) |>
    dplyr::distinct(.data$sample, .data$group)

  sample_group <- tibble::tibble(sample = sample_columns) |>
    dplyr::left_join(deftable, by = "sample")

  if (any(is.na(sample_group$group))) {
    missing_samples <- sample_group |>
      dplyr::filter(is.na(.data$group)) |>
      dplyr::pull(.data$sample)
    stop(
      sprintf(
        "The following samples were not found in deftable: %s",
        paste(missing_samples, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  group_levels <- unique(sample_group$group)
  sample_group |>
    dplyr::mutate(group = factor(.data$group, levels = group_levels))
}

build_plot_table <- function(normcount, sample_group, target_genes) {
  sample_columns <- extract_sample_columns(normcount)

  gene_subset <- normcount |>
    dplyr::filter(.data$ens_gene %in% target_genes | (!is.na(.data$ext_gene) & .data$ext_gene %in% target_genes)) |>
    dplyr::mutate(
      gene = dplyr::if_else(!is.na(.data$ext_gene) & .data$ext_gene != "", .data$ext_gene, .data$ens_gene)
    )

  if (nrow(gene_subset) == 0) {
    stop("No requested genes were found in the normalized count table.", call. = FALSE)
  }

  plot_df <- gene_subset |>
    dplyr::select("gene", dplyr::all_of(sample_columns)) |>
    tidyr::pivot_longer(
      cols = dplyr::all_of(sample_columns),
      names_to = "sample",
      values_to = "normalized_count"
    ) |>
    dplyr::left_join(sample_group, by = "sample") |>
    dplyr::mutate(
      sample = factor(.data$sample, levels = sample_columns),
      gene = factor(.data$gene, levels = unique(.data$gene))
    )

  plot_df
}

build_gene_page_plot <- function(plot_df, panel_labels) {
  vcfg <- pic_plot_spec()$plot$violin
  palette_name <- pic_plot_spec()$plot$palette_d3

  ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$group, y = .data$normalized_count, color = .data$group, fill = .data$group)) +
    ggplot2::geom_boxplot(
      width = vcfg$box_width,
      outlier.shape = NA,
      alpha = vcfg$box_alpha,
      linewidth = vcfg$box_linewidth
    ) +
    ggbeeswarm::geom_quasirandom(
      method = "pseudorandom",
      size = vcfg$point_size,
      alpha = vcfg$point_alpha
    ) +
    ggplot2::facet_wrap(
      ~ gene_panel,
      scales = "free_y",
      nrow = vcfg$facet_nrow,
      ncol = vcfg$facet_ncol,
      drop = FALSE,
      labeller = ggplot2::as_labeller(panel_labels)
    ) +
    ggsci::scale_color_d3(palette_name, name = "Group") +
    ggsci::scale_fill_d3(palette_name, name = "Group") +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.title.x = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = vcfg$legend_title_size),
      legend.text = ggplot2::element_text(size = vcfg$legend_text_size),
      strip.background = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = vcfg$box_linewidth)
    ) +
    ggplot2::labs(y = "Normalized count")
}

save_normalized_count_pdf <- function(plot_df, output_path) {
  vcfg <- pic_plot_spec()$plot$violin
  gene_levels <- levels(plot_df$gene)
  gene_pages <- split(gene_levels, ceiling(seq_along(gene_levels) / vcfg$panels_per_page))
  panel_slots <- vcfg$panels_per_page

  grDevices::pdf(output_path, width = vcfg$pdf_width, height = vcfg$pdf_height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)

  for (gene_batch in gene_pages) {
    panel_levels <- as.character(gene_batch)
    if (length(panel_levels) < panel_slots) {
      panel_levels <- c(
        panel_levels,
        sprintf("%s%02d", vcfg$empty_panel_prefix, seq_len(panel_slots - length(panel_levels)))
      )
    }

    panel_labels <- stats::setNames(panel_levels, panel_levels)
    empty_idx <- startsWith(panel_levels, vcfg$empty_panel_prefix)
    panel_labels[empty_idx] <- ""

    page_df <- plot_df |>
      dplyr::filter(.data$gene %in% gene_batch) |>
      dplyr::mutate(gene_panel = as.character(.data$gene))

    if (any(empty_idx)) {
      dummy_df <- tidyr::expand_grid(
        gene_panel = panel_levels[empty_idx],
        group = levels(plot_df$group)
      ) |>
        dplyr::mutate(
          normalized_count = NA_real_,
          sample = NA_character_,
          gene = NA_character_
        )

      page_df <- dplyr::bind_rows(page_df, dummy_df)
    }

    page_df <- page_df |>
      dplyr::mutate(
        gene_panel = factor(.data$gene_panel, levels = panel_levels),
        group = factor(.data$group, levels = levels(plot_df$group))
      )

    print(build_gene_page_plot(page_df, panel_labels))
  }
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  load_plot_packages()

  target_genes <- read_gene_list(args$gene_list)
  normcount <- read_normalized_count(args$normalized_count)
  sample_group <- build_sample_group_table(args$deftable, extract_sample_columns(normcount))
  plot_df <- build_plot_table(normcount, sample_group, target_genes)
  save_normalized_count_pdf(plot_df, args$out)
}

tryCatch(
  main(),
  error = function(err) {
    message(sprintf("Error: %s", conditionMessage(err)))
    cat("\n")
    print_usage()
    quit(save = "no", status = 1)
  }
)
