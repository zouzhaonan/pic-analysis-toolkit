sanitize_plot_label <- function(x) {
  x <- as.character(x)
  y <- iconv(x, from = "UTF-8", to = "ASCII//TRANSLIT", sub = "")
  y[is.na(y)] <- x[is.na(y)]
  y
}

build_gsea_plot <- function(gsea_df, numerator, denominator) {
  ecfg <- pic_plot_spec()$plot$enrichment
  if (nrow(gsea_df) == 0) {
    return(NULL)
  }

  required_cols <- c("Description", "NES", "setSize", "p.adjust", "direction")
  if (!all(required_cols %in% colnames(gsea_df))) {
    return(NULL)
  }

  direction_levels <- c(
    denominator,
    numerator,
    setdiff(unique(as.character(gsea_df$direction)), c(denominator, numerator))
  )

  plot_df <- gsea_df |>
    dplyr::filter(
      !is.na(.data$Description), .data$Description != "",
      is.finite(.data$NES),
      is.finite(.data$setSize),
      is.finite(.data$p.adjust)
    ) |>
    dplyr::mutate(direction = factor(.data$direction, levels = unique(direction_levels))) |>
    dplyr::group_by(direction) |>
    dplyr::slice_max(order_by = abs(.data$NES), n = ecfg$top_n_per_direction, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(term_id = paste0("term_", dplyr::row_number()))

  if (nrow(plot_df) == 0) {
    return(NULL)
  }

  term_levels <- unlist(lapply(levels(plot_df$direction), function(dir_i) {
    x <- plot_df |>
      dplyr::filter(.data$direction == dir_i) |>
      dplyr::arrange(.data$NES)
    x$term_id
  }))

  label_map <- stats::setNames(sanitize_plot_label(plot_df$Description), plot_df$term_id)

  plot_df <- plot_df |>
    dplyr::mutate(term_label = factor(.data$term_id, levels = term_levels))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$NES, y = .data$term_label)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$setSize, color = .data$p.adjust),
      alpha = ecfg$point_alpha
    ) +
    ggplot2::facet_grid(. ~ direction, scales = "free_x", space = "fixed") +
    ggplot2::scale_y_discrete(labels = function(x) substr(unname(label_map[x]), 1, ecfg$term_label_max_chars)) +
    ggplot2::scale_color_gradientn(
      colours = c(
        ecfg$gradient_low, ecfg$gradient_high
      ),
      name = "p.adjust"
    ) +
    ggplot2::guides(
      color = ggplot2::guide_colorbar(order = 1),
      size = ggplot2::guide_legend(order = 2)
    ) +
    ggplot2::labs(x = "NES", y = NULL) +
    ggplot2::theme_bw()

  attr(p, "pic_n_y_labels") <- length(unique(as.character(plot_df$term_label)))
  p
}

build_cluster_ora_plot <- function(ora_df) {
  ecfg <- pic_plot_spec()$plot$enrichment
  if (nrow(ora_df) == 0) return(NULL)

  required_cols <- c("cluster", "Description", "GeneRatio")
  if (!all(required_cols %in% colnames(ora_df))) return(NULL)

  parse_ratio_num <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) {
      if (length(p) != 2) return(NA_real_)
      n <- suppressWarnings(as.numeric(p[[1]]))
      d <- suppressWarnings(as.numeric(p[[2]]))
      if (!is.finite(n) || !is.finite(d) || d == 0) return(NA_real_)
      n / d
    }, numeric(1))
  }
  parse_ratio_den <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) {
      if (length(p) != 2) return(NA_real_)
      suppressWarnings(as.numeric(p[[2]]))
    }, numeric(1))
  }

  plot_df <- ora_df |>
    dplyr::mutate(
      cluster = as.character(.data$cluster),
      GeneRatio_num = parse_ratio_num(.data$GeneRatio),
      GeneRatio_den = parse_ratio_den(.data$GeneRatio),
      p_adj_plot = dplyr::coalesce(
        suppressWarnings(as.numeric(.data$p.adjust)),
        suppressWarnings(as.numeric(.data$pvalue))
      )
    ) |>
    dplyr::filter(
      !is.na(.data$cluster), .data$cluster != "",
      !is.na(.data$Description), .data$Description != "",
      is.finite(.data$GeneRatio_num),
      is.finite(.data$p_adj_plot)
    ) |>
    dplyr::group_by(.data$cluster, .data$Description) |>
    dplyr::slice_min(.data$p_adj_plot, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$cluster) |>
    dplyr::slice_min(.data$p_adj_plot, n = 10, with_ties = FALSE) |>
    dplyr::ungroup()

  if (nrow(plot_df) == 0) return(NULL)

  cluster_den <- plot_df |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(
      den = {
        vals <- unique(.data$GeneRatio_den[is.finite(.data$GeneRatio_den)])
        if (length(vals) == 0) NA_real_ else vals[[1]]
      },
      .groups = "drop"
    )

  cluster_levels <- unique(plot_df$cluster)
  cluster_labels <- stats::setNames(
    vapply(cluster_levels, function(cl) {
      den <- cluster_den$den[cluster_den$cluster == cl]
      if (length(den) == 0 || !is.finite(den[[1]])) {
        sprintf("%s (N=NA)", cl)
      } else {
        sprintf("%s (N=%s)", cl, format(den[[1]], scientific = FALSE, trim = TRUE))
      }
    }, character(1)),
    cluster_levels
  )

  plot_df <- plot_df |>
    dplyr::mutate(
      cluster_label = factor(cluster_labels[.data$cluster], levels = cluster_labels[cluster_levels])
    )

  plot_df <- plot_df |>
    dplyr::mutate(Description_plot = sanitize_plot_label(.data$Description))

  term_levels <- plot_df |>
    dplyr::group_by(.data$Description_plot) |>
    dplyr::summarise(best_padj = min(.data$p_adj_plot, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$best_padj, .data$Description_plot) |>
    dplyr::pull(.data$Description_plot)

  plot_df <- plot_df |>
    dplyr::mutate(term_label = factor(.data$Description_plot, levels = rev(term_levels)))

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$cluster_label, y = .data$term_label)
  ) +
    ggplot2::geom_point(
      ggplot2::aes(size = .data$GeneRatio_num, color = .data$p_adj_plot),
      alpha = ecfg$point_alpha
    ) +
    ggplot2::scale_y_discrete(labels = function(x) substr(x, 1, ecfg$term_label_max_chars)) +
    ggplot2::scale_color_gradientn(
      colours = c(ecfg$gradient_low, ecfg$gradient_high),
      name = "p.adjust"
    ) +
    ggplot2::labs(x = "Cluster", y = NULL, size = "GeneRatio") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  attr(p, "pic_n_y_labels") <- length(unique(as.character(plot_df$term_label)))
  p
}

save_plot <- function(plot_obj, out_dir, out_name) {
  ecfg <- pic_plot_spec()$plot$enrichment
  out_file <- file.path(out_dir, out_name)
  n_labels <- suppressWarnings(as.numeric(attr(plot_obj, "pic_n_y_labels")))
  if (!is.finite(n_labels)) n_labels <- NA_real_

  height_auto <- if (is.finite(n_labels)) {
    max(ecfg$save_height_min, min(ecfg$save_height_max, 2 + n_labels * ecfg$save_height_per_label))
  } else {
    ecfg$save_height
  }

  ggplot2::ggsave(out_file, plot = plot_obj, width = ecfg$save_width, height = height_auto)
}
