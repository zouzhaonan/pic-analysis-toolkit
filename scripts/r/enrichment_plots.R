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

# ---------------------------------------------------------------------------
# plotly 用 dot plot spec (report のインタラクティブ描画)
# ---------------------------------------------------------------------------

# 連続値を plotly マーカーの px サイズに線形写像する。
plotly_marker_size <- function(v, min_px = 5, max_px = 17) {
  v <- suppressWarnings(as.numeric(v))
  fin <- v[is.finite(v)]
  if (length(fin) == 0) return(rep((min_px + max_px) / 2, length(v)))
  lo <- min(fin); hi <- max(fin)
  if (!is.finite(lo) || !is.finite(hi) || hi == lo) return(rep((min_px + max_px) / 2, length(v)))
  out <- min_px + (v - lo) / (hi - lo) * (max_px - min_px)
  out[!is.finite(out)] <- min_px
  out
}

# p.adjust を色 (小さい=有意=赤 -> 大きい=青) に写像する colorscale。
enrich_colorscale <- function() {
  ecfg <- pic_plot_spec()$plot$enrichment
  list(list(0, ecfg$gradient_low), list(1, ecfg$gradient_high))
}

# GSEA を plotly の dot plot spec (list(data, layout, n_terms)) に変換する。
# 元の facet_grid(. ~ direction, scales="free_x") と同様に、direction ごとに
# 独立した x 軸 (free_x) のパネルを左右に並べ、y (term ラベル) は共有して左端に揃える。
# color=p.adjust, size=setSize。
build_gsea_plotly <- function(gsea_df, numerator, denominator) {
  ecfg <- pic_plot_spec()$plot$enrichment
  required_cols <- c("Description", "NES", "setSize", "p.adjust", "direction")
  if (nrow(gsea_df) == 0 || !all(required_cols %in% colnames(gsea_df))) return(NULL)

  direction_levels <- c(denominator, numerator,
    setdiff(unique(as.character(gsea_df$direction)), c(denominator, numerator)))

  plot_df <- gsea_df |>
    dplyr::filter(!is.na(.data$Description), .data$Description != "",
                  is.finite(.data$NES), is.finite(.data$setSize), is.finite(.data$p.adjust)) |>
    dplyr::mutate(direction = factor(.data$direction, levels = unique(direction_levels))) |>
    dplyr::group_by(.data$direction) |>
    dplyr::slice_max(order_by = abs(.data$NES), n = ecfg$top_n_per_direction, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::mutate(term_id = paste0("term_", dplyr::row_number()))
  if (nrow(plot_df) == 0) return(NULL)

  # y 順 (下→上): direction level ごとに NES 昇順
  term_levels <- unlist(lapply(levels(plot_df$direction), function(dir_i) {
    x <- plot_df |> dplyr::filter(.data$direction == dir_i) |> dplyr::arrange(.data$NES)
    x$term_id
  }))
  plot_df <- plot_df[match(term_levels, plot_df$term_id), , drop = FALSE]

  desc_full <- as.character(plot_df$Description)
  ylab <- substr(sanitize_plot_label(plot_df$Description), 1, ecfg$term_label_max_chars)
  ypos_all <- seq_len(nrow(plot_df))
  padj_all <- as.numeric(plot_df$p.adjust)
  size_all <- plotly_marker_size(as.numeric(plot_df$setSize))
  cmin <- min(padj_all); cmax <- max(padj_all)

  ht <- paste0("<b>%{customdata[0]}</b><br>enriched in: %{customdata[3]}<br>",
               "NES: %{x:.2f}<br>setSize: %{customdata[1]}<br>p.adjust: %{customdata[2]:.3g}<extra></extra>")

  dirs <- levels(plot_df$direction)
  dirs <- dirs[dirs %in% unique(as.character(plot_df$direction))]
  ndir <- length(dirs)
  gap <- 0.035
  panelw <- (1 - gap * (ndir - 1)) / ndir

  traces <- list()
  annotations <- list()
  layout <- list(
    margin = list(l = 10, r = 10, t = 30, b = 46), hovermode = "closest",
    yaxis = list(tickmode = "array", tickvals = as.list(ypos_all), ticktext = as.list(ylab),
                 automargin = TRUE, range = list(0.5, nrow(plot_df) + 0.5),
                 domain = list(0, 1), anchor = "x")
  )
  for (i in seq_len(ndir)) {
    d <- dirs[[i]]
    idx <- which(as.character(plot_df$direction) == d)
    x0 <- (i - 1) * (panelw + gap); x1 <- x0 + panelw
    axref <- if (i == 1) "x" else paste0("x", i)
    axname <- if (i == 1) "xaxis" else paste0("xaxis", i)
    layout[[axname]] <- list(title = list(text = "NES"), domain = list(x0, x1),
                             anchor = "y", zeroline = FALSE, automargin = TRUE)
    showscale <- (i == ndir)
    mk <- list(size = as.list(size_all[idx]), color = as.list(padj_all[idx]),
               colorscale = enrich_colorscale(), cmin = cmin, cmax = cmax,
               showscale = showscale, line = list(width = 0.5, color = "rgba(0,0,0,.35)"),
               opacity = ecfg$point_alpha)
    if (showscale) mk$colorbar <- list(title = list(text = "p.adjust", side = "right"),
                                        thickness = 12, len = 0.6)
    traces[[length(traces) + 1L]] <- list(
      x = as.list(as.numeric(plot_df$NES[idx])), y = as.list(ypos_all[idx]),
      xaxis = axref, yaxis = "y",
      customdata = lapply(idx, function(j) list(desc_full[j], as.numeric(plot_df$setSize[j]), padj_all[j], d)),
      mode = "markers", type = "scatter", marker = mk, hovertemplate = ht
    )
    # facet strip: パネル上部に direction (enriched group) 名
    annotations[[length(annotations) + 1L]] <- list(
      text = html_escape(d), xref = "paper", yref = "paper",
      x = (x0 + x1) / 2, y = 1.012, xanchor = "center", yanchor = "bottom",
      showarrow = FALSE, font = list(size = 12, color = "#1f2933"),
      bgcolor = "#eef3f9", bordercolor = "#cdd7e2", borderwidth = 1, borderpad = 3
    )
  }
  layout$annotations <- annotations
  list(data = traces, layout = layout, n_terms = nrow(plot_df))
}

# ORA (cluster) を plotly の dot plot spec に変換する。
# x=cluster, y=term, color=p.adjust, size=GeneRatio。
build_ora_plotly <- function(ora_df) {
  ecfg <- pic_plot_spec()$plot$enrichment
  required_cols <- c("cluster", "Description", "GeneRatio")
  if (nrow(ora_df) == 0 || !all(required_cols %in% colnames(ora_df))) return(NULL)

  parse_ratio_num <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) {
      if (length(p) != 2) return(NA_real_)
      n <- suppressWarnings(as.numeric(p[[1]])); d <- suppressWarnings(as.numeric(p[[2]]))
      if (!is.finite(n) || !is.finite(d) || d == 0) return(NA_real_)
      n / d
    }, numeric(1))
  }
  parse_ratio_den <- function(x) {
    parts <- strsplit(as.character(x), "/", fixed = TRUE)
    vapply(parts, function(p) if (length(p) != 2) NA_real_ else suppressWarnings(as.numeric(p[[2]])), numeric(1))
  }

  plot_df <- ora_df |>
    dplyr::mutate(
      cluster = as.character(.data$cluster),
      GeneRatio_num = parse_ratio_num(.data$GeneRatio),
      GeneRatio_den = parse_ratio_den(.data$GeneRatio),
      p_adj_plot = dplyr::coalesce(suppressWarnings(as.numeric(.data$p.adjust)),
                                   suppressWarnings(as.numeric(.data$pvalue)))
    ) |>
    dplyr::filter(!is.na(.data$cluster), .data$cluster != "",
                  !is.na(.data$Description), .data$Description != "",
                  is.finite(.data$GeneRatio_num), is.finite(.data$p_adj_plot)) |>
    dplyr::group_by(.data$cluster, .data$Description) |>
    dplyr::slice_min(.data$p_adj_plot, n = 1, with_ties = FALSE) |>
    dplyr::ungroup() |>
    dplyr::group_by(.data$cluster) |>
    dplyr::slice_min(.data$p_adj_plot, n = 10, with_ties = FALSE) |>
    dplyr::ungroup()
  if (nrow(plot_df) == 0) return(NULL)

  cluster_den <- plot_df |>
    dplyr::group_by(.data$cluster) |>
    dplyr::summarise(den = {
      vals <- unique(.data$GeneRatio_den[is.finite(.data$GeneRatio_den)])
      if (length(vals) == 0) NA_real_ else vals[[1]]
    }, .groups = "drop")

  cluster_levels <- unique(plot_df$cluster)
  cluster_labels <- stats::setNames(
    vapply(cluster_levels, function(cl) {
      den <- cluster_den$den[cluster_den$cluster == cl]
      if (length(den) == 0 || !is.finite(den[[1]])) sprintf("%s (N=NA)", cl)
      else sprintf("%s (N=%s)", cl, format(den[[1]], scientific = FALSE, trim = TRUE))
    }, character(1)), cluster_levels)

  plot_df <- plot_df |>
    dplyr::mutate(cluster_label = factor(cluster_labels[.data$cluster], levels = cluster_labels[cluster_levels]),
                  Description_plot = sanitize_plot_label(.data$Description))

  term_levels <- plot_df |>
    dplyr::group_by(.data$Description_plot) |>
    dplyr::summarise(best_padj = min(.data$p_adj_plot, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$best_padj, .data$Description_plot) |>
    dplyr::pull(.data$Description_plot)

  plot_df <- plot_df |>
    dplyr::mutate(term_label = factor(.data$Description_plot, levels = rev(term_levels)))

  ylab_levels <- substr(levels(plot_df$term_label), 1, ecfg$term_label_max_chars)
  yy <- as.integer(plot_df$term_label)
  xx <- as.character(plot_df$cluster_label)
  gr <- as.numeric(plot_df$GeneRatio_num)
  padj <- as.numeric(plot_df$p_adj_plot)

  cd <- lapply(seq_len(nrow(plot_df)), function(i)
    list(as.character(plot_df$Description_plot[i]), gr[i], padj[i], xx[i]))
  ht <- paste0("<b>%{customdata[0]}</b><br>cluster: %{customdata[3]}<br>",
               "GeneRatio: %{customdata[1]:.3g}<br>p.adjust: %{customdata[2]:.3g}<extra></extra>")

  trace <- list(
    x = as.list(xx), y = as.list(yy), customdata = cd,
    mode = "markers", type = "scatter",
    marker = list(
      size = as.list(plotly_marker_size(gr)),
      color = as.list(padj), colorscale = enrich_colorscale(),
      colorbar = list(title = list(text = "p.adjust", side = "right"), thickness = 12, len = 0.6),
      showscale = TRUE, line = list(width = 0.5, color = "rgba(0,0,0,.35)"),
      opacity = ecfg$point_alpha
    ),
    hovertemplate = ht
  )
  layout <- list(
    xaxis = list(type = "category", categoryorder = "array",
                 categoryarray = as.list(levels(plot_df$cluster_label)),
                 tickangle = -40, automargin = TRUE),
    yaxis = list(tickmode = "array", tickvals = as.list(seq_along(ylab_levels)),
                 ticktext = as.list(ylab_levels), automargin = TRUE,
                 range = list(0.5, length(ylab_levels) + 0.5)),
    margin = list(l = 10, r = 10, t = 16, b = 80), hovermode = "closest"
  )
  list(data = list(trace), layout = layout, n_terms = length(ylab_levels))
}

# dot plot の term 数からプロット高さ (px) を決める。
enrich_plot_height <- function(n_terms) {
  as.integer(max(300, min(1500, 110 + n_terms * 20)))
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
