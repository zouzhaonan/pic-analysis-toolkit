# 役割: DEG (MA/Volcano) + Aggregation + Gene Expression + Cluster profile セクション。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。

# ---------------------------------------------------------------------------
# DESeq2: contrast ごとの MA / volcano (plotly)
# ---------------------------------------------------------------------------

# 背景 (非有意) 点を上限まで間引く。有意点は全て残す。
downsample_idx <- function(is_sig, cap = 3500L) {
  sig_idx <- which(is_sig)
  bg_idx <- which(!is_sig)
  if (length(bg_idx) > cap) {
    set.seed(42L)
    bg_idx <- sort(sample(bg_idx, cap))
  }
  sort(c(sig_idx, bg_idx))
}

scatter_traces_by_dir <- function(df, x, y, numerator, denominator, hovertemplate, group_pal = NULL) {
  # df は idx で絞り込み済み。direction 列 (up/down/ns) でトレース分割。
  # df には cd1 (padj), cd2 (pvalue) 列を含み、customdata として hover に渡す。
  # 凡例の group 名は group 色付き文字。
  # 凡例の group 名はドット色 (up=赤 / down=青) に合わせる
  cats <- list(
    up = list(name = sprintf('<span style="color:#d7301f">%s ↑</span>', html_escape(numerator)), color = "#d7301f"),
    down = list(name = sprintf('<span style="color:#2166ac">%s ↑</span>', html_escape(denominator)), color = "#2166ac"),
    ns = list(name = "not significant", color = "#bdbdbd")
  )
  traces <- list()
  for (k in c("ns", "down", "up")) {
    sub <- df[df$dir == k, , drop = FALSE]
    if (nrow(sub) == 0) next
    cd <- lapply(seq_len(nrow(sub)), function(i) list(sub$cd1[[i]], sub$cd2[[i]]))
    traces[[length(traces) + 1]] <- list(
      x = as.list(as.numeric(sub[[x]])),
      y = as.list(as.numeric(sub[[y]])),
      text = as.list(as.character(sub$gene)),
      customdata = cd,
      name = cats[[k]]$name,
      mode = "markers",
      type = "scatter",
      showlegend = (k != "ns"),  # not significant は凡例に出さない (点は描画)
      marker = list(size = if (k == "ns") 4 else 6, color = cats[[k]]$color,
                    opacity = if (k == "ns") 0.45 else 0.85),
      hovertemplate = hovertemplate
    )
  }
  traces
}

build_contrast_plots <- function(reg, stats, deg_counts, fdr, id_prefix = "", stats_src = "", group_order = NULL, group_pal = NULL) {
  aspects <- unique(stats$aspect)
  # contrast の並び順: sample_sheet の group 順 (なければ DEG 数で降順)
  if (!is.null(group_order)) {
    aspects <- pic_order_contrasts(aspects, group_order)
  } else {
    deg_per <- vapply(aspects, function(a) {
      if (!is.null(deg_counts) && a %in% names(deg_counts)) deg_counts[[a]] else {
        sum(stats$aspect == a & !is.na(stats$padj) & stats$padj < fdr)
      }
    }, numeric(1))
    aspects <- aspects[order(deg_per, decreasing = TRUE)]
  }
  deg_per <- vapply(aspects, function(a) {
    if (!is.null(deg_counts) && a %in% names(deg_counts)) deg_counts[[a]] else {
      sum(stats$aspect == a & !is.na(stats$padj) & stats$padj < fdr)
    }
  }, numeric(1))

  items <- list()
  for (k in seq_along(aspects)) {
    a <- aspects[[k]]
    nd <- deg_per[[k]]
    cs <- stats[stats$aspect == a, , drop = FALSE]
    sp <- strsplit(a, " / ", fixed = TRUE)[[1]]
    numerator <- if (length(sp) >= 1) sp[[1]] else "num"
    denominator <- if (length(sp) >= 2) sp[[2]] else "den"
    # 表示用に deftable どおりの casing へ (stats$aspect は小文字化されている)
    numerator <- group_name_of(numerator, group_pal); denominator <- group_name_of(denominator, group_pal)
    flab <- format_contrast_file_label(a)

    gene <- if ("ext_gene" %in% colnames(cs)) {
      g <- as.character(cs$ext_gene); g[is.na(g) | g == ""] <- as.character(cs$ens_gene)[is.na(g) | g == ""]; g
    } else as.character(cs$ens_gene)

    lfc <- suppressWarnings(as.numeric(cs$log2FoldChange))
    padj <- suppressWarnings(as.numeric(cs$padj))
    pval <- suppressWarnings(as.numeric(cs$pvalue))
    base <- suppressWarnings(as.numeric(cs$baseMean))
    is_sig <- !is.na(padj) & padj < fdr & !is.na(lfc)
    dir <- ifelse(is_sig & lfc > 0, "up", ifelse(is_sig & lfc < 0, "down", "ns"))

    legend_top <- list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0)
    # 縦軸 (fold change) ラベル: log2(numerator / denominator) を group 色で
    grp_col <- function(g) group_color_or(g, group_pal)
    lfc_lab <- sprintf('log<sub>2</sub>(%s / %s)', html_escape(numerator), html_escape(denominator))
    cells <- character(0)

    # ---- M-A (hover に padj) ----
    keepMA <- is.finite(base) & base > 0 & is.finite(lfc)
    idxMA <- downsample_idx(is_sig & keepMA)
    idxMA <- idxMA[keepMA[idxMA]]
    if (length(idxMA) > 0) {
      dfMA <- data.frame(x = base[idxMA], y = lfc[idxMA], gene = gene[idxMA], dir = dir[idxMA],
                         cd1 = padj[idxMA], cd2 = pval[idxMA], stringsAsFactors = FALSE)
      ht <- "<b>%{text}</b><br>baseMean: %{x:.3g}<br>log<sub>2</sub>FC: %{y:.2f}<br>p.adjust: %{customdata[0]:.3g}<extra></extra>"
      tr <- scatter_traces_by_dir(dfMA, "x", "y", numerator, denominator, ht, group_pal)
      id <- sprintf("%sma_%s", id_prefix, flab)
      layout <- list(
        xaxis = list(title = "mean expression (baseMean)", type = "log"),
        yaxis = list(title = lfc_lab, zeroline = TRUE),
        hovermode = "closest", margin = list(l = 64, r = 20, t = 28, b = 46),
        legend = legend_top
      )
      register_plot(reg, id, list(data = tr, layout = layout))
      cells <- c(cells, sprintf('<div class="pic-plot-cell"><h4>M-A</h4><div id="%s" class="pic-plot"></div></div>', id))
    }

    # ---- Volcano: y=-log10(p-value) 一本化、hover に padj も表示 ----
    keepV <- is.finite(pval) & pval > 0 & is.finite(lfc)
    idxV <- downsample_idx(is_sig & keepV)
    idxV <- idxV[keepV[idxV]]
    if (length(idxV) > 0) {
      dfV <- data.frame(x = lfc[idxV], y = -log10(pval[idxV]), gene = gene[idxV], dir = dir[idxV],
                        cd1 = padj[idxV], cd2 = pval[idxV], stringsAsFactors = FALSE)
      ht <- "<b>%{text}</b><br>log<sub>2</sub>FC: %{x:.2f}<br>p-value: %{customdata[1]:.3g}<br>p.adjust: %{customdata[0]:.3g}<extra></extra>"
      tr <- scatter_traces_by_dir(dfV, "x", "y", numerator, denominator, ht, group_pal)
      id <- sprintf("%svolcano_%s", id_prefix, flab)
      layout <- list(
        xaxis = list(title = lfc_lab, zeroline = TRUE),
        yaxis = list(title = "−log<sub>10</sub>(p-value)"),
        hovermode = "closest", margin = list(l = 60, r = 20, t = 28, b = 46),
        legend = legend_top
      )
      register_plot(reg, id, list(data = tr, layout = layout))
      cells <- c(cells, sprintf('<div class="pic-plot-cell"><h4>Volcano</h4><div id="%s" class="pic-plot"></div></div>', id))
    }

    item_html <- paste0('<div class="pic-ma-vol">', paste(cells, collapse = ""), '</div>')
    items[[length(items) + 1L]] <- list(
      id = sprintf("%ssel_%s", id_prefix, flab),
      label = a, aspect = a, count = nd,
      html = item_html,
      checked = (k == 1L)
    )
  }

  # 説明文 (info バッジに集約される)。<div class="pic-degsel"> 内 (M-A/Volcano チャンク) に置く。
  note <- paste0(
    '<p class="pic-note">Each gene is one point. The <b>M-A plot</b> shows the log<sub>2</sub> fold change vs. mean expression; ',
    'the <b>Volcano plot</b> shows it vs. statistical significance (&minus;log<sub>10</sub> p-value). ',
    '<b style="color:#d7301f">Red</b> / <b style="color:#2166ac">blue</b> points are significantly up / down (p.adjust below the FDR); grey is not significant. ',
    'Hover a point for its gene name, log<sub>2</sub> fold change and p.adjust.</p>')
  # 右パネル先頭: 説明文 (info へ集約) + タイトル行 (Download csv を隣に)
  view_head <- paste0(note, sub_head('M-A &amp; Volcano Plot', src_note(stats_src)))

  # contrast 選択: group×group 行列 (各セルに radio)。1 つだけなら行列不要。
  groups <- if (!is.null(group_order) && length(group_order) > 0) group_order else
    unique(unlist(lapply(aspects, function(a) trimws(strsplit(a, " / ", fixed = TRUE)[[1]]))))
  entries <- lapply(items, function(it) list(
    aspect = it$aspect, count = it$count, checked = it$checked,
    attr = sprintf('data-target="%s"', it$id)))
  matrix_html <- build_group_matrix(entries, groups, group_pal, show_count = FALSE,
                                    radio_name = paste0(id_prefix, "degcon"))
  blocks <- vapply(items, function(it) {
    hid <- if (isTRUE(it$checked)) "" else " hidden"
    sprintf('<div class="pic-select-item" id="%s"%s>%s</div>', it$id, hid, it$html)
  }, character(1))
  sel_html <- if (length(items) <= 1) {
    view_body <- if (length(items) == 1) sprintf('<div class="pic-select-item">%s</div>', items[[1]]$html) else ""
    sprintf('<div class="pic-degsel"><div class="pic-degsel-ctrl"><!--SUBNAV--></div><div class="pic-degsel-view">%s%s</div></div>', view_head, view_body)
  } else if (nzchar(matrix_html)) {
    # 左 = サブタブ + contrast 行列 (radio)、右 = タイトル行 + MA/volcano ブロック
    paste0('<div class="pic-degsel"><div class="pic-degsel-ctrl"><!--SUBNAV--><h4>Comparison</h4>',
           matrix_html, '</div><div class="pic-degsel-view">', view_head,
           paste(blocks, collapse = ""), '</div></div>')
  } else {
    sprintf('<div class="pic-degsel"><div class="pic-degsel-ctrl"><!--SUBNAV--></div><div class="pic-degsel-view">%s%s</div></div>',
            view_head, build_select_group(items))
  }
  sel_html
}

# ---------------------------------------------------------------------------
# enrich: GSEA (contrast ごと) / ORA (method ごと) を再生成して内包
# ---------------------------------------------------------------------------

png_data_uri <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  paste0("data:image/png;base64,", base64enc::base64encode(raw))
}

# <base_dir>/aggregate_profile_<run>_<genome>.csv からインタラクティブな
# メタジーン折れ線 (plotly) を作る。
# (平坦化: aggregate/ フォルダも PNG も廃し、プロジェクト直下の CSV を suffix で識別)
build_aggregate_html <- function(reg, base_dir, id_prefix = "", proj_dir = NULL, sample_order = NULL, group_pal = NULL) {
  csvs <- list.files(base_dir, pattern = "^aggregate_profile_.*\\.csv$", full.names = TRUE)
  blocks <- character(0)
  for (csv in sort(csvs)) {
    proj <- sub("\\.csv$", "", sub("^aggregate_profile_", "", basename(csv)))
    {
      d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE, progress = FALSE))
      d <- as.data.frame(d, check.names = FALSE)
      if (nrow(d) == 0 || !all(c("sample", "group", "pos", "value") %in% colnames(d))) next
      pal <- ensure_palette(group_pal, unique(as.character(d$group)))
      traces <- list()
      # 凡例は group 単位 (各 group の最初の trace のみ凡例に表示、他は legendgroup で連動)。
      # サンプルは sample_sheet の順に描画し、hover にサンプル名を表示。
      seen <- character(0)
      for (s in pic_reorder_vec(unique(as.character(d$sample)), sample_order)) {
        sub <- d[as.character(d$sample) == s, , drop = FALSE]
        sub <- sub[order(sub$pos), , drop = FALSE]
        g <- as.character(sub$group[[1]])
        first <- !(g %in% seen); seen <- c(seen, g)
        traces[[length(traces) + 1]] <- list(
          x = as.list(as.numeric(sub$pos)), y = as.list(as.numeric(sub$value)),
          name = if (first) group_span(g, pal) else html_escape(s), legendgroup = g, showlegend = first,
          mode = "lines", type = "scatter",
          line = list(color = unname(pal[[g]]), width = 1.4),
          hovertemplate = sprintf("%s<br>CPM: %%{y:.2f}<extra></extra>", html_escape(s)))
      }
      id <- sprintf("%saggregate_%s", id_prefix, proj)
      # フランク端ラベル (±Nkb / ±Nbp)
      fb <- if ("flank_bp" %in% colnames(d)) suppressWarnings(as.numeric(d$flank_bp[[1]])) else NA
      flab <- if (is.finite(fb)) { if (fb >= 1000) sprintf("%gkb", fb / 1000) else sprintf("%dbp", as.integer(fb)) } else "flank"
      xmin <- min(as.numeric(d$pos)); xmax <- max(as.numeric(d$pos))
      tickvals <- list(xmin, 0, 100, xmax)
      ticktext <- list(paste0("-", flab), "TSS", "TES", paste0("+", flab))
      vline <- function(x) list(type = "line", x0 = x, x1 = x, yref = "paper", y0 = 0, y1 = 1,
                                line = list(color = "#9aa4af", width = 1, dash = "dash"))
      layout <- list(
        xaxis = list(title = "gene body (TSS → TES)", tickmode = "array", tickvals = tickvals, ticktext = ticktext, zeroline = FALSE),
        yaxis = list(title = "mean CPM", autorange = TRUE),
        shapes = list(vline(0), vline(100)),
        hovermode = "closest", margin = list(l = 60, r = 20, t = 30, b = 46),
        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0))
      register_plot(reg, id, list(data = traces, layout = layout))
      blocks <- c(blocks, sprintf(
        '<div class="pic-plot-cell">%s<div id="%s" class="pic-plot"></div></div>',
        src_note(report_rel_path(csv, proj_dir)), id))
    }
  }
  if (length(blocks) == 0) return("")
  paste0('<p class="pic-note">Average read coverage along the gene body, scaled from the transcription start (TSS) to the end (TES) ',
         'with short flanks on each side. Each line is one sample, colored by group &mdash; click a group in the legend to show or hide it. ',
         'This shows where reads accumulate along genes (e.g. a 3&prime; bias).</p><div class="pic-plot-grid">',
         paste(blocks, collapse = ""), '</div>')
}

# 正規化カウント (normalizedCountTable) を id 付きで登録する (HTML に埋め込む用)。
# stats があれば pvalue 最小の 5 遺伝子を default として保持する。
register_expr_data <- function(reg, id, deseq2_dir, project, group_map, group_pal, stats = NULL, sample_order = NULL, group_order = NULL) {
  f <- file.path(deseq2_dir, sprintf("normalizedCountTable_%s.csv", project))
  if (!file.exists(f)) return(FALSE)
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  d <- as.data.frame(d, check.names = FALSE)
  ann <- c("ens_gene", "ext_gene", "biotype", "chr")
  sample_cols <- setdiff(colnames(d), ann)
  if (length(sample_cols) < 1 || nrow(d) < 1) return(FALSE)
  # サンプルを group 順 -> sample_sheet 順に並べる
  grp <- vapply(sample_cols, function(s) if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else "all", character(1))
  grank <- if (!is.null(group_order)) { r <- match(grp, group_order); r[is.na(r)] <- length(group_order) + 1L; r } else match(grp, sort(unique(grp)))
  srank <- if (!is.null(sample_order)) { r <- match(sample_cols, sample_order); r[is.na(r)] <- length(sample_order) + 1L; r } else seq_along(sample_cols)
  ord <- order(grank, srank, sample_cols)
  sample_cols <- sample_cols[ord]; grp <- unname(grp[ord])
  ug <- unique(grp)
  pal <- ensure_palette(group_pal, ug)
  mat <- as.matrix(d[, sample_cols, drop = FALSE]); storage.mode(mat) <- "double"
  mat <- round(mat, 2)
  ens_all <- as.character(d$ens_gene)
  ext <- if ("ext_gene" %in% colnames(d)) as.character(d$ext_gene) else rep(NA_character_, nrow(d))

  # default = pvalue 最小の 5 遺伝子 (contrast 横断で遺伝子ごとに min pvalue)、
  # および contrast ごとの padj (q-value, hover 表示用)。
  default_ens <- character(0)
  contrasts <- character(0)
  qval <- NULL
  if (!is.null(stats) && all(c("ens_gene", "aspect") %in% colnames(stats))) {
    stats$ens_gene <- as.character(stats$ens_gene)
    if ("pvalue" %in% colnames(stats)) {
      st <- stats[!is.na(stats$pvalue), c("ens_gene", "pvalue")]
      if (nrow(st) > 0) {
        agg <- stats::aggregate(pvalue ~ ens_gene, data = st, FUN = min)
        cand <- as.character(agg$ens_gene[order(agg$pvalue)])
        default_ens <- head(cand[cand %in% ens_all], 5)
      }
    }
    if ("padj" %in% colnames(stats)) {
      contrasts <- unique(as.character(stats$aspect))
      if (!is.null(group_order)) contrasts <- pic_order_contrasts(contrasts, group_order)
      qmat <- matrix(NA_real_, nrow = length(ens_all), ncol = length(contrasts))
      for (ci in seq_along(contrasts)) {
        sub <- stats[as.character(stats$aspect) == contrasts[ci], c("ens_gene", "padj")]
        idx <- match(ens_all, sub$ens_gene)
        qmat[, ci] <- signif(as.numeric(sub$padj)[idx], 3)
      }
      qval <- lapply(seq_len(nrow(qmat)), function(i) as.numeric(qmat[i, ]))
      # 表示用 contrast 名は deftable どおりの casing に (stats$aspect は小文字化されている)
      # 注: 名前付きベクトルにすると JSON が配列でなくオブジェクトになり JS 側で map できないため unname
      contrasts <- unname(vapply(contrasts, function(cc) paste(vapply(strsplit(cc, " / ", fixed = TRUE)[[1]], function(g) group_name_of(g, group_pal), character(1)), collapse = " / "), character(1)))
    }
  }

  reg$expr[[id]] <- list(
    samples = as.list(sample_cols),
    groups = as.list(grp),
    palette = as.list(pal[ug]),
    ens = as.list(ens_all),
    ext = as.list(ext),
    counts = lapply(seq_len(nrow(mat)), function(i) as.numeric(mat[i, ])),
    default = as.list(default_ens),
    contrasts = as.list(contrasts),
    qval = qval,
    fdr = pic_plot_spec()$defaults$fdr
  )
  TRUE
}

# 遺伝子発現セクション: 入力した Gene ID/symbol を box+beeswarm で表示する UI。
build_expression_section <- function(reg, deseq2_dir, project, group_map, group_pal, stats, id, section_id, heading, proj_dir = NULL, sample_order = NULL, group_order = NULL) {
  if (!register_expr_data(reg, id, deseq2_dir, project, group_map, group_pal, stats, sample_order, group_order)) return("")
  countf <- file.path(deseq2_dir, sprintf("normalizedCountTable_%s.csv", project))
  fdr <- pic_plot_spec()$defaults$fdr
  sprintf(paste0(
    '<section id="%s"><h2>%s</h2>',
    '<p class="pic-note">Compare the expression of individual genes across sample groups. ',
    'Start typing a gene ID or symbol and pick it from the list to add it as a tag; ',
    'click <b>&times;</b> on a tag (or press Backspace) to remove it. The plot updates automatically. ',
    'Each box shows the distribution of normalized counts for one group; dots are individual samples.</p>',
    '<p class="pic-note">A gene name turns <b style="color:#d7301f">red</b> when its adjusted p-value (padj) ',
    'is below the significance threshold (FDR = %s) in at least one comparison. ',
    'Hover over a gene name to see its padj per comparison; comparisons with padj&nbsp;=&nbsp;1 are omitted, ',
    'and if none remain the gene is flagged as not significant in any comparison.</p>',
    '%s',
    '<div class="pic-expr" data-expr="%s">',
    '<div class="pic-expr-side"><h4>Genes</h4>',
    '<div class="pic-expr-typeahead"><div class="pic-expr-tagbox">',
    '<input class="pic-expr-input" type="text" autocomplete="off" placeholder="Add a gene…"></div>',
    '<div class="pic-expr-ac"></div></div>',
    '<div class="pic-expr-msg"></div></div>',
    '<div class="pic-expr-plots"></div></div></section>'),
    section_id, html_escape(heading), format(fdr, trim = TRUE),
    src_note(report_rel_path(countf, proj_dir)), id)
}

# DEG クラスタの挙動 (group ごとの rlog 分布) を facet で示す図。
# 各 cluster_N が何を意味するか (どの group で高い/低い) を ORA の前に提示する。
build_cluster_profile_png <- function(deseq2_dir, project, tmp_dir, group_pal = NULL) {
  f <- file.path(degcluster_dir_of(deseq2_dir, project), sprintf("DEGCluster_profile_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  prof <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  if (nrow(prof) == 0 || !all(c("rlog_expr", "cluster_id", "group") %in% colnames(prof))) return(NULL)

  # cluster ラベルに遺伝子数を付与
  summ_f <- file.path(degcluster_dir_of(deseq2_dir, project), sprintf("DEGCluster_summary_%s.csv", project))
  if (file.exists(summ_f)) {
    summ <- suppressMessages(readr::read_csv(summ_f, show_col_types = FALSE, progress = FALSE))
    nmap <- stats::setNames(summ$gene_count, summ$cluster_id)
    prof$cluster_lab <- ifelse(
      prof$cluster_id %in% names(nmap),
      sprintf("%s (n=%s)", prof$cluster_id, nmap[prof$cluster_id]),
      as.character(prof$cluster_id)
    )
  } else {
    prof$cluster_lab <- as.character(prof$cluster_id)
  }
  prof$group <- factor(prof$group, levels = unique(as.character(prof$group)))

  p <- ggplot2::ggplot(prof, ggplot2::aes(x = .data$group, y = .data$rlog_expr, color = .data$group, fill = .data$group)) +
    ggplot2::geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.5, linewidth = 0.3) +
    ggplot2::facet_wrap(~cluster_lab, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      legend.position = "none",
      strip.background = ggplot2::element_blank()
    ) +
    ggplot2::labs(x = NULL, y = "rlog expression")
  p <- tryCatch(
    p + ggbeeswarm::geom_quasirandom(method = "pseudorandom", size = 0.25, alpha = 0.6),
    error = function(e) p
  )
  if (!is.null(group_pal)) {
    gp <- group_pal[names(group_pal) %in% levels(prof$group)]
    p <- p + ggplot2::scale_color_manual(values = gp) + ggplot2::scale_fill_manual(values = gp)
  }

  ncl <- length(unique(prof$cluster_id))
  ncol <- min(3, max(1, ncl))
  nrow_f <- ceiling(ncl / ncol)
  h <- max(2.6, min(12, 1.2 + nrow_f * 2.2))
  png_path <- file.path(tmp_dir, "cluster_profile.png")
  ggplot2::ggsave(png_path, plot = p, width = 10, height = h, dpi = 120, limitsize = FALSE)
  png_data_uri(png_path)
}

# hex 色を rgba() 文字列に変換 (box の半透明塗り用)。
hex_to_rgba <- function(hex, alpha) {
  hex <- sub("^#", "", hex)
  if (nchar(hex) != 6) return(sprintf("rgba(136,136,136,%s)", alpha))
  r <- strtoi(substr(hex, 1, 2), 16L); g <- strtoi(substr(hex, 3, 4), 16L); b <- strtoi(substr(hex, 5, 6), 16L)
  sprintf("rgba(%d,%d,%d,%s)", r, g, b, alpha)
}

# Cluster expression profiles を plotly のインタラクティブ版に変換する。
# cluster ごとにパネル (free_y) を横に並べ、各パネルで group ごとの box + beeswarm。
# 戻り値: list(data, layout, height) もしくは NULL。
build_cluster_profile_plotly <- function(deseq2_dir, project, group_pal = NULL) {
  f <- file.path(degcluster_dir_of(deseq2_dir, project), sprintf("DEGCluster_profile_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  prof <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  if (nrow(prof) == 0 || !all(c("rlog_expr", "cluster_id", "group") %in% colnames(prof))) return(NULL)
  prof$value <- suppressWarnings(as.numeric(prof$rlog_expr))
  prof <- prof[is.finite(prof$value), , drop = FALSE]
  if (nrow(prof) == 0) return(NULL)

  summ_f <- file.path(degcluster_dir_of(deseq2_dir, project), sprintf("DEGCluster_summary_%s.csv", project))
  nmap <- NULL
  if (file.exists(summ_f)) {
    summ <- suppressMessages(readr::read_csv(summ_f, show_col_types = FALSE, progress = FALSE))
    if (all(c("cluster_id", "gene_count") %in% colnames(summ))) nmap <- stats::setNames(summ$gene_count, summ$cluster_id)
  }
  clusters <- unique(as.character(prof$cluster_id))
  clab <- function(cl) if (!is.null(nmap) && cl %in% names(nmap)) sprintf("%s (n=%s)", cl, nmap[[cl]]) else cl
  groups <- unique(as.character(prof$group))
  pal <- ensure_palette(group_pal, groups)
  # group を palette (= sample_sheet) の順に並べる
  if (!is.null(group_pal)) groups <- pic_reorder_vec(groups, names(group_pal))

  ncl <- length(clusters)
  gap <- 0.045
  panelw <- (1 - gap * (ncl - 1)) / ncl

  traces <- list()
  annotations <- list()
  # 凡例を strip ラベルより上に置き、重なりを避ける。boxgap を詰めて箱を太くする。
  layout <- list(margin = list(l = 54, r = 10, t = 74, b = 64), hovermode = "closest",
                 showlegend = TRUE, boxgap = 0.2, boxgroupgap = 0,
                 legend = list(orientation = "h", yanchor = "bottom", y = 1.14, x = 0))
  for (i in seq_len(ncl)) {
    cl <- clusters[[i]]
    x0 <- (i - 1) * (panelw + gap); x1 <- x0 + panelw
    xref <- if (i == 1) "x" else paste0("x", i)
    yref <- if (i == 1) "y" else paste0("y", i)
    xname <- if (i == 1) "xaxis" else paste0("xaxis", i)
    yname <- if (i == 1) "yaxis" else paste0("yaxis", i)
    layout[[xname]] <- list(domain = list(x0, x1), anchor = yref, type = "category",
                            categoryorder = "array", categoryarray = as.list(groups),
                            tickangle = -40, automargin = TRUE, showline = TRUE, mirror = TRUE,
                            linecolor = "#cdd7e2")
    # 全パネル 0 起点 (free-y だが下端は 0)。y=0 の黒線 (zeroline) は消す。
    layout[[yname]] <- list(domain = list(0, 1), anchor = xref, automargin = TRUE, rangemode = "tozero",
                            zeroline = FALSE,
                            title = if (i == 1) list(text = "rlog expression") else list(text = ""),
                            showline = TRUE, mirror = TRUE, linecolor = "#cdd7e2")
    sub <- prof[as.character(prof$cluster_id) == cl, , drop = FALSE]
    for (g in groups) {
      gv <- as.numeric(sub$value[as.character(sub$group) == g])
      if (length(gv) == 0) next
      col <- if (g %in% names(pal)) pal[[g]] else "#888888"
      traces[[length(traces) + 1L]] <- list(
        y = as.list(gv), x = as.list(rep(g, length(gv))),
        type = "box", name = group_span(g, pal), legendgroup = g, showlegend = (i == 1),
        boxpoints = "all", jitter = 1, pointpos = 0, boxmean = FALSE, whiskerwidth = 0.6,
        fillcolor = hex_to_rgba(col, "0.45"),
        line = list(color = col, width = 1.2),
        marker = list(color = col, size = 3, opacity = 0.55),
        xaxis = xref, yaxis = yref,
        hovertemplate = sprintf("%s<br>rlog: %%{y:.2f}<extra></extra>", html_escape(g))
      )
    }
    # strip ラベルは枠なし (テキストのみ)
    annotations[[length(annotations) + 1L]] <- list(
      text = html_escape(clab(cl)), xref = "paper", yref = "paper",
      x = (x0 + x1) / 2, y = 1.0, xanchor = "center", yanchor = "bottom",
      showarrow = FALSE, font = list(size = 13, color = "#1f2933")
    )
  }
  layout$annotations <- annotations
  list(data = traces, layout = layout, height = 380L,
       config = list(responsive = TRUE, displaylogo = FALSE, displayModeBar = TRUE))
}

