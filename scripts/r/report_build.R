# 役割:
#   解析出力ディレクトリ (mapping / deseq2 / enrich) を読み取り、
#   自己完結型の HTML レポートを 1 ファイル生成する補助関数群。
# 入力:
#   <out-dir> 配下の mapping_sum__*.tsv, deseq2/**/stats_*.csv 等。
# 出力:
#   <out-dir>/report_<project>.html (plotly + 画像を base64 で内包)。
# 注記:
#   PCA / heatmap / MA / volcano は plotly でインタラクティブ (hover でラベル)。
#   mapping QC は CSS で 100% 積み上げ棒 + データバー表。enrich plot は再生成して内包。

# ---------------------------------------------------------------------------
# 小物
# ---------------------------------------------------------------------------

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

fmt_int <- function(x) {
  formatC(round(as.numeric(x)), format = "d", big.mark = ",")
}

fmt_ratio <- function(x, digits = 1) {
  formatC(as.numeric(x), format = "f", digits = digits)
}

# plot 仕様を共有レジストリ env に登録する
pic_report_registry <- function() {
  reg <- new.env(parent = emptyenv())
  reg$plots <- list()
  reg
}

register_plot <- function(reg, id, spec) {
  reg$plots[[id]] <- spec
  invisible(id)
}

# ---------------------------------------------------------------------------
# プロジェクト探索
# ---------------------------------------------------------------------------

# <out_dir> 配下から deseq2 の stats_*.csv を探し、プロジェクト単位の記述子を返す。
pic_report_discover_projects <- function(out_dir) {
  stats_files <- list.files(
    out_dir,
    pattern = "^stats_.*\\.csv$",
    recursive = TRUE,
    full.names = TRUE
  )
  # deseq2 ディレクトリ配下のものだけ採用
  stats_files <- stats_files[grepl("deseq2", stats_files, fixed = TRUE)]
  if (length(stats_files) == 0) {
    return(list())
  }

  projects <- list()
  for (sf in stats_files) {
    deseq2_dir <- dirname(sf)
    project <- sub("^stats_", "", basename(sf))
    project <- sub("\\.csv$", "", project)

    enrich_dir <- pic_report_find_enrich_dir(out_dir, project)

    projects[[length(projects) + 1]] <- list(
      project = project,
      deseq2_dir = deseq2_dir,
      enrich_dir = enrich_dir,
      stats_csv = sf
    )
  }
  projects
}

# enrich ディレクトリ (csv/GSEA を含み project 名に一致するもの) を探す
pic_report_find_enrich_dir <- function(out_dir, project) {
  all_csvs <- list.files(out_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  gsea_csvs <- all_csvs[
    grepl("enrich", all_csvs, fixed = TRUE) &
      grepl(file.path("csv", "GSEA"), all_csvs, fixed = TRUE) &
      grepl(paste0(project, ".csv"), basename(all_csvs), fixed = TRUE)
  ]
  if (length(gsea_csvs) == 0) {
    return(NA_character_)
  }
  # enrich/csv/GSEA/<METHOD>/file -> enrich ルートは 3 階層上
  enrich_root <- dirname(dirname(dirname(dirname(gsea_csvs[[1]]))))
  enrich_root
}

# ---------------------------------------------------------------------------
# マッピング QC (CSS) : 100% 積み上げ棒 + データバー表
# ---------------------------------------------------------------------------

# 添付フォーマットの配色 (Office パレット準拠)
PIC_FATE_COLORS <- c(
  trimmed      = "#1F3864",
  unmapped     = "#ED7D31",
  multimapping = "#375623",
  nofeatures   = "#5B9BD5",
  ambiguity    = "#7030A0",
  assigned     = "#70AD47"
)
PIC_FATE_LABELS <- c(
  trimmed = "trimmed", unmapped = "unmapped", multimapping = "multimapping",
  nofeatures = "nofeatures", ambiguity = "ambiguity", assigned = "assigned"
)

read_mapping_sum <- function(path) {
  df <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  df <- as.data.frame(df, check.names = FALSE)
  df
}

section_mapping_qc <- function(msum, section_id = "qc", heading = "1. Mapping QC") {
  fate_cols <- names(PIC_FATE_COLORS)
  have_fate <- all(fate_cols %in% colnames(msum))
  parts <- c(sprintf('<section id="%s"><h2>%s</h2>', section_id, heading))

  # ---- 100% 積み上げ棒 ----
  if (have_fate) {
    parts <- c(parts, '<h3>Read distribution</h3>')
    # 目盛り (0-100%)
    ruler <- '<div class="pic-bar-ruler"><span>0%</span><span>20%</span><span>40%</span><span>60%</span><span>80%</span><span>100%</span></div>'
    rows <- character(0)
    for (i in seq_len(nrow(msum))) {
      sample <- as.character(msum$sample[[i]])
      vals <- vapply(fate_cols, function(cc) suppressWarnings(as.numeric(msum[[cc]][[i]])), numeric(1))
      vals[!is.finite(vals)] <- 0
      tot <- sum(vals)
      if (tot <= 0) tot <- 1
      segs <- character(0)
      for (cc in fate_cols) {
        pct <- 100 * vals[[cc]] / tot
        if (pct <= 0) next
        tip <- sprintf("%s: %s (%.1f%%)", PIC_FATE_LABELS[[cc]], fmt_int(vals[[cc]]), pct)
        segs <- c(segs, sprintf(
          '<div class="pic-seg" style="width:%.4f%%;background:%s" title="%s — %s"></div>',
          pct, PIC_FATE_COLORS[[cc]], html_escape(sample), html_escape(tip)
        ))
      }
      rows <- c(rows, sprintf(
        '<div class="pic-bar-row"><div class="pic-bar-label">%s</div><div class="pic-bar-track">%s</div></div>',
        html_escape(sample), paste(segs, collapse = "")
      ))
    }
    legend_items <- vapply(fate_cols, function(cc) sprintf(
      '<span class="pic-legend-item"><span class="pic-swatch" style="background:%s"></span>%s</span>',
      PIC_FATE_COLORS[[cc]], PIC_FATE_LABELS[[cc]]
    ), character(1))
    parts <- c(parts,
      '<div class="pic-bars">',
      ruler,
      paste(rows, collapse = "\n"),
      '</div>',
      sprintf('<div class="pic-legend">%s</div>', paste(legend_items, collapse = ""))
    )
  }

  # ---- データバー表 ----
  parts <- c(parts, '<h3>Sequencing depth</h3>')
  bar_cols <- list(
    total = list(label = "total", grad = c("#9DC3E6", "#D9E7F5")),
    umis  = list(label = "umis",  grad = c("#63C384", "#D6EFDD")),
    genes = list(label = "genes", grad = c("#FFC000", "#FFE9AE"))
  )
  plain_cols <- list()
  if ("umis/genes" %in% colnames(msum)) plain_cols[["umis/genes"]] <- "umis/genes"
  if ("assigned/umis" %in% colnames(msum)) plain_cols[["assigned/umis"]] <- "assigned/umi"

  maxima <- lapply(names(bar_cols), function(cc) {
    if (cc %in% colnames(msum)) max(suppressWarnings(as.numeric(msum[[cc]])), na.rm = TRUE) else NA_real_
  })
  names(maxima) <- names(bar_cols)

  thead <- '<tr><th>sample</th>'
  for (cc in names(bar_cols)) if (cc %in% colnames(msum)) thead <- paste0(thead, sprintf("<th>%s</th>", bar_cols[[cc]]$label))
  for (cc in names(plain_cols)) thead <- paste0(thead, sprintf("<th>%s</th>", plain_cols[[cc]]))
  thead <- paste0(thead, "</tr>")

  trows <- character(0)
  for (i in seq_len(nrow(msum))) {
    cells <- sprintf('<td class="pic-td-sample">%s</td>', html_escape(as.character(msum$sample[[i]])))
    for (cc in names(bar_cols)) {
      if (!(cc %in% colnames(msum))) next
      v <- suppressWarnings(as.numeric(msum[[cc]][[i]]))
      if (!is.finite(v)) v <- 0
      mx <- maxima[[cc]]
      pct <- if (is.finite(mx) && mx > 0) 100 * v / mx else 0
      grad <- bar_cols[[cc]]$grad
      bar <- sprintf(
        'background:linear-gradient(90deg,%s 0%%,%s %.3f%%,transparent %.3f%%)',
        grad[[1]], grad[[2]], pct, pct
      )
      cells <- paste0(cells, sprintf(
        '<td class="pic-td-bar"><div class="pic-databar" style="%s"></div><span class="pic-databar-val">%s</span></td>',
        bar, fmt_int(v)
      ))
    }
    for (cc in names(plain_cols)) {
      v <- suppressWarnings(as.numeric(msum[[cc]][[i]]))
      cells <- paste0(cells, sprintf('<td class="pic-td-num">%s</td>', if (is.finite(v)) fmt_ratio(v) else "NA"))
    }
    trows <- c(trows, sprintf("<tr>%s</tr>", cells))
  }

  parts <- c(parts,
    '<table class="pic-table pic-databar-table"><thead>',
    thead,
    '</thead><tbody>',
    paste(trows, collapse = "\n"),
    '</tbody></table>'
  )

  parts <- c(parts, '</section>')
  paste(parts, collapse = "\n")
}

# ---------------------------------------------------------------------------
# DESeq2: PCA / heatmap (plotly)
# ---------------------------------------------------------------------------

pca_variance_pct <- function(deseq2_dir, project) {
  f <- file.path(deseq2_dir, "PCA", sprintf("PCA_Variance_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  v <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  v <- as.data.frame(v, check.names = FALSE)
  row <- v[v[[1]] == "Proportion of Variance", , drop = FALSE]
  if (nrow(row) == 0) return(NULL)
  pcs <- setdiff(colnames(row), colnames(row)[[1]])
  out <- as.numeric(row[1, pcs]) * 100
  names(out) <- pcs
  out
}

group_palette <- function(groups) {
  pal <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
           "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#393b79", "#637939",
           "#8c6d31", "#843c39", "#7b4173", "#5254a3", "#8ca252", "#bd9e39")
  ug <- unique(groups)
  cols <- pal[((seq_along(ug) - 1) %% length(pal)) + 1]
  stats::setNames(cols, ug)
}

build_pca_plots <- function(reg, deseq2_dir, project, group_pal = NULL, id_prefix = "") {
  f <- file.path(deseq2_dir, "PCA", sprintf("PCA_RegLog_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  d <- as.data.frame(d, check.names = FALSE)
  if (!("group" %in% colnames(d))) d$group <- "all"
  varpct <- pca_variance_pct(deseq2_dir, project)
  pal <- if (!is.null(group_pal)) group_pal else group_palette(d$group)
  # パレットに無いグループがあれば補完
  miss <- setdiff(unique(d$group), names(pal))
  if (length(miss) > 0) pal <- c(pal, group_palette(miss))

  # ---- 寄与率 (scree) ----
  if (!is.null(varpct)) {
    pcs <- names(varpct)
    npc <- min(length(pcs), 15L)
    pcs <- pcs[seq_len(npc)]
    cum <- cumsum(varpct[pcs])
    scree <- list(
      list(x = as.list(pcs), y = as.list(unname(varpct[pcs])), type = "bar",
           name = "variance", marker = list(color = "#5B9BD5"),
           hovertemplate = "%{x}: %{y:.1f}%<extra></extra>"),
      list(x = as.list(pcs), y = as.list(unname(cum)), type = "scatter", mode = "lines+markers",
           name = "cumulative", line = list(color = "#d7301f"),
           hovertemplate = "%{x} cumulative: %{y:.1f}%<extra></extra>")
    )
    register_plot(reg, paste0(id_prefix, "pca_scree"), list(data = scree, layout = list(
      xaxis = list(title = ""),
      yaxis = list(title = "% variance", rangemode = "tozero"),
      hovermode = "x", margin = list(l = 56, r = 16, t = 28, b = 40),
      legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0)
    )))
  }
  scree_block <- sprintf('<div class="pic-plot-cell"><h4>Variance explained (scree)</h4><div id="%spca_scree" class="pic-plot"></div></div>', id_prefix)

  pair_blocks <- character(0)
  pairs <- list(c("PC1", "PC2"), c("PC2", "PC3"))
  for (pr in pairs) {
    xc <- pr[[1]]; yc <- pr[[2]]
    if (!all(c(xc, yc) %in% colnames(d))) next
    traces <- list()
    for (g in unique(d$group)) {
      sub <- d[d$group == g, , drop = FALSE]
      traces[[length(traces) + 1]] <- list(
        x = as.list(as.numeric(sub[[xc]])),
        y = as.list(as.numeric(sub[[yc]])),
        text = as.list(as.character(sub$sample)),
        name = g,
        mode = "markers",
        type = "scatter",
        marker = list(size = 11, color = unname(pal[[g]]), line = list(width = 1, color = "#ffffff")),
        hovertemplate = paste0("%{text}<br>", xc, ": %{x:.2f}<br>", yc, ": %{y:.2f}<extra>", html_escape(g), "</extra>")
      )
    }
    xlab <- if (!is.null(varpct) && xc %in% names(varpct)) sprintf("%s (%.1f%%)", xc, varpct[[xc]]) else xc
    ylab <- if (!is.null(varpct) && yc %in% names(varpct)) sprintf("%s (%.1f%%)", yc, varpct[[yc]]) else yc
    id <- sprintf("%spca_%s_%s", id_prefix, xc, yc)
    layout <- list(
      xaxis = list(title = xlab, zeroline = TRUE),
      yaxis = list(title = ylab, zeroline = TRUE),
      hovermode = "closest",
      margin = list(l = 60, r = 20, t = 28, b = 46),
      legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0)
    )
    register_plot(reg, id, list(data = traces, layout = layout))
    pair_blocks <- c(pair_blocks, sprintf('<div class="pic-plot-cell"><h4>%s vs %s</h4><div id="%s" class="pic-plot"></div></div>', xc, yc, id))
  }
  if (is.null(varpct) && length(pair_blocks) == 0) return(NULL)
  # 1 行目: 寄与率のみ / 2 行目: PC1-2, PC2-3
  paste0(
    if (!is.null(varpct)) sprintf('<div class="pic-plot-grid pic-pca-scree-row">%s</div>', scree_block) else "",
    if (length(pair_blocks) > 0) sprintf('<div class="pic-plot-grid pic-pca-pair-row">%s</div>', paste(pair_blocks, collapse = "")) else ""
  )
}

# z 値 -> 発散カラー (青-白-赤) の hex
zcolor_vec <- function(v, zmax) {
  ramp <- grDevices::colorRamp(c("#2166ac", "#ffffff", "#b2182b"))
  t <- (v + zmax) / (2 * zmax)
  t[!is.finite(t)] <- 0.5
  t <- pmin(1, pmax(0, t))
  rgb <- ramp(t)
  grDevices::rgb(rgb[, 1], rgb[, 2], rgb[, 3], maxColorValue = 255)
}

# DEG ヒートマップを、サンプル行ヘッダ固定・縦スクロール可能な HTML 表として生成する。
build_heatmap_html <- function(deseq2_dir, project, group_map = NULL, group_pal = NULL) {
  f <- file.path(deseq2_dir, "DEG", sprintf("DEG_normalizedCountTable_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  d <- as.data.frame(d, check.names = FALSE)
  ann <- c("ens_gene", "ext_gene", "biotype", "chr")
  sample_cols <- setdiff(colnames(d), ann)
  if (length(sample_cols) < 2 || nrow(d) < 2) return(NULL)

  labels <- if ("ext_gene" %in% colnames(d)) {
    lab <- as.character(d$ext_gene)
    lab[is.na(lab) | lab == ""] <- as.character(d$ens_gene)[is.na(lab) | lab == ""]
    lab
  } else {
    as.character(d$ens_gene)
  }
  mat <- as.matrix(d[, sample_cols, drop = FALSE])
  storage.mode(mat) <- "double"
  z <- t(scale(t(mat)))                 # 行方向 z-score
  z[!is.finite(z)] <- 0

  # 行のみ hclust で並べ替え。列 (サンプル) は元の順番を維持。
  order_rows <- tryCatch(stats::hclust(stats::dist(z))$order, error = function(e) seq_len(nrow(z)))
  z <- z[order_rows, , drop = FALSE]
  labels <- labels[order_rows]
  scols <- sample_cols

  zmax <- max(1, stats::quantile(abs(z), 0.98, names = FALSE, na.rm = TRUE))
  colmat <- matrix(zcolor_vec(as.numeric(z), zmax), nrow = nrow(z))

  # ヘッダ (サンプル名・グループ色, sticky)
  ths <- vapply(seq_along(scols), function(j) {
    s <- scols[[j]]
    g <- if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s
    col <- if (!is.null(group_pal) && g %in% names(group_pal)) unname(group_pal[[g]]) else "#333333"
    sprintf('<th class="s"><span style="color:%s">%s</span></th>', col, html_escape(s))
  }, character(1))
  header <- sprintf('<tr><th class="corner"></th>%s</tr>', paste(ths, collapse = ""))

  rows <- vapply(seq_len(nrow(z)), function(i) {
    tds <- vapply(seq_len(ncol(z)), function(j) sprintf(
      '<td style="background:%s" title="Gene: %s&#10;Sample: %s&#10;Z: %.2f"></td>',
      colmat[i, j], html_escape(labels[[i]]), html_escape(scols[[j]]), z[i, j]
    ), character(1))
    sprintf('<tr><th class="g">%s</th>%s</tr>', html_escape(labels[[i]]), paste(tds, collapse = ""))
  }, character(1))

  legend_html <- ""
  if (!is.null(group_pal)) {
    items <- vapply(names(group_pal), function(g) sprintf(
      '<span class="pic-legend-item"><span class="pic-swatch" style="background:%s"></span>%s</span>',
      unname(group_pal[[g]]), html_escape(g)
    ), character(1))
    legend_html <- sprintf('<div class="pic-legend" style="margin-left:0">%s</div>', paste(items, collapse = ""))
  }

  sprintf(
    '<h3>Heatmap (row z-score, %d genes)</h3><p class="pic-note">Hover a cell to show gene / sample / z.</p>%s<div class="pic-hm-scroll"><table class="pic-hm"><thead>%s</thead><tbody>%s</tbody></table></div>',
    nrow(z), legend_html, header, paste(rows, collapse = "")
  )
}

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

scatter_traces_by_dir <- function(df, x, y, numerator, denominator, hovertemplate) {
  # df は idx で絞り込み済み。direction 列 (up/down/ns) でトレース分割。
  # df には cd1 (padj), cd2 (pvalue) 列を含み、customdata として hover に渡す。
  cats <- list(
    up = list(name = sprintf("%s ↑", numerator), color = "#d7301f"),
    down = list(name = sprintf("%s ↑", denominator), color = "#2166ac"),
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
      marker = list(size = if (k == "ns") 4 else 6, color = cats[[k]]$color,
                    opacity = if (k == "ns") 0.45 else 0.85),
      hovertemplate = hovertemplate
    )
  }
  traces
}

build_contrast_plots <- function(reg, stats, deg_counts, fdr, id_prefix = "") {
  aspects <- unique(stats$aspect)
  # contrast を DEG 数で降順に並べる
  deg_per <- vapply(aspects, function(a) {
    if (!is.null(deg_counts) && a %in% names(deg_counts)) deg_counts[[a]] else {
      sum(stats$aspect == a & !is.na(stats$padj) & stats$padj < fdr)
    }
  }, numeric(1))
  ord <- order(deg_per, decreasing = TRUE)
  aspects <- aspects[ord]
  deg_per <- deg_per[ord]

  blocks <- character(0)
  for (k in seq_along(aspects)) {
    a <- aspects[[k]]
    nd <- deg_per[[k]]
    cs <- stats[stats$aspect == a, , drop = FALSE]
    sp <- strsplit(a, " / ", fixed = TRUE)[[1]]
    numerator <- if (length(sp) >= 1) sp[[1]] else "num"
    denominator <- if (length(sp) >= 2) sp[[2]] else "den"
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
    meta_html <- sprintf('<p class="pic-contrast-meta">DEG (padj &lt; %.2g): <b>%d</b></p>', fdr, as.integer(nd))
    cells <- character(0)

    # ---- MA (hover に padj) ----
    keepMA <- is.finite(base) & base > 0 & is.finite(lfc)
    idxMA <- downsample_idx(is_sig & keepMA)
    idxMA <- idxMA[keepMA[idxMA]]
    if (length(idxMA) > 0) {
      dfMA <- data.frame(x = base[idxMA], y = lfc[idxMA], gene = gene[idxMA], dir = dir[idxMA],
                         cd1 = padj[idxMA], cd2 = pval[idxMA], stringsAsFactors = FALSE)
      ht <- "<b>%{text}</b><br>baseMean: %{x:.3g}<br>log2FC: %{y:.2f}<br>padj: %{customdata[0]:.3g}<extra></extra>"
      tr <- scatter_traces_by_dir(dfMA, "x", "y", numerator, denominator, ht)
      id <- sprintf("%sma_%s", id_prefix, flab)
      layout <- list(
        xaxis = list(title = "baseMean", type = "log"),
        yaxis = list(title = "log2 fold change", zeroline = TRUE),
        hovermode = "closest", margin = list(l = 60, r = 20, t = 28, b = 46),
        legend = legend_top
      )
      register_plot(reg, id, list(data = tr, layout = layout))
      cells <- c(cells, sprintf('<div class="pic-plot-cell"><h4>MA</h4><div id="%s" class="pic-plot"></div></div>', id))
    }

    # ---- volcano: y=-log10(pvalue) 一本化、hover に padj も表示 ----
    keepV <- is.finite(pval) & pval > 0 & is.finite(lfc)
    idxV <- downsample_idx(is_sig & keepV)
    idxV <- idxV[keepV[idxV]]
    if (length(idxV) > 0) {
      dfV <- data.frame(x = lfc[idxV], y = -log10(pval[idxV]), gene = gene[idxV], dir = dir[idxV],
                        cd1 = padj[idxV], cd2 = pval[idxV], stringsAsFactors = FALSE)
      ht <- "<b>%{text}</b><br>log2FC: %{x:.2f}<br>pvalue: %{customdata[1]:.3g}<br>padj: %{customdata[0]:.3g}<extra></extra>"
      tr <- scatter_traces_by_dir(dfV, "x", "y", numerator, denominator, ht)
      id <- sprintf("%svolcano_%s", id_prefix, flab)
      layout <- list(
        xaxis = list(title = "log2 fold change", zeroline = TRUE),
        yaxis = list(title = "-log10(pvalue)"),
        hovermode = "closest", margin = list(l = 60, r = 20, t = 28, b = 46),
        legend = legend_top
      )
      register_plot(reg, id, list(data = tr, layout = layout))
      cells <- c(cells, sprintf('<div class="pic-plot-cell"><h4>volcano</h4><div id="%s" class="pic-plot"></div></div>', id))
    }

    open_attr <- if (k <= 3) " open" else ""
    blocks <- c(blocks, sprintf(
      '<details class="pic-contrast"%s><summary>%s — DEG %d</summary>%s<div class="pic-ma-vol">%s</div></details>',
      open_attr, html_escape(a), as.integer(nd), meta_html, paste(cells, collapse = "")
    ))
  }
  paste0('<h3>MA / volcano</h3>',
         '<p class="pic-note">Volcano y-axis is -log10(pvalue); point color indicates significance by padj (hover shows padj).</p>',
         paste(blocks, collapse = "\n"))
}

# ---------------------------------------------------------------------------
# enrich: GSEA (contrast ごと) / ORA (method ごと) を再生成して内包
# ---------------------------------------------------------------------------

png_data_uri <- function(path) {
  if (!file.exists(path)) return(NA_character_)
  raw <- readBin(path, what = "raw", n = file.info(path)$size)
  paste0("data:image/png;base64,", base64enc::base64encode(raw))
}

# <base_dir>/aggregate/<genome>/aggregate_profile.csv からインタラクティブな
# メタジーン折れ線 (plotly) を作る。CSV が無ければ profile.png を埋め込む。
build_aggregate_html <- function(reg, base_dir, id_prefix = "") {
  agg_root <- file.path(base_dir, "aggregate")
  if (!dir.exists(agg_root)) return("")
  gdirs <- list.dirs(agg_root, recursive = FALSE, full.names = TRUE)
  blocks <- character(0)
  for (gd in gdirs) {
    genome <- basename(gd)
    csv <- file.path(gd, "aggregate_profile.csv")
    if (file.exists(csv)) {
      d <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE, progress = FALSE))
      d <- as.data.frame(d, check.names = FALSE)
      if (nrow(d) == 0 || !all(c("sample", "group", "pos", "value") %in% colnames(d))) next
      pal <- group_palette(unique(as.character(d$group)))
      traces <- list()
      for (s in unique(as.character(d$sample))) {
        sub <- d[as.character(d$sample) == s, , drop = FALSE]
        sub <- sub[order(sub$pos), , drop = FALSE]
        g <- as.character(sub$group[[1]])
        traces[[length(traces) + 1]] <- list(
          x = as.list(as.numeric(sub$pos)), y = as.list(as.numeric(sub$value)),
          name = s, legendgroup = g, mode = "lines", type = "scatter",
          line = list(color = unname(pal[[g]]), width = 1.4),
          hovertemplate = paste0("%{fullData.name}<br>CPM: %{y:.2f}<extra>", html_escape(g), "</extra>"))
      }
      id <- sprintf("%saggregate_%s", id_prefix, genome)
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
        hovermode = "x unified", margin = list(l = 60, r = 20, t = 30, b = 46),
        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0))
      register_plot(reg, id, list(data = traces, layout = layout))
      blocks <- c(blocks, sprintf(
        '<div class="pic-plot-cell"><h4>%s — gene-body profile (TSS→TES)</h4><div id="%s" class="pic-plot"></div></div>',
        html_escape(genome), id))
    } else {
      png <- file.path(gd, "aggregate_profile.png")
      if (file.exists(png)) {
        uri <- png_data_uri(png)
        if (!is.na(uri)) blocks <- c(blocks, sprintf(
          '<div class="pic-enrich-item"><h4>%s — gene-body profile (TSS→TES)</h4><img loading="lazy" src="%s"></div>',
          html_escape(genome), uri))
      }
    }
  }
  if (length(blocks) == 0) return("")
  paste0('<p class="pic-note">Mean CPM over the scaled gene body (TSS→TES). Hover for sample/position/value; drag to zoom.</p><div class="pic-plot-grid">',
         paste(blocks, collapse = ""), '</div>')
}

# DEG クラスタの挙動 (group ごとの rlog 分布) を facet で示す図。
# 各 cluster_N が何を意味するか (どの group で高い/低い) を ORA の前に提示する。
build_cluster_profile_png <- function(deseq2_dir, project, tmp_dir, group_pal = NULL) {
  f <- file.path(deseq2_dir, "DEG", "DEGCluster", sprintf("DEGCluster_profile_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  prof <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  if (nrow(prof) == 0 || !all(c("rlog_expr", "cluster_id", "group") %in% colnames(prof))) return(NULL)

  # cluster ラベルに遺伝子数を付与
  summ_f <- file.path(deseq2_dir, "DEG", "DEGCluster", sprintf("DEGCluster_summary_%s.csv", project))
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

# enrichment セクションの中身 (h3 GSEA / ORA) を返す。<section> ラッパは付けない。
enrichment_blocks <- function(enrich_dir, project, tmp_dir, deg_counts = NULL, deseq2_dir = NULL, group_pal = NULL) {
  if (is.null(enrich_dir) || is.na(enrich_dir) || !dir.exists(enrich_dir)) return("")
  ecfg <- pic_plot_spec()$plot$enrichment
  parts <- character(0)

  method_order <- strsplit(pic_plot_spec()$defaults$enrich_methods_csv, ",", fixed = TRUE)[[1]]
  order_methods <- function(ms) {
    c(method_order[method_order %in% ms], sort(setdiff(ms, method_order)))
  }

  # ---- GSEA: contrast ごとにグループ化し、配下に method (GO_BP, KEGG, ...) ----
  gsea_root <- file.path(enrich_dir, "csv", "GSEA")
  gsea_blocks <- character(0)
  if (dir.exists(gsea_root)) {
    mdirs <- list.dirs(gsea_root, recursive = FALSE, full.names = TRUE)
    # records: 各 (method, contrast) -> png uri
    by_contrast <- list()
    for (mdir in mdirs) {
      method <- basename(mdir)
      csvs <- list.files(mdir, pattern = paste0(project, "\\.csv$"), full.names = TRUE)
      for (csv in sort(csvs)) {
        df <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE, progress = FALSE))
        if (nrow(df) == 0 || !("contrast" %in% colnames(df))) next
        contrast <- as.character(df$contrast[[1]])
        sp <- strsplit(contrast, " / ", fixed = TRUE)[[1]]
        numerator <- if (length(sp) >= 1) sp[[1]] else ""
        denominator <- if (length(sp) >= 2) sp[[2]] else ""
        p <- tryCatch(build_gsea_plot(df, numerator, denominator), error = function(e) NULL)
        if (is.null(p)) next
        p <- p + ggplot2::theme(plot.margin = ggplot2::margin(4, 6, 4, 4))
        n_labels <- suppressWarnings(as.numeric(attr(p, "pic_n_y_labels")))
        h <- if (is.finite(n_labels)) max(2.2, min(11, 0.9 + n_labels * 0.16)) else 5
        png_path <- file.path(tmp_dir, sprintf("gsea_%s_%s.png", method, format_contrast_file_label(contrast)))
        ggplot2::ggsave(png_path, plot = p, width = 9, height = h, dpi = 120, limitsize = FALSE)
        uri <- png_data_uri(png_path)
        if (is.na(uri)) next
        rec <- list(method = method, uri = uri)
        by_contrast[[contrast]] <- c(by_contrast[[contrast]], list(rec))
      }
    }
    # contrast の並び順: DEG 数降順 (なければ名前順)
    contrasts <- names(by_contrast)
    if (!is.null(deg_counts)) {
      keyv <- vapply(contrasts, function(a) if (a %in% names(deg_counts)) deg_counts[[a]] else 0, numeric(1))
      contrasts <- contrasts[order(keyv, decreasing = TRUE)]
    } else {
      contrasts <- sort(contrasts)
    }
    for (contrast in contrasts) {
      recs <- by_contrast[[contrast]]
      ms <- vapply(recs, function(r) r$method, character(1))
      ord <- match(order_methods(ms), ms)
      recs <- recs[ord[!is.na(ord)]]
      imgs <- vapply(recs, function(r) sprintf(
        '<div class="pic-enrich-item"><h4>%s</h4><img loading="lazy" src="%s"></div>',
        html_escape(r$method), r$uri
      ), character(1))
      open_attr <- if (length(gsea_blocks) < 2) " open" else ""
      gsea_blocks <- c(gsea_blocks, sprintf(
        '<details class="pic-enrich-method"%s><summary>%s (%d method)</summary>%s</details>',
        open_attr, html_escape(contrast), length(imgs), paste(imgs, collapse = "")
      ))
    }
  }
  if (length(gsea_blocks) > 0) {
    parts <- c(parts, sprintf(
      '<details class="pic-enrich-group" open><summary>GSEA</summary><p class="pic-note">Expand a contrast to see GSEA dot plots for GO_BP / KEGG / REACTOME ...</p>%s</details>',
      paste(gsea_blocks, collapse = "\n")
    ))
  }

  # ---- ORA: method ごと (cluster をまとめて 1 枚) ----
  ora_root <- file.path(enrich_dir, "csv", "ORA")
  ora_blocks <- character(0)
  if (dir.exists(ora_root)) {
    methods <- list.dirs(ora_root, recursive = FALSE, full.names = TRUE)
    for (mdir in methods) {
      method <- basename(mdir)
      csvs <- list.files(mdir, pattern = paste0(project, "\\.csv$"), full.names = TRUE)
      if (length(csvs) == 0) next
      ora_all <- tryCatch(
        purrr::map_dfr(csvs, function(fp) suppressMessages(readr::read_csv(fp, show_col_types = FALSE, progress = FALSE))),
        error = function(e) NULL
      )
      if (is.null(ora_all) || nrow(ora_all) == 0) next
      p <- tryCatch(build_cluster_ora_plot(ora_all), error = function(e) NULL)
      if (is.null(p)) next
      p <- p + ggplot2::theme(plot.margin = ggplot2::margin(4, 6, 4, 4))
      n_labels <- suppressWarnings(as.numeric(attr(p, "pic_n_y_labels")))
      h <- if (is.finite(n_labels)) max(2.5, min(13, 0.9 + n_labels * 0.17)) else 7
      png_path <- file.path(tmp_dir, sprintf("ora_%s.png", method))
      ggplot2::ggsave(png_path, plot = p, width = 9, height = h, dpi = 120, limitsize = FALSE)
      uri <- png_data_uri(png_path)
      if (is.na(uri)) next
      ora_blocks <- c(ora_blocks, sprintf(
        '<div class="pic-enrich-item"><h4>ORA — %s</h4><img loading="lazy" src="%s"></div>',
        html_escape(method), uri
      ))
    }
  }
  if (length(ora_blocks) > 0) {
    ora_inner <- character(0)
    # Cluster expression profiles first (so cluster_N labels are interpretable)
    if (!is.null(deseq2_dir)) {
      prof_uri <- tryCatch(build_cluster_profile_png(deseq2_dir, project, tmp_dir, group_pal), error = function(e) NULL)
      if (!is.null(prof_uri) && !is.na(prof_uri)) {
        ora_inner <- c(ora_inner, sprintf(
          '<div class="pic-enrich-item"><h4>Cluster expression profiles (per-cluster rlog by group)</h4><img loading="lazy" src="%s"></div>',
          prof_uri
        ))
      }
    }
    ora_inner <- c(ora_inner,
      sprintf('<details class="pic-enrich-method" open><summary>ORA dot plots (%d method)</summary>%s</details>',
              length(ora_blocks), paste(ora_blocks, collapse = "")))
    parts <- c(parts, sprintf(
      '<details class="pic-enrich-group" open><summary>ORA (DEG clusters)</summary>%s</details>',
      paste(ora_inner, collapse = "")
    ))
  }

  paste(parts, collapse = "\n")
}

# 通常レポート用の enrichment セクション (<section id="enrich"> でラップ)。
section_enrichment <- function(enrich_dir, project, tmp_dir, deg_counts = NULL, deseq2_dir = NULL, group_pal = NULL, heading = "4. Enrichment") {
  inner <- enrichment_blocks(enrich_dir, project, tmp_dir, deg_counts, deseq2_dir, group_pal)
  if (!nzchar(inner)) inner <- '<p>No enrichment plots were generated.</p>'
  sprintf('<section id="enrich"><h2>%s</h2>%s</section>', heading, inner)
}

# ---------------------------------------------------------------------------
# CSS / JS
# ---------------------------------------------------------------------------

report_css <- function() {
  paste(readLines(file.path(pic_report_asset_dir(), "report.css"), warn = FALSE), collapse = "\n")
}

# ---------------------------------------------------------------------------
# メイン: 1 プロジェクトの HTML を生成
# ---------------------------------------------------------------------------

pic_report_asset_dir <- function() {
  getOption("pic.report.asset_dir", default = ".")
}

build_report_for_project <- function(desc, out_dir, msum, asset_dir) {
  options(pic.report.asset_dir = asset_dir)
  fdr <- pic_plot_spec()$defaults$fdr
  project <- desc$project
  reg <- pic_report_registry()

  stats <- suppressMessages(readr::read_csv(desc$stats_csv, show_col_types = FALSE, progress = FALSE))
  stats <- as.data.frame(stats, check.names = FALSE)

  # DEG 数 (contrast -> 合計) を DEG_count.csv から取得 (なければ stats から)
  deg_counts <- pic_load_deg_counts(desc$deseq2_dir, project)
  if (is.null(deg_counts)) {
    deg_counts <- list()
    for (a in unique(stats$aspect)) {
      deg_counts[[a]] <- sum(stats$aspect == a & !is.na(stats$padj) & stats$padj < fdr)
    }
  }

  tmp_dir <- file.path(tempdir(), paste0("picreport_", project))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  # サンプル -> グループの対応とグループ配色 (PCA / heatmap で共有)
  group_map <- NULL
  group_pal <- NULL
  if (!is.null(msum) && all(c("sample", "group") %in% colnames(msum))) {
    group_map <- stats::setNames(as.character(msum$group), as.character(msum$sample))
    group_pal <- group_palette(as.character(msum$group))
  }

  # セクション組み立て (順序: QC -> Aggregation -> PCA -> DEG -> Enrichment)
  sec_qc <- if (!is.null(msum)) section_mapping_qc(msum) else ""

  # 2. Aggregation (TSS-TES)
  agg_html <- build_aggregate_html(reg, out_dir)
  sec_agg <- if (nzchar(agg_html)) paste0('<section id="aggregate"><h2>2. Aggregation (TSS-TES)</h2>', agg_html, '</section>') else ""

  # 3. PCA
  pca_html <- build_pca_plots(reg, desc$deseq2_dir, project, group_pal)
  sec_pca <- paste0('<section id="pca"><h2>3. PCA</h2>',
                    if (!is.null(pca_html)) pca_html else "<p>No PCA data.</p>",
                    '</section>')

  # 4. DEG (heatmap + MA/volcano)
  hm_html <- build_heatmap_html(desc$deseq2_dir, project, group_map, group_pal)
  contrast_html <- build_contrast_plots(reg, stats, deg_counts, fdr)
  sec_deg <- paste0(sprintf('<section id="deg"><h2>4. DEG (DESeq2; FDR = %s)</h2>', format(fdr, trim = TRUE)),
                    if (!is.null(hm_html)) hm_html else "",
                    contrast_html,
                    '</section>')

  # 5. Enrichment
  sec_enrich <- section_enrichment(desc$enrich_dir, project, tmp_dir, deg_counts, desc$deseq2_dir, group_pal, "5. Enrichment")
  if (is.null(sec_enrich)) sec_enrich <- ""

  nav <- paste0('<nav class="pic-nav"><a href="#qc">QC</a>',
                if (nzchar(sec_agg)) '<a href="#aggregate">Aggregation</a>' else "",
                '<a href="#pca">PCA</a><a href="#deg">DEG</a><a href="#enrich">Enrichment</a></nav>')
  body <- paste0(sec_qc, sec_agg, sec_pca, sec_deg, sec_enrich)
  out_html <- file.path(out_dir, sprintf("report_%s.html", project))
  render_report_page(project, nav, body, reg, asset_dir, out_html)
  unlink(tmp_dir, recursive = TRUE)
  out_html
}

# HTML ページを組み立てて書き出す (plotly を内包)。
render_report_page <- function(title, nav_html, body_html, reg, asset_dir, out_html) {
  plots_json <- jsonlite::toJSON(reg$plots, auto_unbox = TRUE, null = "null", na = "null", digits = 6)
  plotly_js <- paste(readLines(file.path(asset_dir, "plotly.min.js"), warn = FALSE), collapse = "\n")
  page <- paste0(
    '<!DOCTYPE html><html lang="ja"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>pic report — ', html_escape(title), '</title>',
    '<style>', report_css(), '</style>',
    '</head><body>',
    '<header class="pic-header"><h1>pic analysis report</h1>',
    '<div class="pic-sub">project: <b>', html_escape(title), '</b></div></header>',
    nav_html,
    '<main>', body_html, '</main>',
    '<script>', plotly_js, '</script>',
    '<script>var PIC_PLOTS=', plots_json, ';</script>',
    '<script>', report_runtime_js(), '</script>',
    '</body></html>'
  )
  writeLines(page, out_html, useBytes = TRUE)
  out_html
}

# ---------------------------------------------------------------------------
# xenograft 統合レポート (分類 QC + graft/host 2 画分を 1 ファイルに)
# ---------------------------------------------------------------------------

# 分類カテゴリの配色 (host/graft/both/neither/ambiguous)
PIC_XENO_COLORS <- c(
  graft = "#70AD47", host = "#ED7D31", both = "#5B9BD5",
  ambiguous = "#7030A0", neither = "#A6A6A6"
)

# 分類サマリ TSV から 100% 積み上げ棒の QC セクションを作る。
section_xenograft_qc <- function(summary_df) {
  cats <- names(PIC_XENO_COLORS)
  have <- all(cats %in% colnames(summary_df))
  parts <- c('<section id="classification"><h2>1. Xenograft classification (xengsort)</h2>')
  if (!have) {
    parts <- c(parts, '<p>classification summary not found.</p></section>')
    return(paste(parts, collapse = "\n"))
  }
  ruler <- '<div class="pic-bar-ruler"><span>0%</span><span>20%</span><span>40%</span><span>60%</span><span>80%</span><span>100%</span></div>'
  rows <- character(0)
  for (i in seq_len(nrow(summary_df))) {
    sample <- as.character(summary_df$sample[[i]])
    vals <- vapply(cats, function(cc) suppressWarnings(as.numeric(summary_df[[cc]][[i]])), numeric(1))
    vals[!is.finite(vals)] <- 0
    tot <- sum(vals); if (tot <= 0) tot <- 1
    segs <- character(0)
    for (cc in cats) {
      pct <- 100 * vals[[cc]] / tot
      if (pct <= 0) next
      tip <- sprintf("%s: %s (%.1f%%)", cc, fmt_int(vals[[cc]]), pct)
      segs <- c(segs, sprintf('<div class="pic-seg" style="width:%.4f%%;background:%s" title="%s — %s"></div>',
                              pct, PIC_XENO_COLORS[[cc]], html_escape(sample), html_escape(tip)))
    }
    rows <- c(rows, sprintf('<div class="pic-bar-row"><div class="pic-bar-label">%s</div><div class="pic-bar-track">%s</div></div>',
                            html_escape(sample), paste(segs, collapse = "")))
  }
  legend <- vapply(cats, function(cc) sprintf('<span class="pic-legend-item"><span class="pic-swatch" style="background:%s"></span>%s</span>',
                                              PIC_XENO_COLORS[[cc]], cc), character(1))
  parts <- c(parts,
    '<h3>Read classification (per sample, 100% stacked)</h3>',
    '<div class="pic-bars">', ruler, paste(rows, collapse = "\n"), '</div>',
    sprintf('<div class="pic-legend">%s</div>', paste(legend, collapse = "")),
    '</section>')
  paste(parts, collapse = "\n")
}

# 画分の genome 名を解決する (見出しに含めるため)。
fraction_genome <- function(out_dir, frac_key, desc) {
  sheet <- file.path(out_dir, sprintf("sample_sheet_%s.tsv", frac_key))
  if (file.exists(sheet)) {
    df <- tryCatch(suppressMessages(readr::read_tsv(sheet, show_col_types = FALSE, progress = FALSE)), error = function(e) NULL)
    if (!is.null(df) && "genome" %in% colnames(df) && nrow(df) > 0) {
      g <- as.character(df$genome[[1]])
      if (!is.na(g) && nzchar(g)) return(g)
    }
  }
  b <- basename(desc$deseq2_dir)
  if (nzchar(b) && b != "deseq2") return(b)
  ""
}

# 1 画分 (graft/host) を、通常レポートと同じ章立て (Mapping / PCA / DEG /
# Enrichment) の複数 <section> として組み立てる。
# 戻り値: list(html, nav, n) — n は最後に使った章番号。
build_fraction_sections <- function(reg, frac_key, out_dir, frac_dir, tmp_dir, start_num) {
  base_label <- if (frac_key == "graft") "Graft" else "Host"
  projs <- pic_report_discover_projects(frac_dir)
  if (length(projs) == 0) {
    sid <- paste0(frac_key, "_na")
    return(list(
      html = sprintf('<section id="%s"><h2>%s</h2><p>No DESeq2 output.</p></section>', sid, base_label),
      nav = sprintf('<span class="pic-nav-group">%s</span>', base_label),
      n = start_num
    ))
  }
  desc <- projs[[1]]
  project <- desc$project
  fdr <- pic_plot_spec()$defaults$fdr
  id_prefix <- paste0(frac_key, "__")

  genome <- fraction_genome(out_dir, frac_key, desc)
  label <- if (nzchar(genome)) sprintf("%s (%s)", base_label, genome) else base_label

  stats <- suppressMessages(readr::read_csv(desc$stats_csv, show_col_types = FALSE, progress = FALSE))
  stats <- as.data.frame(stats, check.names = FALSE)
  deg_counts <- pic_load_deg_counts(desc$deseq2_dir, project)
  if (is.null(deg_counts)) {
    deg_counts <- list()
    for (a in unique(stats$aspect)) deg_counts[[a]] <- sum(stats$aspect == a & !is.na(stats$padj) & stats$padj < fdr)
  }

  # 画分の group 配色は、その画分の mapping_sum から取得
  msum <- NULL; group_map <- NULL; group_pal <- NULL
  msum_files <- list.files(frac_dir, pattern = "^mapping_sum__.*\\.tsv$", full.names = TRUE)
  if (length(msum_files) > 0) {
    msum <- read_mapping_sum(msum_files[[1]])
    if (all(c("sample", "group") %in% colnames(msum))) {
      group_map <- stats::setNames(as.character(msum$group), as.character(msum$sample))
      group_pal <- group_palette(as.character(msum$group))
    }
  }

  secs <- character(0)
  navs <- c(sprintf('<span class="pic-nav-group">%s</span>', html_escape(label)))
  n <- start_num

  # Mapping QC
  if (!is.null(msum)) {
    n <- n + 1L
    sid <- paste0(frac_key, "_qc")
    secs <- c(secs, section_mapping_qc(msum, sid, sprintf("%d. %s · Mapping QC", n, label)))
    navs <- c(navs, sprintf('<a href="#%s">Mapping</a>', sid))
  }
  # Aggregation (TSS-TES)
  agg_html <- build_aggregate_html(reg, frac_dir, id_prefix)
  if (nzchar(agg_html)) {
    n <- n + 1L
    sid <- paste0(frac_key, "_agg")
    secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · Aggregation (TSS-TES)</h2>%s</section>', sid, n, label, agg_html))
    navs <- c(navs, sprintf('<a href="#%s">Aggregation</a>', sid))
  }
  # PCA
  n <- n + 1L
  sid <- paste0(frac_key, "_pca")
  pca_html <- build_pca_plots(reg, desc$deseq2_dir, project, group_pal, id_prefix)
  secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · PCA</h2>%s</section>',
                          sid, n, label, if (!is.null(pca_html)) pca_html else "<p>No PCA data.</p>"))
  navs <- c(navs, sprintf('<a href="#%s">PCA</a>', sid))
  # DEG (heatmap + MA/volcano)
  n <- n + 1L
  sid <- paste0(frac_key, "_deg")
  hm_html <- build_heatmap_html(desc$deseq2_dir, project, group_map, group_pal)
  contrast_html <- build_contrast_plots(reg, stats, deg_counts, fdr, id_prefix)
  secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · DEG (DESeq2; FDR = %s)</h2>%s%s</section>',
                          sid, n, label, format(fdr, trim = TRUE),
                          if (!is.null(hm_html)) hm_html else "", contrast_html))
  navs <- c(navs, sprintf('<a href="#%s">DESeq2</a>', sid))
  # Enrichment
  enrich_inner <- enrichment_blocks(desc$enrich_dir, project, tmp_dir, deg_counts, desc$deseq2_dir, group_pal)
  if (nzchar(enrich_inner)) {
    n <- n + 1L
    sid <- paste0(frac_key, "_enrich")
    secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · Enrichment</h2>%s</section>', sid, n, label, enrich_inner))
    navs <- c(navs, sprintf('<a href="#%s">Enrichment</a>', sid))
  }

  list(html = paste(secs, collapse = ""), nav = paste(navs, collapse = ""), n = n)
}

# out_dir 配下に xenograft 分類結果 (classify summary + graft/host) があるか
pic_is_xenograft_out <- function(out_dir) {
  s <- list.files(out_dir, pattern = "^xenograft_classify_summary__.*\\.tsv$", full.names = TRUE)
  length(s) > 0 && (dir.exists(file.path(out_dir, "graft")) || dir.exists(file.path(out_dir, "host")))
}

build_xenograft_report <- function(out_dir, asset_dir) {
  options(pic.report.asset_dir = asset_dir)
  reg <- pic_report_registry()

  summary_files <- list.files(out_dir, pattern = "^xenograft_classify_summary__.*\\.tsv$", full.names = TRUE)
  run <- sub("^xenograft_classify_summary__", "", basename(summary_files[[1]]))
  run <- sub("\\.tsv$", "", run)
  summary_df <- suppressMessages(readr::read_tsv(summary_files[[1]], show_col_types = FALSE, progress = FALSE))
  summary_df <- as.data.frame(summary_df, check.names = FALSE)

  tmp_dir <- file.path(tempdir(), paste0("picxenoreport_", run))

  sec_qc <- section_xenograft_qc(summary_df)
  nav_items <- c('<a href="#classification">Classification</a>')
  body <- sec_qc
  n <- 1L
  for (frac in c("graft", "host")) {
    if (!dir.exists(file.path(out_dir, frac))) next
    tf <- file.path(tmp_dir, frac); dir.create(tf, recursive = TRUE, showWarnings = FALSE)
    res <- build_fraction_sections(reg, frac, out_dir, file.path(out_dir, frac), tf, n)
    body <- paste0(body, res$html)
    nav_items <- c(nav_items, res$nav)
    n <- res$n
  }
  nav <- sprintf('<nav class="pic-nav">%s</nav>', paste(nav_items, collapse = ""))

  out_html <- file.path(out_dir, sprintf("report_%s.html", run))
  render_report_page(paste0(run, " (xenograft)"), nav, body, reg, asset_dir, out_html)
  unlink(tmp_dir, recursive = TRUE)
  out_html
}

pic_load_deg_counts <- function(deseq2_dir, project) {
  f <- file.path(deseq2_dir, "DEG", sprintf("DEG_count_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  if (!all(c("contrast", "deg_count") %in% colnames(d))) return(NULL)
  agg <- stats::aggregate(deg_count ~ contrast, data = d, FUN = sum)
  out <- as.list(stats::setNames(agg$deg_count, agg$contrast))
  out
}

report_runtime_js <- function() {
  paste(readLines(file.path(pic_report_asset_dir(), "report.js"), warn = FALSE), collapse = "\n")
}
