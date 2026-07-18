# 役割: プロジェクト探索 + Mapping QC + Correlation / PCA / Heatmap セクション。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。

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

# enrich ディレクトリ (GSEA/<method>/*.csv を含み project 名に一致するもの) を探す
pic_report_find_enrich_dir <- function(out_dir, project) {
  all_csvs <- list.files(out_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  gsea_csvs <- all_csvs[
    grepl("enrich", all_csvs, fixed = TRUE) &
      grepl(file.path("GSEA", ""), all_csvs, fixed = TRUE) &
      grepl(paste0(project, ".csv"), basename(all_csvs), fixed = TRUE)
  ]
  if (length(gsea_csvs) == 0) {
    return(NA_character_)
  }
  # enrich/GSEA/<METHOD>/file -> enrich は 3 階層上 (xenograft の <frac>/enrich/... も同様)
  enrich_root <- dirname(dirname(dirname(gsea_csvs[[1]])))
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
  trimmed = "Trimmed", unmapped = "Unmapped", multimapping = "Multi-mapping",
  nofeatures = "No feature", ambiguity = "Ambiguous", assigned = "Assigned"
)

read_mapping_sum <- function(path) {
  df <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  df <- as.data.frame(df, check.names = FALSE)
  df
}

section_mapping_qc <- function(msum, section_id = "qc", heading = "1. Mapping QC", src = "", group_map = NULL, group_pal = NULL) {
  fate_cols <- names(PIC_FATE_COLORS)
  have_fate <- all(fate_cols %in% colnames(msum))
  rd_cap <- paste0(section_id, "_rd_cap")
  sd_cap <- paste0(section_id, "_sd_cap")
  parts <- c(sprintf('<section id="%s"><h2>%s</h2>', section_id, heading))

  # サンプル名を group 配色で (heatmap と同じ)
  scol <- function(s) sample_color(s, group_map, group_pal)

  # ---- 100% 積み上げ棒 ----
  if (have_fate) {
    parts <- c(parts,
               sub_head("Read distribution", src_note(src), png_button(rd_cap, "read_distribution")),
               paste0('<p class="pic-note">Each trimmed read is classified into one category: uniquely <b>Assigned</b> to a gene, ',
                      'or lost as <b>No feature</b> (overlaps no gene), <b>Ambiguous</b> (overlaps several genes), ',
                      '<b>Multi-mapping</b>, or <b>Unmapped</b>. A too-low Assigned fraction &mdash; a large ',
                      '<b>No feature</b> / <b>Ambiguous</b> / <b>Unmapped</b> share &mdash; flags a poor-quality sample. ',
                      'Bars are scaled to 100% so samples are comparable regardless of depth.</p>'),
               sprintf('<div class="pic-capbox" id="%s">', rd_cap))
    # 目盛り (0-100%)
    ruler <- '<div class="pic-bar-ruler"><span>0%</span><span>20%</span><span>40%</span><span>60%</span><span>80%</span><span>100%</span></div>'
    rows <- character(0)
    grp_prev <- NULL
    for (i in seq_len(nrow(msum))) {
      sample <- as.character(msum$sample[[i]])
      g_cur <- if ("group" %in% colnames(msum)) as.character(msum$group[[i]]) else NA
      sep_cls <- if (!is.null(grp_prev) && !identical(g_cur, grp_prev)) " pic-grp-start" else ""
      grp_prev <- g_cur
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
        '<div class="pic-bar-row%s" data-group="%s"><div class="pic-bar-label" style="color:%s;font-weight:600">%s</div><div class="pic-bar-track">%s</div></div>',
        sep_cls, html_escape(if (is.na(g_cur)) sample else g_cur), scol(sample), html_escape(sample), paste(segs, collapse = "")
      ))
    }
    legend_items <- vapply(fate_cols, function(cc) sprintf(
      '<span class="pic-legend-item"><span class="pic-swatch" style="background:%s"></span>%s</span>',
      PIC_FATE_COLORS[[cc]], PIC_FATE_LABELS[[cc]]
    ), character(1))
    parts <- c(parts,
      sprintf('<div class="pic-legend" style="margin-left:0;margin-bottom:6px">%s</div>', paste(legend_items, collapse = "")),
      '<div class="pic-bars">',
      ruler,
      paste(rows, collapse = "\n"),
      '</div></div>'   # close pic-bars, pic-capbox(rd)
    )
  }

  # ---- データバー表 ----
  parts <- c(parts,
             sub_head("Sequencing depth", src_note(src), png_button(sd_cap, "sequencing_depth")),
             paste0('<p class="pic-note">A <b>UMI</b> (unique molecular identifier) is a random barcode attached to each cDNA ',
                    'during reverse transcription, so PCR duplicates can be collapsed and each original transcript counted once. ',
                    'How to read each column:</p>',
                    '<p class="pic-note"><b>Total reads</b> &mdash; raw sequenced reads for the sample (before de-duplication). ',
                    '<b>UMIs</b> &mdash; distinct transcript molecules after collapsing duplicates; this is the effective depth ',
                    'and the main measure of sensitivity. <b>Genes</b> &mdash; genes detected with at least one UMI (transcriptome ',
                    'breadth). <b>UMIs/gene</b> &mdash; average molecules per detected gene (sampling depth per gene). ',
                    '<b>Assigned/UMI</b> &mdash; reads per UMI, i.e. sequencing saturation: values near 1 mean little redundancy ',
                    '(sequence deeper to gain UMIs), high values mean molecules are already read many times.</p>'),
             sprintf('<div class="pic-capbox" id="%s">', sd_cap))
  bar_cols <- list(
    total = list(label = "Total reads", grad = c("#9DC3E6", "#D9E7F5")),
    umis  = list(label = "UMIs",  grad = c("#63C384", "#D6EFDD")),
    genes = list(label = "Genes", grad = c("#FFC000", "#FFE9AE"))
  )
  plain_cols <- list()
  if ("umis/genes" %in% colnames(msum)) plain_cols[["umis/genes"]] <- "UMIs/gene"
  if ("assigned/umis" %in% colnames(msum)) plain_cols[["assigned/umis"]] <- "Assigned/UMI"

  maxima <- lapply(names(bar_cols), function(cc) {
    if (cc %in% colnames(msum)) max(suppressWarnings(as.numeric(msum[[cc]])), na.rm = TRUE) else NA_real_
  })
  names(maxima) <- names(bar_cols)

  thead <- '<tr><th>Sample</th>'
  for (cc in names(bar_cols)) if (cc %in% colnames(msum)) thead <- paste0(thead, sprintf("<th>%s</th>", bar_cols[[cc]]$label))
  for (cc in names(plain_cols)) thead <- paste0(thead, sprintf("<th>%s</th>", plain_cols[[cc]]))
  thead <- paste0(thead, "</tr>")

  trows <- character(0)
  grp_prev_t <- NULL
  for (i in seq_len(nrow(msum))) {
    g_cur <- if ("group" %in% colnames(msum)) as.character(msum$group[[i]]) else NA
    sep_cls <- if (!is.null(grp_prev_t) && !identical(g_cur, grp_prev_t)) " pic-grp-start" else ""
    grp_prev_t <- g_cur
    cells <- sprintf('<td class="pic-td-sample" style="color:%s">%s</td>',
                     scol(as.character(msum$sample[[i]])), html_escape(as.character(msum$sample[[i]])))
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
    g_row <- if ("group" %in% colnames(msum)) as.character(msum$group[[i]]) else as.character(msum$sample[[i]])
    trows <- c(trows, sprintf('<tr class="%s" data-group="%s">%s</tr>', trimws(sep_cls), html_escape(g_row), cells))
  }

  parts <- c(parts,
    '<table class="pic-table pic-databar-table"><thead>',
    thead,
    '</thead><tbody>',
    paste(trows, collapse = "\n"),
    '</tbody></table>'
  )

  parts <- c(parts, '</div></section>')
  paste(parts, collapse = "\n")
}

# ---------------------------------------------------------------------------
# DESeq2: PCA / heatmap (plotly)
# ---------------------------------------------------------------------------

pca_variance_pct <- function(deseq2_dir, project) {
  f <- file.path(pca_dir_of(deseq2_dir, project), sprintf("PCA_Variance_%s.csv", project))
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

# sample_sheet(_<frac>).tsv を読み、レポート全体で使う sample / group の正順を返す。
# 見つからなければ NULL (この場合は従来どおり mapping_sum の順に従う)。
pic_read_sample_sheet <- function(out_dir, frac_key = NULL) {
  cand <- if (!is.null(frac_key)) file.path(out_dir, sprintf("sample_sheet_%s.tsv", frac_key))
          else file.path(out_dir, "sample_sheet.tsv")
  if (!file.exists(cand)) return(NULL)
  df <- tryCatch(suppressMessages(readr::read_tsv(cand, show_col_types = FALSE, progress = FALSE)),
                 error = function(e) NULL)
  if (is.null(df) || !("sample" %in% colnames(df))) return(NULL)
  df <- as.data.frame(df, check.names = FALSE)
  list(samples = as.character(df$sample),
       groups = if ("group" %in% colnames(df)) unique(as.character(df$group)) else character(0),
       genome = if ("genome" %in% colnames(df) && nrow(df) > 0) as.character(df$genome[[1]]) else NA_character_)
}

# base_dir の mapping_sum__<run>.tsv から run 名を取り出す (無ければ "")。
pic_report_run_name <- function(base_dir) {
  mf <- list.files(base_dir, pattern = "^mapping_sum__.*\\.tsv$", full.names = FALSE)
  if (length(mf) == 0) return("")
  sub("\\.tsv$", "", sub("^mapping_sum__", "", mf[[1]]))
}

# project (= <run>_<genome>) と base_dir から genome を解決する。
# 優先: sample_sheet の genome 列 → 次善: project から run 接頭辞を除去。
pic_report_genome_of <- function(base_dir, project, frac_key = NULL) {
  sheet <- pic_read_sample_sheet(base_dir, frac_key)
  if (!is.null(sheet) && !is.na(sheet$genome) && nzchar(sheet$genome)) return(sheet$genome)
  run <- pic_report_run_name(base_dir)
  if (nzchar(run) && startsWith(project, paste0(run, "_"))) return(substring(project, nchar(run) + 2L))
  project
}

# data.frame を col の値が order に一致する順に並べ替える。order にない値は末尾に元順で残す。
pic_reorder_rows <- function(x, col, order) {
  if (is.null(order) || length(order) == 0 || !(col %in% colnames(x))) return(x)
  v <- as.character(x[[col]])
  rank <- match(v, order)
  na <- is.na(rank)
  rank[na] <- length(order) + seq_len(sum(na))
  x[order(rank), , drop = FALSE]
}

# 文字列ベクトルを order の順に並べ替える (order にないものは末尾)。
pic_reorder_vec <- function(v, order) {
  if (is.null(order) || length(order) == 0) return(v)
  rank <- match(v, order)
  na <- is.na(rank); rank[na] <- length(order) + seq_len(sum(na))
  v[order(rank)]
}

# contrast (aspect "num / den") を group の正順に並べ替える (大文字小文字は無視)。
pic_order_contrasts <- function(aspects, group_order) {
  if (is.null(group_order) || length(group_order) == 0) return(aspects)
  go <- tolower(group_order)
  keyf <- function(a) {
    sp <- strsplit(a, " / ", fixed = TRUE)[[1]]
    num <- if (length(sp) >= 1) tolower(sp[[1]]) else ""
    den <- if (length(sp) >= 2) tolower(sp[[2]]) else ""
    ni <- match(num, go); di <- match(den, go)
    if (is.na(ni)) ni <- length(go) + 1L
    if (is.na(di)) di <- length(go) + 1L
    ni * 1000 + di
  }
  ks <- vapply(aspects, keyf, numeric(1))
  aspects[order(ks)]
}

# 相関値 (lo..1) を青-白-赤の hex に写像する。
cor_color <- function(v, lo) {
  ramp <- grDevices::colorRamp(c("#2166ac", "#f7f7f7", "#b2182b"))
  t <- (as.numeric(v) - lo) / max(1e-9, (1 - lo))
  t[!is.finite(t)] <- 0.5; t <- pmax(0, pmin(1, t))
  mm <- ramp(t)
  grDevices::rgb(mm[, 1], mm[, 2], mm[, 3], maxColorValue = 255)
}
# 背景 hex の輝度から文字色 (白 or 濃) を決める。
contrast_text <- function(hex) {
  h <- sub("^#", "", hex)
  r <- strtoi(substr(h, 1, 2), 16L); g <- strtoi(substr(h, 3, 4), 16L); b <- strtoi(substr(h, 5, 6), 16L)
  lum <- 0.299 * r + 0.587 * g + 0.114 * b
  ifelse(lum < 140, "#ffffff", "#1f2933")
}

# サンプル間相関ヒートマップ (deseq2/correlation_<project>.csv を HTML テーブルで図示)。
# 赤青の発散配色、セルに相関値 (小数 2 桁)、ラベルは group 配色 (列は上部・横書き)、
# group 境界に区切り線 (ラベルも含め、テーブル内に収まりはみ出さない)。
build_correlation_html <- function(reg, deseq2_dir, project, group_map = NULL, group_pal = NULL,
                                   sample_order = NULL, id_prefix = "", proj_dir = NULL) {
  f <- file.path(deseq2_dir, sprintf("correlation_%s.csv", project))
  if (!file.exists(f)) return("")
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  d <- as.data.frame(d, check.names = FALSE)
  if (ncol(d) < 3) return("")
  rn <- as.character(d[[1]]); d <- d[, -1, drop = FALSE]
  cn <- colnames(d)
  samples <- pic_reorder_vec(intersect(cn, rn), sample_order)
  if (length(samples) < 2) return("")
  m <- as.matrix(d)[match(samples, rn), match(samples, cn), drop = FALSE]
  storage.mode(m) <- "double"
  n <- length(samples)

  grpof <- function(s) if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s
  gcol <- function(s) { g <- grpof(s); if (!is.null(group_pal) && g %in% names(group_pal)) unname(group_pal[[g]]) else "#1f2933" }
  grp <- vapply(samples, grpof, character(1))
  gstart <- c(FALSE, grp[-1] != grp[-length(grp)])   # 各サンプルが group の先頭か
  zmin <- suppressWarnings(min(m[is.finite(m)])); if (!is.finite(zmin)) zmin <- 0

  # ヘッダ (列ラベル: 上部・縦書き・group 配色)。data-gcol で group 列を識別。
  ths <- vapply(seq_len(n), function(j)
    sprintf('<th class="pic-cor-ch%s" data-gcol="%s"><span style="color:%s">%s</span></th>',
            if (gstart[[j]]) " gsep-l" else "", html_escape(grpof(samples[[j]])), gcol(samples[[j]]), html_escape(samples[[j]])),
    character(1))
  header <- sprintf('<tr><th class="pic-cor-corner"></th>%s</tr>', paste(ths, collapse = ""))

  rows <- vapply(seq_len(n), function(i) {
    rowsep <- if (gstart[[i]]) " gsep-t" else ""
    rh <- sprintf('<th class="pic-cor-rh%s"><span style="color:%s">%s</span></th>',
                  rowsep, gcol(samples[[i]]), html_escape(samples[[i]]))
    tds <- vapply(seq_len(n), function(j) {
      v <- m[i, j]
      bg <- if (is.finite(v)) cor_color(v, zmin) else "#ffffff"
      fc <- contrast_text(bg)
      cls <- paste0("val", if (gstart[[j]]) " gsep-l" else "", if (gstart[[i]]) " gsep-t" else "")
      sprintf('<td class="%s" data-gcol="%s" style="background:%s;color:%s">%s</td>', cls, html_escape(grpof(samples[[j]])), bg, fc,
              if (is.finite(v)) sprintf("%.2f", v) else "")
    }, character(1))
    sprintf('<tr data-grow="%s">%s%s</tr>', html_escape(grpof(samples[[i]])), rh, paste(tds, collapse = ""))
  }, character(1))

  cap_id <- paste0(id_prefix, "cor_cap")
  paste0(
    src_note(report_rel_path(f, proj_dir)),
    sprintf('<span class="pic-headact">%s</span>', png_button(cap_id, "sample_correlation")),
    '<p class="pic-note">Sample-to-sample correlation. Replicates of the same group should correlate most strongly ',
    '(<b style="color:#b2182b">red</b>); a sample that stands out from its group (more <b style="color:#2166ac">blue</b>) may be an outlier.</p>',
    sprintf('<div class="pic-capbox" id="%s"><div class="pic-cor-wrap"><table class="pic-cor"><thead>%s</thead><tbody>%s</tbody></table></div></div>',
            cap_id, header, paste(rows, collapse = "")))
}

build_pca_plots <- function(reg, deseq2_dir, project, group_pal = NULL, id_prefix = "", proj_dir = NULL) {
  f <- file.path(pca_dir_of(deseq2_dir, project), sprintf("PCA_RegLog_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  d <- as.data.frame(d, check.names = FALSE)
  if (!("group" %in% colnames(d))) d$group <- "all"
  varpct <- pca_variance_pct(deseq2_dir, project)
  pal <- ensure_palette(group_pal, d$group)
  # group (凡例) の順は palette (= sample_sheet) の順に従う
  groups_ord <- pic_reorder_vec(unique(as.character(d$group)), names(pal))

  # ---- 寄与率 (scree) ----
  if (!is.null(varpct)) {
    pcs <- names(varpct)
    npc <- min(length(pcs), 15L)
    pcs <- pcs[seq_len(npc)]
    scree <- list(
      list(x = as.list(pcs), y = as.list(unname(varpct[pcs])), type = "bar",
           name = "variance", marker = list(color = "#5B9BD5"),
           hovertemplate = "%{x}: %{y:.1f}%<extra></extra>")
    )
    register_plot(reg, paste0(id_prefix, "pca_scree"), list(data = scree, layout = list(
      xaxis = list(title = ""),
      yaxis = list(title = "% variance", rangemode = "tozero"),  # 上限は自動 (100% までは不要)
      hovermode = "x", margin = list(l = 56, r = 16, t = 8, b = 28), showlegend = FALSE  # 1 色なので凡例不要
    )))
  }
  scree_block <- sprintf('<div class="pic-plot-cell"><div id="%spca_scree" class="pic-plot"></div></div>', id_prefix)

  pair_blocks <- character(0)
  pairs <- list(c("PC1", "PC2"), c("PC2", "PC3"))
  for (pr in pairs) {
    xc <- pr[[1]]; yc <- pr[[2]]
    if (!all(c(xc, yc) %in% colnames(d))) next
    traces <- list()
    for (g in groups_ord) {
      sub <- d[d$group == g, , drop = FALSE]
      traces[[length(traces) + 1]] <- list(
        x = as.list(as.numeric(sub[[xc]])),
        y = as.list(as.numeric(sub[[yc]])),
        text = as.list(as.character(sub$sample)),
        name = group_span(g, pal),
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
    pair_blocks <- c(pair_blocks, sprintf('<div class="pic-plot-cell"><div id="%s" class="pic-plot"></div></div>', id))
  }
  if (is.null(varpct) && length(pair_blocks) == 0) return(NULL)
  # 1 行目: 寄与率のみ / 2 行目: PC1-2, PC2-3
  paste0(
    src_note(report_rel_path(f, proj_dir)),
    '<p class="pic-note">PCA summarizes each sample&rsquo;s overall expression into a few axes so you can see how samples relate. ',
    'Samples close together are similar; replicates of the same group should cluster. The scree plot shows how much variance each axis captures.</p>',
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
build_heatmap_html <- function(deseq2_dir, project, group_map = NULL, group_pal = NULL, proj_dir = NULL, sample_order = NULL, id_prefix = "") {
  f <- file.path(deg_dir_of(deseq2_dir, project), sprintf("DEG_normalizedCountTable_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  d <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  d <- as.data.frame(d, check.names = FALSE)
  ann <- c("ens_gene", "ext_gene", "biotype", "chr")
  sample_cols <- setdiff(colnames(d), ann)
  if (length(sample_cols) < 2 || nrow(d) < 2) return(NULL)
  # 列 (サンプル) を sample_sheet の順に並べる
  sample_cols <- pic_reorder_vec(sample_cols, sample_order)

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

  # group 境界 (列がグループの先頭なら区切り線)
  gvec <- vapply(scols, function(s) if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s, character(1))
  gstart <- c(FALSE, gvec[-1] != gvec[-length(gvec)])

  # ヘッダ (サンプル名・グループ色, sticky)
  ths <- vapply(seq_along(scols), function(j) {
    s <- scols[[j]]
    g <- if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s
    col <- if (!is.null(group_pal) && g %in% names(group_pal)) unname(group_pal[[g]]) else "#333333"
    sprintf('<th class="s%s" data-gcol="%s"><span style="color:%s">%s</span></th>',
            if (gstart[[j]]) " gsep" else "", html_escape(g), col, html_escape(s))
  }, character(1))
  header <- sprintf('<tr><th class="corner"></th>%s</tr>', paste(ths, collapse = ""))

  rows <- vapply(seq_len(nrow(z)), function(i) {
    tds <- vapply(seq_len(ncol(z)), function(j) sprintf(
      '<td class="%s" data-gcol="%s" style="background:%s" title="Gene: %s&#10;Sample: %s&#10;Z: %.2f"></td>',
      if (gstart[[j]]) "gsep" else "", html_escape(gvec[[j]]), colmat[i, j], html_escape(labels[[i]]), html_escape(scols[[j]]), z[i, j]
    ), character(1))
    sprintf('<tr><th class="g">%s</th>%s</tr>', html_escape(labels[[i]]), paste(tds, collapse = ""))
  }, character(1))

  cap_id <- paste0(id_prefix, "heatmap_cap")
  # 左 = group 表示切替 / 右 = ヒートマップ本体 (タイトル右に Download csv + PNG ボタン)
  ctrl <- group_toggle_panel(group_pal)
  view <- paste0(
    sub_head(sprintf('Expression heatmap (%d DEGs)', nrow(z)),
             src_note(report_rel_path(f, proj_dir)), png_button(cap_id, "deg_heatmap")),
    '<p class="pic-note">Each row is a gene, each column a sample (grouped by color). ',
    'Cell color is the gene&rsquo;s z-score across samples &mdash; <b style="color:#b2182b">red</b> is higher than that gene&rsquo;s average, ',
    '<b style="color:#2166ac">blue</b> is lower. Genes are ordered by similarity. Hover a cell for the gene, sample, and z-score.</p>',
    sprintf('<div class="pic-capbox" id="%s">', cap_id),
    sprintf('<div class="pic-hm-scroll"><table class="pic-hm"><thead>%s</thead><tbody>%s</tbody></table></div>',
            header, paste(rows, collapse = "")),
    '</div>'
  )
  # 左 (control) の先頭に <!--SUBNAV--> を置き、pic_tab_panel がサブタブ選択を差し込む。
  sprintf('<div class="pic-hmsel"><aside class="pic-hmsel-ctrl"><!--SUBNAV-->%s</aside><div class="pic-hmsel-view">%s</div></div>',
          ctrl, view)
}

