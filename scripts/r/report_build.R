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

# project フォルダ (report 出力先) から見た相対パス。
report_rel_path <- function(path, base) {
  if (is.null(base) || !nzchar(base) || is.null(path) || !nzchar(path)) return("")
  p <- tryCatch(normalizePath(path, winslash = "/", mustWork = FALSE), error = function(e) path)
  b <- tryCatch(normalizePath(base, winslash = "/", mustWork = FALSE), error = function(e) base)
  b2 <- paste0(b, "/")
  if (startsWith(p, b2)) substring(p, nchar(b2) + 1) else basename(p)
}

# 「ここで選択」を明示する目立つ案内バッジ。
pick_hint <- function(text) {
  sprintf('<div class="pic-pick"><span class="pic-pick-ic">&#128071;</span><span>%s</span></div>', text)
}

# HTML 要素 (id=cap_id) を PNG 画像でダウンロードするボタン。
png_button <- function(cap_id, name) {
  sprintf('<button class="pic-png-btn" type="button" data-cap="%s" data-name="%s">Download plot as a png.</button>',
          html_escape(cap_id), html_escape(name))
}

# 元データの場所を示すキャプション (グレー背景 + 黄色文字 + コピーボタン)。
src_note <- function(rel) {
  if (is.null(rel) || !nzchar(rel)) return("")
  e <- html_escape(rel)
  sprintf(paste0('<p class="pic-src"><span class="pic-src-label">source</span>',
                 '<code class="pic-src-path">%s</code>',
                 '<button class="pic-src-copy" type="button" data-copy="%s">copy path</button></p>'),
          e, e)
}

# 複数ブロックをチェックボックスで表示選択する UI。既定は先頭のみ表示。
# items: list(list(id=, label=, html=, checked=logical), ...)
# ブロックが 1 つだけなら選択 UI は付けずにそのまま表示する。
build_select_group <- function(items) {
  items <- Filter(function(it) !is.null(it$html) && nzchar(it$html), items)
  if (length(items) == 0) return("")
  if (length(items) == 1) {
    return(sprintf('<div class="pic-select-item">%s</div>', items[[1]]$html))
  }
  bar <- vapply(items, function(it) {
    ck <- if (isTRUE(it$checked)) " checked" else ""
    sprintf('<label class="pic-select-chk"><input type="checkbox" data-target="%s"%s>%s</label>',
            html_escape(it$id), ck, html_escape(it$label))
  }, character(1))
  blk <- vapply(items, function(it) {
    hid <- if (isTRUE(it$checked)) "" else " hidden"
    sprintf('<div class="pic-select-item" id="%s"%s>%s</div>', html_escape(it$id), hid, it$html)
  }, character(1))
  paste0('<div class="pic-select-bar">', paste(bar, collapse = ""), '</div>',
         paste(blk, collapse = ""))
}

# contrast を group×group の行列で選ばせる UI。
# 各セル (group i × group j) に、その 2 群の比較があれば DEG 数 + チェックボックス。
# entries: list(list(aspect=, count=, checked=, attr=))
#   attr は checkbox に付与する属性 (例: 'data-target="id"' や 'data-axis="r" data-key="k"')。
# groups: 表示順の group ベクトル (contrast は大文字小文字を無視して照合)。
build_group_matrix <- function(entries, groups, group_pal = NULL, show_count = TRUE) {
  if (length(groups) < 2 || length(entries) == 0) return("")
  gl <- tolower(groups)
  cellmap <- list()
  for (e in entries) {
    sp <- strsplit(e$aspect, " / ", fixed = TRUE)[[1]]
    a <- tolower(trimws(sp[[1]])); b <- if (length(sp) >= 2) tolower(trimws(sp[[2]])) else ""
    ia <- match(a, gl); ib <- match(b, gl)
    if (is.na(ia) || is.na(ib)) next
    cellmap[[paste0(min(ia, ib), "_", max(ia, ib))]] <- e
  }
  if (length(cellmap) == 0) return("")
  hcol <- function(g) if (!is.null(group_pal) && g %in% names(group_pal)) unname(group_pal[[g]]) else "#1f2933"
  n <- length(groups)
  ths <- paste(vapply(seq_len(n), function(j)
    sprintf('<th class="pic-cmx-h"><span style="color:%s">%s</span></th>', hcol(groups[[j]]), html_escape(groups[[j]])),
    character(1)), collapse = "")
  header <- sprintf('<tr><th class="pic-cmx-corner"></th>%s</tr>', ths)
  rows_html <- character(0)
  for (i in seq_len(n)) {
    rh <- sprintf('<th class="pic-cmx-rh"><span style="color:%s">%s</span></th>', hcol(groups[[i]]), html_escape(groups[[i]]))
    cs <- character(0)
    for (j in seq_len(n)) {
      # 各比較は上三角 (row < col) の 1 セルにのみ表示し、重複を避ける
      e <- if (i < j) cellmap[[paste0(i, "_", j)]] else NULL
      if (is.null(e)) { cs <- c(cs, '<td class="pic-cmx-empty"></td>'); next }
      ck <- if (isTRUE(e$checked)) " checked" else ""
      cnt <- if (show_count && !is.null(e$count) && is.finite(suppressWarnings(as.numeric(e$count))))
               sprintf('<span class="pic-cmx-n">%s</span>', fmt_int(e$count)) else ""
      cs <- c(cs, sprintf('<td class="pic-cmx-cell" title="%s"><label>%s<input type="checkbox" %s%s></label></td>',
                          html_escape(e$aspect), cnt, e$attr, ck))
    }
    rows_html <- c(rows_html, sprintf('<tr>%s%s</tr>', rh, paste(cs, collapse = "")))
  }
  cls <- if (show_count) "pic-cmatrix" else "pic-cmatrix pic-cmx-nocount"
  sprintf('<div class="pic-cmx-wrap"><table class="%s"><thead>%s</thead><tbody>%s</tbody></table></div>',
          cls, header, paste(rows_html, collapse = ""))
}

# 2 軸 (行 = contrast, 列 = method) をチェックボックスで選び、両方選択のセルのみ表示。
# row_groups を渡すと、行 (contrast) の選択を group×group 行列で提示する。
build_matrix_group <- function(rows, cols, cell_html, row_title = "Contrast", col_title = "Method",
                               row_groups = NULL, row_counts = NULL, group_pal = NULL) {
  if (length(rows) == 0 || length(cols) == 0) return("")
  first_row <- rows[[1]]$key
  first_col <- cols[[1]]
  # 行 (contrast) セレクタ: group 行列 or 線形バー
  gm <- ""
  if (!is.null(row_groups) && length(row_groups) >= 2) {
    entries <- lapply(seq_along(rows), function(i) list(
      aspect = rows[[i]]$label,
      count = if (!is.null(row_counts) && rows[[i]]$label %in% names(row_counts)) row_counts[[rows[[i]]$label]] else NULL,
      checked = (i == 1L),
      attr = sprintf('data-axis="r" data-key="%s"', html_escape(rows[[i]]$key))))
    gm <- build_group_matrix(entries, row_groups, group_pal, show_count = !is.null(row_counts))
  }
  row_sel <- if (nzchar(gm)) gm else {
    row_bar <- vapply(seq_along(rows), function(i) {
      ck <- if (i == 1L) " checked" else ""
      sprintf('<label class="pic-select-chk"><input type="checkbox" data-axis="r"%s data-key="%s">%s</label>',
              ck, html_escape(rows[[i]]$key), html_escape(rows[[i]]$label))
    }, character(1))
    sprintf('<div class="pic-select-bar"><span class="pic-select-lbl">%s:</span>%s</div>',
            html_escape(row_title), paste(row_bar, collapse = ""))
  }
  col_bar <- vapply(seq_along(cols), function(i) {
    ck <- if (i == 1L) " checked" else ""
    sprintf('<label class="pic-select-chk"><input type="checkbox" data-axis="c"%s data-key="%s">%s</label>',
            ck, html_escape(cols[[i]]), html_escape(cols[[i]]))
  }, character(1))
  cells <- character(0)
  for (r in rows) {
    for (co in cols) {
      h <- cell_html(r$key, co)
      if (is.null(h) || !nzchar(h)) next
      vis <- (r$key == first_row && co == first_col)
      hid <- if (vis) "" else " hidden"
      cells <- c(cells, sprintf('<div class="pic-select-item" data-r="%s" data-c="%s"%s>%s</div>',
                                html_escape(r$key), html_escape(co), hid, h))
    }
  }
  if (length(cells) == 0) return("")
  sprintf('<div class="pic-matrix">%s<div class="pic-select-bar"><span class="pic-select-lbl">%s:</span>%s</div><div class="pic-matrix-cells">%s</div></div>',
          row_sel, html_escape(col_title), paste(col_bar, collapse = ""),
          paste(cells, collapse = ""))
}

# plot 仕様を共有レジストリ env に登録する
pic_report_registry <- function() {
  reg <- new.env(parent = emptyenv())
  reg$plots <- list()
  reg$expr <- list()   # 遺伝子発現データ (id -> normalized counts など)
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
  # enrich/<genome>/GSEA/<METHOD>/file -> enrich/<genome> は 3 階層上
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
  trimmed = "trimmed", unmapped = "unmapped", multimapping = "multimapping",
  nofeatures = "nofeatures", ambiguity = "ambiguity", assigned = "assigned"
)

read_mapping_sum <- function(path) {
  df <- suppressMessages(readr::read_tsv(path, show_col_types = FALSE, progress = FALSE))
  df <- as.data.frame(df, check.names = FALSE)
  df
}

section_mapping_qc <- function(msum, section_id = "qc", heading = "1. Mapping QC", src = "", group_map = NULL, group_pal = NULL) {
  fate_cols <- names(PIC_FATE_COLORS)
  have_fate <- all(fate_cols %in% colnames(msum))
  cap_id <- paste0(section_id, "_cap")
  parts <- c(sprintf('<section id="%s"><h2>%s</h2>', section_id, heading), src_note(src),
             png_button(cap_id, "mapping_qc"),
             sprintf('<div class="pic-capbox" id="%s">', cap_id))

  # サンプル名を group 配色で (heatmap と同じ)
  scol <- function(s) {
    g <- if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s
    if (!is.null(group_pal) && g %in% names(group_pal)) unname(group_pal[[g]]) else "#333333"
  }

  # ---- 100% 積み上げ棒 ----
  if (have_fate) {
    parts <- c(parts, '<h3>Read distribution</h3>')
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
        '<div class="pic-bar-row%s"><div class="pic-bar-label" style="color:%s;font-weight:600">%s</div><div class="pic-bar-track">%s</div></div>',
        sep_cls, scol(sample), html_escape(sample), paste(segs, collapse = "")
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
    trows <- c(trows, sprintf('<tr class="%s">%s</tr>', trimws(sep_cls), cells))
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
       groups = if ("group" %in% colnames(df)) unique(as.character(df$group)) else character(0))
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

  # ヘッダ (列ラベル: 上部・横書き・group 配色)
  ths <- vapply(seq_len(n), function(j)
    sprintf('<th class="pic-cor-ch%s"><span style="color:%s">%s</span></th>',
            if (gstart[[j]]) " gsep-l" else "", gcol(samples[[j]]), html_escape(samples[[j]])),
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
      sprintf('<td class="%s" style="background:%s;color:%s">%s</td>', cls, bg, fc,
              if (is.finite(v)) sprintf("%.2f", v) else "")
    }, character(1))
    sprintf('<tr>%s%s</tr>', rh, paste(tds, collapse = ""))
  }, character(1))

  cap_id <- paste0(id_prefix, "cor_cap")
  paste0(
    src_note(report_rel_path(f, proj_dir)),
    '<p class="pic-note">Sample-to-sample correlation. Replicates of the same group should correlate most strongly ',
    '(<b style="color:#b2182b">red</b>); a sample that stands out from its group (more <b style="color:#2166ac">blue</b>) may be an outlier.</p>',
    png_button(cap_id, "sample_correlation"),
    sprintf('<div class="pic-capbox" id="%s"><div class="pic-cor-wrap"><table class="pic-cor"><thead>%s</thead><tbody>%s</tbody></table></div></div>',
            cap_id, header, paste(rows, collapse = "")))
}

build_pca_plots <- function(reg, deseq2_dir, project, group_pal = NULL, id_prefix = "", proj_dir = NULL) {
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
  # group (凡例) の順は palette (= sample_sheet) の順に従う
  groups_ord <- pic_reorder_vec(unique(as.character(d$group)), names(pal))

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
    for (g in groups_ord) {
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
  f <- file.path(deseq2_dir, "DEG", sprintf("DEG_normalizedCountTable_%s.csv", project))
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
    sprintf('<th class="s%s"><span style="color:%s">%s</span></th>', if (gstart[[j]]) " gsep" else "", col, html_escape(s))
  }, character(1))
  header <- sprintf('<tr><th class="corner"></th>%s</tr>', paste(ths, collapse = ""))

  rows <- vapply(seq_len(nrow(z)), function(i) {
    tds <- vapply(seq_len(ncol(z)), function(j) sprintf(
      '<td class="%s" style="background:%s" title="Gene: %s&#10;Sample: %s&#10;Z: %.2f"></td>',
      if (gstart[[j]]) "gsep" else "", colmat[i, j], html_escape(labels[[i]]), html_escape(scols[[j]]), z[i, j]
    ), character(1))
    sprintf('<tr><th class="g">%s</th>%s</tr>', html_escape(labels[[i]]), paste(tds, collapse = ""))
  }, character(1))

  cap_id <- paste0(id_prefix, "heatmap_cap")
  paste0(
    sprintf('<h3>Expression heatmap (%d differentially expressed genes)</h3>', nrow(z)),
    src_note(report_rel_path(f, proj_dir)),
    '<p class="pic-note">Each row is a gene, each column a sample (grouped by color). ',
    'Cell color is the gene&rsquo;s z-score across samples &mdash; <b style="color:#b2182b">red</b> is higher than that gene&rsquo;s average, ',
    '<b style="color:#2166ac">blue</b> is lower. Genes are ordered by similarity. Hover a cell for the gene, sample, and z-score.</p>',
    png_button(cap_id, "deg_heatmap"),
    sprintf('<div class="pic-capbox" id="%s">', cap_id),
    sprintf('<div class="pic-hm-scroll"><table class="pic-hm"><thead>%s</thead><tbody>%s</tbody></table></div>',
            header, paste(rows, collapse = "")),
    '</div>'
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
    meta_html <- sprintf('<p class="pic-contrast-meta"><b>%s</b></p>', html_escape(a))
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

    item_html <- paste0(meta_html, '<div class="pic-ma-vol">', paste(cells, collapse = ""), '</div>')
    items[[length(items) + 1L]] <- list(
      id = sprintf("%ssel_%s", id_prefix, flab),
      label = a, aspect = a, count = nd,
      html = item_html,
      checked = (k == 1L)
    )
  }

  note <- paste0('<h3>Differential expression (MA &amp; volcano)</h3>',
    src_note(stats_src),
    '<p class="pic-note">Each gene is one point. The <b>MA plot</b> shows the fold change vs. average expression; ',
    'the <b>volcano plot</b> shows the fold change vs. statistical significance (&minus;log<sub>10</sub> p-value). ',
    '<b style="color:#d7301f">Red</b> / <b style="color:#2166ac">blue</b> points are significantly up / down (padj below the FDR); grey is not significant. ',
    'Hover a point for its gene name, fold change and padj.</p>')

  # contrast 選択: group×group 行列 (各セルに DEG 数 + チェックボックス)。1 つだけなら行列不要。
  groups <- if (!is.null(group_order) && length(group_order) > 0) group_order else
    unique(unlist(lapply(aspects, function(a) trimws(strsplit(a, " / ", fixed = TRUE)[[1]]))))
  entries <- lapply(items, function(it) list(
    aspect = it$aspect, count = it$count, checked = it$checked,
    attr = sprintf('data-target="%s"', it$id)))
  matrix_html <- build_group_matrix(entries, groups, group_pal, show_count = TRUE)
  blocks <- vapply(items, function(it) {
    hid <- if (isTRUE(it$checked)) "" else " hidden"
    sprintf('<div class="pic-select-item" id="%s"%s>%s</div>', it$id, hid, it$html)
  }, character(1))
  sel_html <- if (length(items) <= 1) {
    if (length(items) == 1) sprintf('<div class="pic-select-item">%s</div>', items[[1]]$html) else ""
  } else if (nzchar(matrix_html)) {
    paste0(pick_hint('Select comparisons to show &mdash; tick a cell in the matrix below.'),
           '<p class="pic-note">Each cell shows the number of differentially expressed genes ',
           '(padj &lt; FDR) between the two groups; the plots appear underneath.</p>',
           matrix_html, paste(blocks, collapse = ""))
  } else {
    build_select_group(items)
  }
  paste0(note, sel_html)
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
build_aggregate_html <- function(reg, base_dir, id_prefix = "", proj_dir = NULL, sample_order = NULL, group_pal = NULL) {
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
      pal <- if (!is.null(group_pal)) group_pal else group_palette(unique(as.character(d$group)))
      miss <- setdiff(unique(as.character(d$group)), names(pal)); if (length(miss) > 0) pal <- c(pal, group_palette(miss))
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
          name = if (first) g else s, legendgroup = g, showlegend = first,
          mode = "lines", type = "scatter",
          line = list(color = unname(pal[[g]]), width = 1.4),
          hovertemplate = sprintf("%s<br>CPM: %%{y:.2f}<extra></extra>", html_escape(s)))
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
        hovermode = "closest", margin = list(l = 60, r = 20, t = 30, b = 46),
        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0))
      register_plot(reg, id, list(data = traces, layout = layout))
      blocks <- c(blocks, sprintf(
        '<div class="pic-plot-cell"><h4>%s — gene-body profile (TSS→TES)</h4>%s<div id="%s" class="pic-plot"></div></div>',
        html_escape(genome), src_note(report_rel_path(csv, proj_dir)), id))
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
  pal <- if (!is.null(group_pal)) group_pal else group_palette(ug)
  miss <- setdiff(ug, names(pal)); if (length(miss) > 0) pal <- c(pal, group_palette(miss))
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
    'Hover over a gene name to see its padj for every comparison.</p>',
    '%s',
    '<div class="pic-expr" data-expr="%s">',
    '<div class="pic-expr-typeahead"><div class="pic-expr-tagbox">',
    '<input class="pic-expr-input" type="text" autocomplete="off" placeholder="Add a gene…"></div>',
    '<div class="pic-expr-ac"></div></div>',
    '<div class="pic-expr-msg"></div><div class="pic-expr-plots"></div></div></section>'),
    section_id, html_escape(heading), format(fdr, trim = TRUE),
    src_note(report_rel_path(countf, proj_dir)), id)
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
  f <- file.path(deseq2_dir, "DEG", "DEGCluster", sprintf("DEGCluster_profile_%s.csv", project))
  if (!file.exists(f)) return(NULL)
  prof <- suppressMessages(readr::read_csv(f, show_col_types = FALSE, progress = FALSE))
  if (nrow(prof) == 0 || !all(c("rlog_expr", "cluster_id", "group") %in% colnames(prof))) return(NULL)
  prof$value <- suppressWarnings(as.numeric(prof$rlog_expr))
  prof <- prof[is.finite(prof$value), , drop = FALSE]
  if (nrow(prof) == 0) return(NULL)

  summ_f <- file.path(deseq2_dir, "DEG", "DEGCluster", sprintf("DEGCluster_summary_%s.csv", project))
  nmap <- NULL
  if (file.exists(summ_f)) {
    summ <- suppressMessages(readr::read_csv(summ_f, show_col_types = FALSE, progress = FALSE))
    if (all(c("cluster_id", "gene_count") %in% colnames(summ))) nmap <- stats::setNames(summ$gene_count, summ$cluster_id)
  }
  clusters <- unique(as.character(prof$cluster_id))
  clab <- function(cl) if (!is.null(nmap) && cl %in% names(nmap)) sprintf("%s (n=%s)", cl, nmap[[cl]]) else cl
  groups <- unique(as.character(prof$group))
  pal <- if (!is.null(group_pal)) group_pal else group_palette(groups)
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
        type = "box", name = g, legendgroup = g, showlegend = (i == 1),
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

# enrichment セクションの中身 (h3 GSEA / ORA) を返す。<section> ラッパは付けない。
enrichment_blocks <- function(enrich_dir, project, tmp_dir, deg_counts = NULL, deseq2_dir = NULL, group_pal = NULL, proj_dir = NULL, id_prefix = "", reg = NULL, group_order = NULL) {
  if (is.null(enrich_dir) || is.na(enrich_dir) || !dir.exists(enrich_dir)) return("")
  ecfg <- pic_plot_spec()$plot$enrichment
  parts <- character(0)

  method_order <- strsplit(pic_plot_spec()$defaults$enrich_methods_csv, ",", fixed = TRUE)[[1]]
  order_methods <- function(ms) {
    c(method_order[method_order %in% ms], sort(setdiff(ms, method_order)))
  }

  # ---- GSEA: contrast × method の 2 軸をチェックボックスで選択 ----
  gsea_root <- file.path(enrich_dir, "GSEA")
  gsea_html <- ""
  if (dir.exists(gsea_root)) {
    mdirs <- list.dirs(gsea_root, recursive = FALSE, full.names = TRUE)
    cellmap <- list()          # cellmap[[contrast]][[method]] = uri
    methods_seen <- character(0)
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
        spec <- tryCatch(build_gsea_plotly(df, numerator, denominator), error = function(e) NULL)
        if (is.null(spec)) next
        pid <- sprintf("%sgseaplot_%s_%s", id_prefix, method, format_contrast_file_label(contrast))
        register_plot(reg, pid, list(data = spec$data, layout = spec$layout, config = spec$config))
        if (is.null(cellmap[[contrast]])) cellmap[[contrast]] <- list()
        cellmap[[contrast]][[method]] <- list(id = pid, h = enrich_plot_height(spec$n_terms))
        methods_seen <- union(methods_seen, method)
      }
    }
    contrasts <- names(cellmap)
    if (length(contrasts) > 0) {
      # contrast の並び順: sample_sheet の group 順 (なければ DEG 数降順)
      if (!is.null(group_order)) {
        contrasts <- pic_order_contrasts(contrasts, group_order)
      } else if (!is.null(deg_counts)) {
        keyv <- vapply(contrasts, function(a) if (a %in% names(deg_counts)) deg_counts[[a]] else 0, numeric(1))
        contrasts <- contrasts[order(keyv, decreasing = TRUE)]
      } else {
        contrasts <- sort(contrasts)
      }
      methods <- order_methods(methods_seen)
      rows <- lapply(contrasts, function(cc) list(key = format_contrast_file_label(cc), label = cc))
      key2contrast <- stats::setNames(contrasts, vapply(contrasts, format_contrast_file_label, character(1)))
      cell_html <- function(rkey, method) {
        contrast <- key2contrast[[rkey]]
        info <- cellmap[[contrast]][[method]]
        if (is.null(info)) return("")
        sprintf('<div class="pic-enrich-item"><h4>%s &mdash; %s</h4><div id="%s" class="pic-plot" style="height:%dpx"></div></div>',
                html_escape(contrast), html_escape(method), info$id, info$h)
      }
      row_groups <- if (!is.null(group_order) && length(group_order) > 0) group_order else
        unique(unlist(lapply(contrasts, function(a) trimws(strsplit(a, " / ", fixed = TRUE)[[1]]))))
      gsea_html <- build_matrix_group(rows, methods, cell_html, row_groups = row_groups,
                                      row_counts = NULL, group_pal = group_pal)
    }
  }
  if (nzchar(gsea_html)) {
    parts <- c(parts, sprintf(
      paste0('<details class="pic-enrich-group" open><summary>GSEA &mdash; gene set enrichment</summary>%s',
             '<p class="pic-note">GSEA asks which biological terms (GO, pathways&hellip;) are shifted up or down in a comparison, using the whole ranked gene list. Each dot is a term: position is its enrichment score (NES), color its significance (padj), size the number of genes.</p>',
             pick_hint('Select what to view &mdash; tick a <b>Contrast</b> cell and a <b>Method</b> below.'),
             '%s</details>'),
      src_note(report_rel_path(gsea_root, proj_dir)),
      gsea_html
    ))
  }

  # ---- ORA: method ごと (cluster をまとめて 1 枚) ----
  ora_root <- file.path(enrich_dir, "ORA")
  ora_items <- list()
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
      spec <- tryCatch(build_ora_plotly(ora_all), error = function(e) NULL)
      if (is.null(spec)) next
      pid <- sprintf("%soraplot_%s", id_prefix, method)
      register_plot(reg, pid, list(data = spec$data, layout = spec$layout, config = spec$config))
      ora_items[[length(ora_items) + 1L]] <- list(
        id = sprintf("%sora_%s", id_prefix, method),
        label = method,
        html = sprintf('<div class="pic-enrich-item"><h4>%s</h4><div id="%s" class="pic-plot" style="height:%dpx"></div></div>',
                       html_escape(method), pid, enrich_plot_height(spec$n_terms)),
        checked = (length(ora_items) == 0L)
      )
    }
  }
  if (length(ora_items) > 0) {
    ora_inner <- character(0)
    # Cluster expression profiles first (so cluster_N labels are interpretable)
    if (!is.null(deseq2_dir)) {
      prof_spec <- tryCatch(build_cluster_profile_plotly(deseq2_dir, project, group_pal), error = function(e) NULL)
      if (!is.null(prof_spec)) {
        prof_csv <- file.path(deseq2_dir, "DEG", "DEGCluster", sprintf("DEGCluster_profile_%s.csv", project))
        prof_id <- sprintf("%sclusterprofile", id_prefix)
        register_plot(reg, prof_id, list(data = prof_spec$data, layout = prof_spec$layout, config = prof_spec$config))
        ora_inner <- c(ora_inner, sprintf(
          paste0('<div class="pic-enrich-item"><h4>Cluster expression profiles</h4>',
                 '<p class="pic-note">The differentially expressed genes are grouped into clusters that share a similar pattern across groups. ',
                 'Each panel is one cluster; the boxes show its genes&rsquo; expression per group, so you can read off what each cluster represents.</p>',
                 '%s<div id="%s" class="pic-plot" style="height:%dpx"></div></div>'),
          src_note(report_rel_path(prof_csv, proj_dir)), prof_id, prof_spec$height
        ))
      }
    }
    ora_inner <- c(ora_inner,
      sprintf(paste0('<details class="pic-enrich-method" open><summary>Enriched terms per cluster</summary>',
                     '<p class="pic-note">For each cluster above, the biological terms over-represented among its genes. ',
                     'Dot color is significance (padj), size the gene ratio. Tick a term set to view it.</p>%s%s</details>'),
              src_note(report_rel_path(ora_root, proj_dir)), build_select_group(ora_items)))
    parts <- c(parts, sprintf(
      paste0('<details class="pic-enrich-group" open><summary>ORA &mdash; enrichment of gene clusters</summary>',
             '%s</details>'),
      paste(ora_inner, collapse = "")
    ))
  }

  paste(parts, collapse = "\n")
}

# 通常レポート用の enrichment セクション (<section id="enrich"> でラップ)。
section_enrichment <- function(enrich_dir, project, tmp_dir, deg_counts = NULL, deseq2_dir = NULL, group_pal = NULL, heading = "4. Enrichment", proj_dir = NULL, id_prefix = "", reg = NULL, group_order = NULL) {
  inner <- enrichment_blocks(enrich_dir, project, tmp_dir, deg_counts, deseq2_dir, group_pal, proj_dir, id_prefix, reg, group_order)
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

  # sample_sheet の順を正順とする (なければ mapping_sum の順)。
  sheet <- pic_read_sample_sheet(out_dir)
  sample_order <- if (!is.null(sheet)) sheet$samples else NULL
  group_order  <- if (!is.null(sheet) && length(sheet$groups) > 0) sheet$groups else NULL
  if (!is.null(msum)) msum <- pic_reorder_rows(msum, "sample", sample_order)

  # サンプル -> グループの対応とグループ配色 (PCA / heatmap で共有)
  group_map <- NULL
  group_pal <- NULL
  if (!is.null(msum) && all(c("sample", "group") %in% colnames(msum))) {
    group_map <- stats::setNames(as.character(msum$group), as.character(msum$sample))
    go <- if (!is.null(group_order)) pic_reorder_vec(unique(as.character(msum$group)), group_order)
          else unique(as.character(msum$group))
    group_pal <- group_palette(go)
  }

  # セクション組み立て (QC -> Correlation -> Aggregation -> PCA -> DEG -> Expr -> Enrichment)
  msum_file <- { mf <- list.files(out_dir, pattern = "^mapping_sum__.*\\.tsv$", full.names = TRUE); if (length(mf) > 0) mf[[1]] else "" }
  sec_qc <- if (!is.null(msum)) section_mapping_qc(msum, "qc", "1. Mapping QC", report_rel_path(msum_file, out_dir), group_map, group_pal) else ""

  # 2. Aggregation (TSS-TES)
  agg_html <- build_aggregate_html(reg, out_dir, "", out_dir, sample_order, group_pal)
  sec_agg <- if (nzchar(agg_html)) paste0('<section id="aggregate"><h2>2. Aggregation (TSS-TES)</h2>', agg_html, '</section>') else ""

  # 3. Sample correlation
  cor_html <- build_correlation_html(reg, desc$deseq2_dir, project, group_map, group_pal, sample_order, "", out_dir)
  sec_cor <- if (nzchar(cor_html)) paste0('<section id="cor"><h2>3. Sample correlation</h2>', cor_html, '</section>') else ""

  # 4. PCA
  pca_html <- build_pca_plots(reg, desc$deseq2_dir, project, group_pal, "", out_dir)
  sec_pca <- paste0('<section id="pca"><h2>4. PCA</h2>',
                    if (!is.null(pca_html)) pca_html else "<p>No PCA data.</p>",
                    '</section>')

  # 5. DEG (heatmap + MA/volcano)
  hm_html <- build_heatmap_html(desc$deseq2_dir, project, group_map, group_pal, out_dir, sample_order)
  contrast_html <- build_contrast_plots(reg, stats, deg_counts, fdr, "", report_rel_path(desc$stats_csv, out_dir), group_order, group_pal)
  sec_deg <- paste0(sprintf('<section id="deg"><h2>5. Differential expression (DESeq2, FDR = %s)</h2>', format(fdr, trim = TRUE)),
                    if (!is.null(hm_html)) hm_html else "",
                    contrast_html,
                    '</section>')

  # 6. Gene expression (normalized counts, interactive)
  sec_expr <- build_expression_section(reg, desc$deseq2_dir, project, group_map, group_pal, stats, "expr", "expr", "6. Gene expression", out_dir, sample_order, group_order)

  # 7. Enrichment
  sec_enrich <- section_enrichment(desc$enrich_dir, project, tmp_dir, deg_counts, desc$deseq2_dir, group_pal, "7. Enrichment", out_dir, "", reg, group_order)
  if (is.null(sec_enrich)) sec_enrich <- ""

  nav <- paste0('<nav class="pic-nav"><a href="#qc">QC</a>',
                if (nzchar(sec_agg)) '<a href="#aggregate">Aggregation</a>' else "",
                if (nzchar(sec_cor)) '<a href="#cor">Correlation</a>' else "",
                '<a href="#pca">PCA</a><a href="#deg">DEG</a>',
                if (nzchar(sec_expr)) '<a href="#expr">Expression</a>' else "",
                '<a href="#enrich">Enrichment</a></nav>')
  body <- paste0(sec_qc, sec_agg, sec_cor, sec_pca, sec_deg, sec_expr, sec_enrich)
  out_html <- file.path(out_dir, sprintf("report_%s.html", project))
  render_report_page(project, nav, body, reg, asset_dir, out_html)
  unlink(tmp_dir, recursive = TRUE)
  out_html
}

# HTML ページを組み立てて書き出す (plotly を内包)。
render_report_page <- function(title, nav_html, body_html, reg, asset_dir, out_html) {
  plots_json <- jsonlite::toJSON(reg$plots, auto_unbox = TRUE, null = "null", na = "null", digits = 6)
  expr_json <- jsonlite::toJSON(if (length(reg$expr) > 0) reg$expr else stats::setNames(list(), character(0)),
                                auto_unbox = TRUE, null = "null", na = "null", digits = 6)
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
    '<script>var PIC_PLOTS=', plots_json, ';var PIC_EXPR=', expr_json, ';</script>',
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
section_xenograft_qc <- function(summary_df, src = "") {
  cats <- names(PIC_XENO_COLORS)
  have <- all(cats %in% colnames(summary_df))
  parts <- c('<section id="classification"><h2>1. Xenograft classification (xengsort)</h2>', src_note(src))
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

  # sample_sheet_<frac>.tsv の順を正順とする
  sheet <- pic_read_sample_sheet(out_dir, frac_key)
  sample_order <- if (!is.null(sheet)) sheet$samples else NULL
  group_order  <- if (!is.null(sheet) && length(sheet$groups) > 0) sheet$groups else NULL

  # 画分の group 配色は、その画分の mapping_sum から取得
  msum <- NULL; group_map <- NULL; group_pal <- NULL
  msum_files <- list.files(frac_dir, pattern = "^mapping_sum__.*\\.tsv$", full.names = TRUE)
  if (length(msum_files) > 0) {
    msum <- pic_reorder_rows(read_mapping_sum(msum_files[[1]]), "sample", sample_order)
    if (all(c("sample", "group") %in% colnames(msum))) {
      group_map <- stats::setNames(as.character(msum$group), as.character(msum$sample))
      go <- if (!is.null(group_order)) pic_reorder_vec(unique(as.character(msum$group)), group_order)
            else unique(as.character(msum$group))
      group_pal <- group_palette(go)
    }
  }

  secs <- character(0)
  navs <- c(sprintf('<span class="pic-nav-group">%s</span>', html_escape(label)))
  n <- start_num

  # Mapping QC
  if (!is.null(msum)) {
    n <- n + 1L
    sid <- paste0(frac_key, "_qc")
    secs <- c(secs, section_mapping_qc(msum, sid, sprintf("%d. %s · Mapping QC", n, label),
                                       report_rel_path(msum_files[[1]], out_dir), group_map, group_pal))
    navs <- c(navs, sprintf('<a href="#%s">Mapping</a>', sid))
  }
  # Aggregation (TSS-TES)
  agg_html <- build_aggregate_html(reg, frac_dir, id_prefix, out_dir, sample_order, group_pal)
  if (nzchar(agg_html)) {
    n <- n + 1L
    sid <- paste0(frac_key, "_agg")
    secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · Aggregation (TSS-TES)</h2>%s</section>', sid, n, label, agg_html))
    navs <- c(navs, sprintf('<a href="#%s">Aggregation</a>', sid))
  }
  # Sample correlation
  cor_html <- build_correlation_html(reg, desc$deseq2_dir, project, group_map, group_pal, sample_order, id_prefix, out_dir)
  if (nzchar(cor_html)) {
    n <- n + 1L
    sid <- paste0(frac_key, "_cor")
    secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · Sample correlation</h2>%s</section>', sid, n, label, cor_html))
    navs <- c(navs, sprintf('<a href="#%s">Correlation</a>', sid))
  }
  # PCA
  n <- n + 1L
  sid <- paste0(frac_key, "_pca")
  pca_html <- build_pca_plots(reg, desc$deseq2_dir, project, group_pal, id_prefix, out_dir)
  secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · PCA</h2>%s</section>',
                          sid, n, label, if (!is.null(pca_html)) pca_html else "<p>No PCA data.</p>"))
  navs <- c(navs, sprintf('<a href="#%s">PCA</a>', sid))
  # DEG (heatmap + MA/volcano)
  n <- n + 1L
  sid <- paste0(frac_key, "_deg")
  hm_html <- build_heatmap_html(desc$deseq2_dir, project, group_map, group_pal, out_dir, sample_order, id_prefix)
  contrast_html <- build_contrast_plots(reg, stats, deg_counts, fdr, id_prefix, report_rel_path(desc$stats_csv, out_dir), group_order, group_pal)
  secs <- c(secs, sprintf('<section id="%s"><h2>%d. %s · Differential expression (DESeq2, FDR = %s)</h2>%s%s</section>',
                          sid, n, label, format(fdr, trim = TRUE),
                          if (!is.null(hm_html)) hm_html else "", contrast_html))
  navs <- c(navs, sprintf('<a href="#%s">DESeq2</a>', sid))
  # Gene expression (normalized counts)
  sid <- paste0(frac_key, "_expr")
  expr_html <- build_expression_section(reg, desc$deseq2_dir, project, group_map, group_pal, stats,
                                        paste0(frac_key, "__expr"), sid, sprintf("%d. %s · Expression", n + 1L, label), out_dir,
                                        sample_order, group_order)
  if (nzchar(expr_html)) {
    n <- n + 1L
    secs <- c(secs, expr_html)
    navs <- c(navs, sprintf('<a href="#%s">Expression</a>', sid))
  }
  # Enrichment
  enrich_inner <- enrichment_blocks(desc$enrich_dir, project, tmp_dir, deg_counts, desc$deseq2_dir, group_pal, out_dir, id_prefix, reg, group_order)
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

  sec_qc <- section_xenograft_qc(summary_df, report_rel_path(summary_files[[1]], out_dir))
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
