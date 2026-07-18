# 役割: 整形/エスケープ・group 色/名・パレット・ファイル登録・選択 UI 等の共通ヘルパ。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。

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

# group 表示切替パネル (plotly 用: trace を on/off)。
group_toggle_panel <- function(group_pal) {
  if (is.null(group_pal) || length(group_pal) == 0) return("")
  items <- vapply(names(group_pal), function(g) sprintf(
    '<label class="pic-tgl"><input type="checkbox" class="pic-gtoggle" data-group="%s" checked><span class="pic-swatch" style="background:%s"></span>%s</label>',
    html_escape(g), unname(group_pal[[g]]), html_escape(g)), character(1))
  paste0('<h4>Groups</h4><div class="pic-tgl-list">', paste(items, collapse = ""), '</div>')
}

# sample 表示切替パネル (HTML テーブル用: 行/列を隠す)。ラベルは group 配色。
sample_toggle_panel <- function(samples, group_map = NULL, group_pal = NULL) {
  if (length(samples) == 0) return("")
  scol <- function(s) sample_color(s, group_map, group_pal)
  items <- vapply(samples, function(s) sprintf(
    '<label class="pic-tgl"><input type="checkbox" class="pic-stoggle" data-sample="%s" checked><span style="color:%s;font-weight:600">%s</span></label>',
    html_escape(s), scol(s), html_escape(s)), character(1))
  paste0('<h4>Samples</h4><div class="pic-tgl-list">', paste(items, collapse = ""), '</div>')
}

# 「ここで選択」を明示する目立つ案内バッジ。
pick_hint <- function(text) {
  sprintf('<div class="pic-pick"><span class="pic-pick-ic">&#128071;</span><span>%s</span></div>', text)
}

# サブタブ内タイトル行 (左=タイトル / 右=アクションボタン群)。title は HTML 済み。
sub_head <- function(title_html, ...) {
  actions <- paste(c(...), collapse = "")
  sprintf('<div class="pic-sub-hd"><h3>%s</h3><span class="pic-sub-hd-actions">%s</span></div>',
          title_html, actions)
}

# HTML 要素 (id=cap_id) を PNG 画像でダウンロードするボタン。
png_button <- function(cap_id, name) {
  sprintf('<button class="pic-png-btn" type="button" data-cap="%s" data-name="%s">Download PNG</button>',
          html_escape(cap_id), html_escape(name))
}

# ファイル名 (相対パス) から、その中身の短い説明を返す。
file_desc <- function(rel) {
  b <- basename(rel)
  d <- function(...) paste0(...)
  if (grepl("^xenograft_classify_summary", b)) return("Xenograft read classification summary (per-sample genome split: graft/host/both/neither/ambiguous).")
  if (grepl("^mapping_sum", b)) return("Per-sample mapping/QC summary (read fate, total reads, UMIs, genes).")
  if (grepl("^sample_sheet", b)) return("Input sample sheet (fastq prefix, barcode, sample, group).")
  if (grepl("^deftable", b)) return("DESeq2 definition table (count prefix, sample, group).")
  if (grepl("^Num_UMIs_genes", b)) return("Per-sample UMI and detected-gene counts.")
  if (grepl("^UMI_count", b)) return("Raw UMI count matrix (genes × samples).")
  if (grepl("^normalizedCountTable", b)) return("DESeq2 size-factor normalized counts (genes × samples).")
  if (grepl("^DEG_normalizedCountTable", b)) return("Normalized counts restricted to DEGs (heatmap input).")
  if (grepl("^stats", b)) return("DESeq2 differential-expression statistics for all contrasts.")
  if (grepl("^correlation", b)) return("Sample-to-sample correlation matrix.")
  if (grepl("PCA_RegLog", b)) return("PCA sample coordinates (regularized-log transform).")
  if (grepl("PCA.*[Vv]ar", b) || grepl("variance", b)) return("PCA variance explained per component.")
  if (grepl("DEGCluster_profile", b)) return("Per-cluster gene expression profile (values by group).")
  if (grepl("DEGCluster_gene_for_ora", b)) return("Gene-to-cluster assignment used as ORA input.")
  if (grepl("DEGCluster_merge_map", b)) return("Mapping of merged DEG clusters.")
  if (grepl("DEGCluster_summary", b)) return("DEG cluster summary (size and membership per cluster).")
  if (grepl("DEGCluster", b)) return("DEG clustering result.")
  if (grepl("^DEGList|DEG_", b)) return("Differentially expressed gene list.")
  if (grepl("aggregate_profile", b)) return("Gene-body aggregation profile (TSS→TES metagene).")
  m <- regmatches(b, regexpr("^GSEA_([A-Za-z0-9]+)_(.+?)_", b, perl = TRUE))
  if (grepl("^GSEA_", b)) {
    mm <- regmatches(b, regexec("^GSEA_([A-Za-z0-9]+)_(.+)_[^_]+_[^_]+\\.csv$", b))[[1]]
    if (length(mm) >= 3) return(d("GSEA result — ", mm[[2]], ", contrast ", gsub("_vs_", " / ", mm[[3]]), "."))
    return("GSEA enrichment result.")
  }
  if (grepl("^ORA_", b)) {
    mm <- regmatches(b, regexec("^ORA_([A-Za-z0-9]+)_", b))[[1]]
    if (length(mm) >= 2) return(d("ORA result — ", mm[[2]], " (all clusters combined)."))
    return("ORA enrichment result (all clusters).")
  }
  mm <- regmatches(b, regexec("^([A-Za-z0-9]+)_cluster_([0-9]+)_", b))[[1]]
  if (length(mm) >= 3) return(d("ORA result — ", mm[[2]], ", cluster ", mm[[3]], "."))
  if (grepl("\\.txt\\.gz$|counts", rel)) return("Per-sample gene count table.")
  "Source data table."
}

# ダウンロード対象から除外するファイル (冗長・派生的なもの)。
file_excluded <- function(rel) {
  b <- basename(rel)
  grepl("^sample_sheet", b) || grepl("^Num_UMIs_genes", b) ||
    grepl("^DEG_count", b) || grepl("^DEG_normalizedCountTable", b) ||
    grepl("^DEGCluster_summary", b) || grepl("^DEGCluster_gene_for_ora", b) ||
    grepl("^DEGCluster_merge_map", b) ||
    # ORA は method ごとの結合 CSV のみ提供し、cluster ごとの表は含めない
    # (enrich/<genome>/ORA/ や <frac>/enrich/... など genome ネスト構造にも対応)
    (grepl("(^|/)ORA/", rel) && grepl("_cluster_[0-9]", b))
}

# 相対パスからカタログ用のカテゴリ (見出し + 対応タブ id) を返す。タブ順に並べる。
# rel の分類。enrich/aggregate は genome ネスト (enrich/<genome>/GSEA) や
# xenograft の <frac>/enrich/... にも対応するためパス中の GSEA/ORA/aggregate で判定。
# xenograft (先頭が graft//host/) では title に画分を付けて Downloads で区別する。
file_category <- function(rel) {
  b <- basename(rel)
  if (grepl("^xenograft_classify_summary", b)) return(list(key = "split", title = "Genome split", tab = "genomesplit"))
  frac <- if (grepl("^(graft|host)/", rel)) sub("^(graft|host)/.*$", "\\1", rel) else ""
  base <-
    if (grepl("enrich", rel, fixed = TRUE) && grepl("(^|/)GSEA/", rel)) list(key = "gsea", title = "Enrichment — GSEA", tab = "enrich")
    else if (grepl("enrich", rel, fixed = TRUE) && grepl("(^|/)ORA/", rel)) list(key = "ora", title = "Enrichment — ORA", tab = "enrich")
    else if (grepl("(^|/)aggregate/", rel)) list(key = "agg", title = "Aggregation", tab = "aggregate")
    else if (grepl("correlation", b)) list(key = "cor", title = "Sample Correlation", tab = "cor")
    else if (grepl("PCA", b)) list(key = "pca", title = "PCA", tab = "pca")
    else if (grepl("^stats|DEG", b)) list(key = "deg", title = "Differential Expression", tab = "deg")
    else if (grepl("UMI_count|normalizedCountTable", b)) list(key = "expr", title = "Count Table", tab = "expr")
    else list(key = "mapping", title = "Mapping", tab = "qc")
  if (nzchar(frac)) base$title <- sprintf("%s — %s", base$title, frac)
  base
}

# ファイル (相対パス) を gzip+base64 で埋め込みレジストリに登録し、id を返す。
register_file <- function(rel) {
  reg <- getOption("pic.report.reg"); proj <- getOption("pic.report.projdir")
  if (is.null(reg) || is.null(proj) || is.null(rel) || !nzchar(rel)) return(NULL)
  if (file_excluded(rel)) return(NULL)
  if (is.null(reg$files)) reg$files <- list()
  id <- gsub("[^A-Za-z0-9]+", "_", rel)
  if (!is.null(reg$files[[id]])) return(id)
  abs <- file.path(proj, rel)
  if (!file.exists(abs) || dir.exists(abs)) return(NULL)
  gz <- tryCatch({
    raw <- readBin(abs, "raw", n = file.info(abs)$size)
    jsonlite::base64_enc(memCompress(raw, "gzip"))
  }, error = function(e) NULL)
  if (is.null(gz)) return(NULL)
  cat_i <- file_category(rel)
  reg$files[[id]] <- list(name = basename(rel), path = rel, desc = file_desc(rel),
                          cat = cat_i$title, tab = cat_i$tab, gz = gz)
  id
}

# メモリ上の文字列を仮想ファイルとして埋め込み登録し、id を返す (ORA 結合 CSV 用)。
register_virtual_file <- function(rel, content, desc = NULL) {
  reg <- getOption("pic.report.reg")
  if (is.null(reg) || is.null(content)) return(NULL)
  if (is.null(reg$files)) reg$files <- list()
  id <- gsub("[^A-Za-z0-9]+", "_", rel)
  if (!is.null(reg$files[[id]])) return(id)
  gz <- tryCatch(jsonlite::base64_enc(memCompress(charToRaw(content), "gzip")), error = function(e) NULL)
  if (is.null(gz)) return(NULL)
  cat_i <- file_category(rel)
  reg$files[[id]] <- list(name = basename(rel), path = rel,
                          desc = if (!is.null(desc)) desc else file_desc(rel),
                          cat = cat_i$title, tab = cat_i$tab, gz = gz)
  id
}

# 元データのダウンロードボタン。ファイルを HTML に埋め込み、クリックで DL (自己完結)。
src_note <- function(rel) {
  id <- register_file(rel)
  if (is.null(id)) return("")
  ext <- toupper(tools::file_ext(rel))
  label <- if (nzchar(ext)) sprintf("Download %s", ext) else "Download"
  sprintf('<span class="pic-src"><a class="pic-dlcsv" role="button" tabindex="0" data-file="%s">%s</a></span>',
          id, html_escape(label))
}

# group 名をその group 色に対応させる (大文字小文字は無視)。無ければ NULL。
group_color <- function(g, group_pal) {
  if (is.null(group_pal)) return(NULL)
  idx <- match(tolower(g), tolower(names(group_pal)))
  if (!is.na(idx)) unname(group_pal[[idx]]) else NULL
}

# group 色 (大文字小文字は無視)。見つからなければ default。
group_color_or <- function(g, group_pal, default = "#1f2933") {
  col <- group_color(g, group_pal)
  if (!is.null(col)) col else default
}

# group 名 -> deftable どおりの正規名 (大文字小文字は無視)。見つからなければそのまま。
group_name_of <- function(g, group_pal) {
  if (is.null(group_pal)) return(g)
  idx <- match(tolower(g), tolower(names(group_pal)))
  if (!is.na(idx)) names(group_pal)[[idx]] else g
}

# sample -> その group 色 (group_map 経由)。見つからなければ default。
sample_color <- function(s, group_map, group_pal, default = "#333333") {
  g <- if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s
  if (!is.null(group_pal) && g %in% names(group_pal)) unname(group_pal[[g]]) else default
}

# group_pal が無ければ groups から生成し、欠けている group を補完したパレットを返す。
ensure_palette <- function(group_pal, groups) {
  pal <- if (!is.null(group_pal)) group_pal else group_palette(groups)
  miss <- setdiff(groups, names(pal))
  if (length(miss) > 0) pal <- c(pal, group_palette(miss))
  pal
}

# group 名を group 色の <span> で装飾 (plotly の name/legend でも色付き文字になる)。
group_span <- function(g, group_pal) {
  col <- group_color(g, group_pal)
  if (!is.null(col)) sprintf('<span style="color:%s">%s</span>', col, html_escape(g)) else html_escape(g)
}

# "a_vs_b" 形式の contrast を group 色付きの "a / b" にする。
color_vs_contrast <- function(cf, group_pal) {
  parts <- strsplit(cf, "_vs_", fixed = TRUE)[[1]]
  paste(vapply(parts, function(p) group_span(p, group_pal), character(1)),
        collapse = ' <span style="color:var(--muted)">/</span> ')
}

# contrast 名 ("groupA / groupB") の各 group をその group 色で装飾する。
color_contrast <- function(contrast, group_pal) {
  if (is.null(contrast) || !nzchar(contrast)) return(html_escape(contrast))
  parts <- trimws(strsplit(contrast, "/", fixed = TRUE)[[1]])
  gl <- if (!is.null(group_pal)) tolower(names(group_pal)) else character(0)
  colored <- vapply(parts, function(p) {
    idx <- match(tolower(p), gl)
    col <- if (!is.na(idx)) unname(group_pal[[idx]]) else NULL
    if (!is.null(col)) sprintf('<span style="color:%s;font-weight:700">%s</span>', col, html_escape(p))
    else html_escape(p)
  }, character(1))
  paste(colored, collapse = ' <span style="color:var(--muted)">/</span> ')
}

# 複数ブロックをチェックボックスで表示選択する UI。既定は先頭のみ表示。
# items: list(list(id=, label=, html=, checked=logical), ...)
# ブロックが 1 つだけなら選択 UI は付けずにそのまま表示する。
build_select_group <- function(items, prefix = "", ctrl_title = NULL, view_header = "") {
  items <- Filter(function(it) !is.null(it$html) && nzchar(it$html), items)
  if (length(items) == 0) return(if (nzchar(prefix) || nzchar(view_header))
    sprintf('<div class="pic-select-item">%s%s</div>', view_header, prefix) else "")
  rname <- paste0("sel_", gsub("[^A-Za-z0-9_]", "", items[[1]]$id))
  bar <- vapply(items, function(it) {
    ck <- if (isTRUE(it$checked)) " checked" else ""
    sprintf('<label class="pic-select-chk"><input type="radio" name="%s" data-target="%s"%s>%s</label>',
            rname, html_escape(it$id), ck, html_escape(it$label))
  }, character(1))
  blk <- vapply(items, function(it) {
    hid <- if (isTRUE(it$checked)) "" else " hidden"
    sprintf('<div class="pic-select-item" id="%s"%s>%s</div>', html_escape(it$id), hid, it$html)
  }, character(1))
  head_html <- if (!is.null(ctrl_title)) sprintf('<h4>%s</h4>', html_escape(ctrl_title)) else ""
  # 左 = サブタブ + ラジオ選択 / 右 = タイトル行 + (prefix +) 選択ブロック
  paste0('<div class="pic-selgrid"><div class="pic-selgrid-ctrl"><!--SUBNAV-->', head_html,
         '<div class="pic-select-bar">', paste(bar, collapse = ""),
         '</div></div><div class="pic-selgrid-view">', view_header, prefix,
         paste(blk, collapse = ""), '</div></div>')
}

# contrast を group×group の行列で選ばせる UI。
# 各セル (group i × group j) に、その 2 群の比較があれば DEG 数 + チェックボックス。
# entries: list(list(aspect=, count=, checked=, attr=))
#   attr は checkbox に付与する属性 (例: 'data-target="id"' や 'data-axis="r" data-key="k"')。
# groups: 表示順の group ベクトル (contrast は大文字小文字を無視して照合)。
build_group_matrix <- function(entries, groups, group_pal = NULL, show_count = TRUE, radio_name = NULL) {
  if (length(groups) < 2 || length(entries) == 0) return("")
  gl <- tolower(groups)
  cellmap <- list()
  for (e in entries) {
    sp <- strsplit(e$aspect, " / ", fixed = TRUE)[[1]]
    a <- tolower(trimws(sp[[1]])); b <- if (length(sp) >= 2) tolower(trimws(sp[[2]])) else ""
    ia <- match(a, gl); ib <- match(b, gl)
    if (is.na(ia) || is.na(ib)) next
    # tooltip 用の表示名は deftable どおりの casing (contrast データは小文字化されている)
    e$aspect_disp <- if (length(sp) >= 2) paste0(groups[[ia]], " / ", groups[[ib]]) else groups[[ia]]
    cellmap[[paste0(min(ia, ib), "_", max(ia, ib))]] <- e
  }
  if (length(cellmap) == 0) return("")
  # DEG 数を白→赤のヒートマップ色に (数値は出さず色で多寡を示す)。sqrt スケールで低値も視認可能に。
  cnt_of <- function(e) suppressWarnings(as.numeric(e$count))
  vmax <- suppressWarnings(max(vapply(cellmap, function(e) { v <- cnt_of(e); if (is.finite(v)) v else NA_real_ }, numeric(1)), na.rm = TRUE))
  heat_col <- function(v) {
    if (!is.finite(v) || !is.finite(vmax) || vmax <= 0) return(NA_character_)
    t <- sqrt(max(0, min(1, v / vmax)))
    mix <- function(a, b) round(a + t * (b - a))
    sprintf("rgb(%d,%d,%d)", mix(255, 178), mix(255, 24), mix(255, 43))  # #ffffff -> #b2182b
  }
  hcol <- function(g) group_color_or(g, group_pal)
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
      inp <- if (!is.null(radio_name)) sprintf('type="radio" name="%s"', html_escape(radio_name)) else 'type="checkbox"'
      cv <- cnt_of(e); col <- heat_col(cv)
      bg <- if (!is.na(col)) sprintf(' style="background:%s"', col) else ""
      hcls <- if (!is.na(col)) " pic-cmx-heat" else ""
      disp <- if (!is.null(e$aspect_disp)) e$aspect_disp else e$aspect
      ttl <- if (is.finite(cv)) sprintf("%s — %s DEGs", disp, fmt_int(cv)) else disp
      cs <- c(cs, sprintf('<td class="pic-cmx-cell%s" title="%s"%s><label>%s<input %s %s%s></label></td>',
                          hcls, html_escape(ttl), bg, cnt, inp, e$attr, ck))
    }
    rows_html <- c(rows_html, sprintf('<tr>%s%s</tr>', rh, paste(cs, collapse = "")))
  }
  cls <- if (show_count) "pic-cmatrix" else "pic-cmatrix pic-cmx-nocount"
  sprintf('<div class="pic-cmx-wrap"><table class="%s"><thead>%s</thead><tbody>%s</tbody></table></div>',
          cls, header, paste(rows_html, collapse = ""))
}

# 2 軸 (行 = contrast, 列 = method) をチェックボックスで選び、両方選択のセルのみ表示。
# row_groups を渡すと、行 (contrast) の選択を group×group 行列で提示する。
build_matrix_group <- function(rows, cols, cell_html, row_title = "Comparison", col_title = "Method",
                               row_groups = NULL, row_counts = NULL, group_pal = NULL, view_header = "",
                               id_prefix = "") {
  if (length(rows) == 0 || length(cols) == 0) return("")
  first_row <- rows[[1]]$key
  first_col <- cols[[1]]
  # 行 (contrast) セレクタ: group 行列 or 線形バー。radio name は id_prefix で一意化
  # (xenograft で graft/host の GSEA 行列が同名ラジオになり相互干渉するのを防ぐ)
  rname <- paste0(id_prefix, "gseacon_", gsub("[^A-Za-z0-9_]", "", rows[[1]]$key))
  gm <- ""
  if (!is.null(row_groups) && length(row_groups) >= 2) {
    entries <- lapply(seq_along(rows), function(i) list(
      aspect = rows[[i]]$label,
      count = if (!is.null(row_counts)) row_counts[[rows[[i]]$key]] else NULL,
      checked = (i == 1L),
      attr = sprintf('data-axis="r" data-key="%s"', html_escape(rows[[i]]$key))))
    gm <- build_group_matrix(entries, row_groups, group_pal, show_count = FALSE, radio_name = rname)
  }
  row_sel <- if (nzchar(gm)) gm else {
    row_bar <- vapply(seq_along(rows), function(i) {
      ck <- if (i == 1L) " checked" else ""
      sprintf('<label class="pic-select-chk"><input type="radio" name="%s" data-axis="r"%s data-key="%s">%s</label>',
              rname, ck, html_escape(rows[[i]]$key), html_escape(rows[[i]]$label))
    }, character(1))
    sprintf('<div class="pic-select-bar"><span class="pic-select-lbl">%s:</span>%s</div>',
            html_escape(row_title), paste(row_bar, collapse = ""))
  }
  cname <- paste0(id_prefix, "gseameth_", gsub("[^A-Za-z0-9_]", "", cols[[1]]))
  col_bar <- vapply(seq_along(cols), function(i) {
    ck <- if (i == 1L) " checked" else ""
    sprintf('<label class="pic-select-chk"><input type="radio" name="%s" data-axis="c"%s data-key="%s">%s</label>',
            cname, ck, html_escape(cols[[i]]), html_escape(cols[[i]]))
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
  # 左 (.pic-matrix-sel) = サブタブ + Comparison(行列) + Method(見出し+縦並び)、右 = プロット
  comp_hdr <- if (nzchar(gm)) sprintf('<h4>%s</h4>', html_escape(row_title)) else ""
  sprintf(paste0('<div class="pic-matrix"><div class="pic-matrix-sel"><!--SUBNAV-->%s%s',
                 '<h4>%s</h4><div class="pic-select-bar">%s</div></div>',
                 '<div class="pic-matrix-cells">%s%s</div></div>'),
          comp_hdr, row_sel, html_escape(col_title), paste(col_bar, collapse = ""),
          view_header, paste(cells, collapse = ""))
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

