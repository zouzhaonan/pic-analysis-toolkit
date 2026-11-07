# 役割: Enrichment (GSEA / ORA) セクション。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。

# GSEA の CSV の一覧 (method のフォルダ順、その中は名前順)。enrichment_blocks の読む順と同じ。
pic_gsea_csvs <- function(enrich_dir, project) {
  gsea_root <- file.path(enrich_dir, "GSEA")
  if (is.null(enrich_dir) || is.na(enrich_dir) || !dir.exists(gsea_root)) return(list())
  out <- list()
  for (mdir in list.dirs(gsea_root, recursive = FALSE, full.names = TRUE)) {
    csvs <- list.files(mdir, pattern = paste0(project, "\\.csv$"), full.names = TRUE)
    for (csv in sort(csvs)) out[[length(out) + 1L]] <- list(csv = csv, method = basename(mdir))
  }
  out
}

# GSEA の CSV 1 本を読み、図の spec を作る。図にしないもの (空・除外・spec の失敗) は NULL。
# plot は render と同じ JSON にしておく (並べた子から親へ渡すときに軽く、render でも書き直さない)。
pic_gsea_item <- function(csv, method, group_pal) {
  df <- suppressMessages(readr::read_csv(csv, show_col_types = FALSE, progress = FALSE))
  if (nrow(df) == 0 || !("contrast" %in% colnames(df))) return(NULL)
  contrast <- as.character(df$contrast[[1]])
  # データ固有の除外: cntl_nega / cntl_posi の REACTOME プロットは表示しない
  if (identical(tolower(contrast), "cntl_nega / cntl_posi") && identical(toupper(method), "REACTOME")) return(NULL)
  sp <- strsplit(contrast, " / ", fixed = TRUE)[[1]]
  numerator <- if (length(sp) >= 1) sp[[1]] else ""
  denominator <- if (length(sp) >= 2) sp[[2]] else ""
  spec <- tryCatch(build_gsea_plotly(df, numerator, denominator, group_pal), error = function(e) NULL)
  if (is.null(spec)) return(NULL)
  list(contrast = contrast, n_terms = spec$n_terms,
       plot = pic_json_verbatim(list(data = spec$data, layout = spec$layout, config = spec$config), digits = 6))
}

# GSEA の CSV 1 本分の仕事 (並べる子で回す): 図の spec と、Download CSV で埋め込む gzip+base64。
# gz は register_file が見る置き場の鍵 (projdir/<相対パス>) と一緒に返し、親が置き場に入れる。
pic_gsea_job <- function(csv, method, group_pal, proj_dir) {
  it <- pic_gsea_item(csv, method, group_pal)
  gz <- NULL
  proj <- getOption("pic.report.projdir")
  if (!is.null(it) && !is.null(proj)) {
    rel <- report_rel_path(csv, proj_dir)
    abs <- file.path(proj, rel)
    if (nzchar(rel) && !file_excluded(rel) && file.exists(abs) && !dir.exists(abs))
      gz <- list(key = abs, gz = pic_report_file_gz(abs))
  }
  list(item = it, gz = gz)
}

# ORA の method のフォルダ (cluster ごとの CSV をまとめて読む)。
pic_ora_mdirs <- function(enrich_dir) {
  ora_root <- file.path(enrich_dir, "ORA")
  if (is.null(enrich_dir) || is.na(enrich_dir) || !dir.exists(ora_root)) return(character(0))
  list.dirs(ora_root, recursive = FALSE, full.names = TRUE)
}

# ORA の method 1 つ分: cluster の CSV を合わせて読み、図の spec と結合 CSV の中身を作る。図にしないものは NULL。
pic_ora_item <- function(mdir, project) {
  csvs <- list.files(mdir, pattern = paste0(project, "\\.csv$"), full.names = TRUE)
  csvs <- csvs[!grepl("^ORA_", basename(csvs))]   # 結合 CSV は cluster から集計するため除外
  if (length(csvs) == 0) return(NULL)
  ora_all <- tryCatch(
    purrr::map_dfr(csvs, function(fp) suppressMessages(readr::read_csv(fp, show_col_types = FALSE, progress = FALSE))),
    error = function(e) NULL
  )
  if (is.null(ora_all) || nrow(ora_all) == 0) return(NULL)
  spec <- tryCatch(build_ora_plotly(ora_all), error = function(e) NULL)
  if (is.null(spec)) return(NULL)
  list(n_terms = spec$n_terms,
       content = tryCatch(readr::format_csv(ora_all), error = function(e) NULL),
       plot = pic_json_verbatim(list(data = spec$data, layout = spec$layout, config = spec$config), digits = 6))
}

# enrichment セクションの中身 (h3 GSEA / ORA) を返す。<section> ラッパは付けない。
# pre: 先に (並べて) 作った list(gsea = pic_gsea_item の list (pic_gsea_csvs の順), ora = pic_ora_item の list
#      (pic_ora_mdirs の順))。NULL ならここで順に作る。
# tree_html: GSEA / ORA の左パネルの末尾に足す HTML (サンプル除外のパネル)。
enrichment_blocks <- function(enrich_dir, project, tmp_dir, deg_counts = NULL, deseq2_dir = NULL, group_pal = NULL, proj_dir = NULL, id_prefix = "", reg = NULL, group_order = NULL, pre = NULL,
                              tree_html = "") {
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
    cellmap <- list()          # cellmap[[contrast]][[method]] = uri
    methods_seen <- character(0)
    gcsvs <- pic_gsea_csvs(enrich_dir, project)
    gitems <- if (!is.null(pre)) pre$gsea else lapply(gcsvs, function(g) pic_gsea_item(g$csv, g$method, group_pal))
    for (k in seq_along(gcsvs)) {
      it <- gitems[[k]]
      if (is.null(it)) next
      csv <- gcsvs[[k]]$csv; method <- gcsvs[[k]]$method; contrast <- it$contrast
      pid <- sprintf("%sgseaplot_%s_%s", id_prefix, method, format_contrast_file_label(contrast))
      register_plot(reg, pid, it$plot)
      if (is.null(cellmap[[contrast]])) cellmap[[contrast]] <- list()
      cellmap[[contrast]][[method]] <- list(id = pid, h = enrich_plot_height(it$n_terms),
                                            csv = report_rel_path(csv, proj_dir))
      methods_seen <- union(methods_seen, method)
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
        # タイトル行は不要 (左パネルで contrast/method を選択済み)。DL ボタンのみ右寄せ。
        dl <- if (!is.null(info$csv) && nzchar(info$csv)) src_note(info$csv) else ""
        sprintf('<div class="pic-enrich-item"><div class="pic-dl-row">%s</div><div id="%s" class="pic-plot" style="height:%dpx"></div></div>',
                dl, info$id, info$h)
      }
      row_groups <- if (!is.null(group_order) && length(group_order) > 0) group_order else
        unique(unlist(lapply(contrasts, function(a) trimws(strsplit(a, " / ", fixed = TRUE)[[1]]))))
      # DEG 数 (contrast -> key) を渡し、セルをヒートマップ着色 (GSEA も MA/Volcano と共通)
      row_counts <- stats::setNames(
        lapply(contrasts, function(cc) if (!is.null(deg_counts) && cc %in% names(deg_counts)) deg_counts[[cc]] else NA_real_),
        vapply(contrasts, format_contrast_file_label, character(1)))
      gsea_html <- build_matrix_group(rows, methods, cell_html, row_groups = row_groups,
                                      row_counts = row_counts, group_pal = group_pal, view_header = "",
                                      id_prefix = id_prefix, ctrl_extra = tree_html)
    }
  }
  if (nzchar(gsea_html)) {
    parts <- c(parts, sprintf(
      paste0('<div class="pic-enrich-gsea">',
             '<p class="pic-note">GSEA asks which biological terms (GO, pathways&hellip;) are shifted up or down in a comparison, using the whole ranked gene list. Each dot is a term: position is its enrichment score (NES), color its significance (padj), size the number of genes.</p>',
             '%s</div>'),
      gsea_html
    ))
  }

  # ---- ORA: method ごと (cluster をまとめて 1 枚) ----
  ora_root <- file.path(enrich_dir, "ORA")
  ora_items <- list()
  if (dir.exists(ora_root)) {
    methods <- pic_ora_mdirs(enrich_dir)
    oitems <- if (!is.null(pre)) pre$ora else lapply(methods, pic_ora_item, project = project)
    for (k in seq_along(methods)) {
      it <- oitems[[k]]
      if (is.null(it)) next
      method <- basename(methods[[k]])
      # method ごとに cluster をまとめた 1 本の CSV を仮想ファイルとして埋め込み (ワンクリック DL)
      # xenograft では id_prefix ("graft__" 等) から画分ディレクトリを付け、Downloads で画分を区別
      frac_dir <- if (nzchar(id_prefix)) paste0(sub("__$", "", id_prefix), "/") else ""
      vrel <- sprintf("%senrich/ORA/%s/ORA_%s_%s.csv", frac_dir, method, method, project)
      vid <- if (is.null(it$content)) NULL else
        tryCatch(register_virtual_file(vrel, it$content,
                                       desc = sprintf("ORA result — %s (all clusters combined).", method)),
                 error = function(e) NULL)
      dl <- if (!is.null(vid))
        sprintf('<span class="pic-src"><a class="pic-dlcsv" role="button" tabindex="0" data-file="%s">Download CSV</a></span>', vid)
        else ""
      pid <- sprintf("%soraplot_%s", id_prefix, method)
      register_plot(reg, pid, it$plot)
      ora_items[[length(ora_items) + 1L]] <- list(
        id = sprintf("%sora_%s", id_prefix, method),
        label = method,
        html = sprintf('<div class="pic-enrich-item">%s<div id="%s" class="pic-plot" style="height:%dpx"></div></div>',
                       sub_head(html_escape(method), dl), pid, enrich_plot_height(it$n_terms)),
        checked = (length(ora_items) == 0L)
      )
    }
  }
  if (length(ora_items) > 0) {
    # 右ビューの先頭に常時表示する cluster expression profiles (cluster_N ラベルの意味づけ)
    prefix_html <- ""
    if (!is.null(deseq2_dir)) {
      prof_spec <- tryCatch(build_cluster_profile_plotly(deseq2_dir, project, group_pal), error = function(e) NULL)
      if (!is.null(prof_spec)) {
        prof_csv <- file.path(degcluster_dir_of(deseq2_dir, project), sprintf("DEGCluster_profile_%s.csv", project))
        prof_id <- sprintf("%sclusterprofile", id_prefix)
        register_plot(reg, prof_id, list(data = prof_spec$data, layout = prof_spec$layout, config = prof_spec$config))
        prefix_html <- sprintf(
          paste0('<div class="pic-enrich-item">%s',
                 '<p class="pic-note">The differentially expressed genes are grouped into clusters that share a similar pattern across groups. ',
                 'Each panel is one cluster; the boxes show its genes&rsquo; expression per group, so you can read off what each cluster represents.</p>',
                 '<div id="%s" class="pic-plot" style="height:%dpx"></div></div>'),
          sub_head('Cluster expression profiles', src_note(report_rel_path(prof_csv, proj_dir))),
          prof_id, prof_spec$height
        )
      }
    }
    parts <- c(parts, sprintf(
      paste0('<div class="pic-enrich-ora">',
             '<p class="pic-note">ORA takes the differentially expressed genes in each cluster and asks which biological terms are over-represented among them. ',
             'Dot color is significance (p.adjust), size the gene ratio.</p>%s</div>'),
      build_select_group(ora_items, prefix = prefix_html, ctrl_title = "Method", view_header = "", ctrl_extra = tree_html,
                         method_src = TRUE)
    ))
  } else {
    # ORA が無い場合は無言でスキップせず、理由を明記する。
    # 判定材料: DEG クラスタ (ORA の入力) が生成されているか。
    degc <- if (!is.null(deseq2_dir))
      file.path(degcluster_dir_of(deseq2_dir, project), sprintf("DEGCluster_gene_for_ora_%s.csv", project)) else ""
    has_clusters <- nzchar(degc) && file.exists(degc)
    reason <- if (!has_clusters) {
      paste0("No DEG clusters were produced for this dataset &mdash; there were too few differentially expressed genes to form clusters. ",
             "ORA runs on clusters of DEGs, so it was not performed here.")
    } else if (!dir.exists(ora_root)) {
      paste0("DEG clusters were produced, but no ORA results are present. ",
             "ORA was not run for this dataset (it needs the DEG clusters passed to <code>pic enrich --deg-clusters</code>).")
    } else {
      "DEG clusters were produced, but ORA found no over-represented terms at the current significance threshold."
    }
    parts <- c(parts, sprintf(
      paste0('<div class="pic-enrich-ora">%s',
             '<p class="pic-note">ORA (over-representation analysis) asks which biological terms are enriched among the differentially expressed genes in each cluster. ',
             '%s</p></div>'),
      sub_head("ORA", ""), reason
    ))
  }

  paste(parts, collapse = "\n")
}

# 通常レポート用の enrichment セクション (<section id="enrich"> でラップ)。
section_enrichment <- function(enrich_dir, project, tmp_dir, deg_counts = NULL, deseq2_dir = NULL, group_pal = NULL, heading = "4. Enrichment", proj_dir = NULL, id_prefix = "", reg = NULL, group_order = NULL, pre = NULL,
                               tree_html = "") {
  inner <- enrichment_blocks(enrich_dir, project, tmp_dir, deg_counts, deseq2_dir, group_pal, proj_dir, id_prefix, reg, group_order, pre, tree_html)
  if (!nzchar(inner)) inner <- '<p>No enrichment plots were generated.</p>'
  sprintf('<section id="enrich"><h2>%s</h2>%s</section>', heading, inner)
}

