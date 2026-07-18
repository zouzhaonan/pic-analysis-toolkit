# 役割: Enrichment (GSEA / ORA) セクション。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。

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
        # データ固有の除外: cntl_nega / cntl_posi の REACTOME プロットは表示しない
        if (identical(tolower(contrast), "cntl_nega / cntl_posi") && identical(toupper(method), "REACTOME")) next
        sp <- strsplit(contrast, " / ", fixed = TRUE)[[1]]
        numerator <- if (length(sp) >= 1) sp[[1]] else ""
        denominator <- if (length(sp) >= 2) sp[[2]] else ""
        spec <- tryCatch(build_gsea_plotly(df, numerator, denominator, group_pal), error = function(e) NULL)
        if (is.null(spec)) next
        pid <- sprintf("%sgseaplot_%s_%s", id_prefix, method, format_contrast_file_label(contrast))
        register_plot(reg, pid, list(data = spec$data, layout = spec$layout, config = spec$config))
        if (is.null(cellmap[[contrast]])) cellmap[[contrast]] <- list()
        cellmap[[contrast]][[method]] <- list(id = pid, h = enrich_plot_height(spec$n_terms),
                                              csv = report_rel_path(csv, proj_dir))
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
                                      row_counts = row_counts, group_pal = group_pal, view_header = "")
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
    methods <- list.dirs(ora_root, recursive = FALSE, full.names = TRUE)
    for (mdir in methods) {
      method <- basename(mdir)
      csvs <- list.files(mdir, pattern = paste0(project, "\\.csv$"), full.names = TRUE)
      csvs <- csvs[!grepl("^ORA_", basename(csvs))]   # 結合 CSV は cluster から集計するため除外
      if (length(csvs) == 0) next
      ora_all <- tryCatch(
        purrr::map_dfr(csvs, function(fp) suppressMessages(readr::read_csv(fp, show_col_types = FALSE, progress = FALSE))),
        error = function(e) NULL
      )
      if (is.null(ora_all) || nrow(ora_all) == 0) next
      spec <- tryCatch(build_ora_plotly(ora_all), error = function(e) NULL)
      if (is.null(spec)) next
      # method ごとに cluster をまとめた 1 本の CSV を仮想ファイルとして埋め込み (ワンクリック DL)
      vrel <- sprintf("enrich/ORA/%s/ORA_%s_%s.csv", method, method, project)
      vid <- tryCatch(register_virtual_file(vrel, readr::format_csv(ora_all),
                                            desc = sprintf("ORA result — %s (all clusters combined).", method)),
                      error = function(e) NULL)
      dl <- if (!is.null(vid))
        sprintf('<span class="pic-src"><a class="pic-dlcsv" role="button" tabindex="0" data-file="%s">Download CSV</a></span>', vid)
        else ""
      pid <- sprintf("%soraplot_%s", id_prefix, method)
      register_plot(reg, pid, list(data = spec$data, layout = spec$layout, config = spec$config))
      ora_items[[length(ora_items) + 1L]] <- list(
        id = sprintf("%sora_%s", id_prefix, method),
        label = method,
        html = sprintf('<div class="pic-enrich-item">%s<div id="%s" class="pic-plot" style="height:%dpx"></div></div>',
                       sub_head(html_escape(method), dl), pid, enrich_plot_height(spec$n_terms)),
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
        prof_csv <- file.path(deseq2_dir, "DEG", "DEGCluster", sprintf("DEGCluster_profile_%s.csv", project))
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
             'Dot color is significance (p.adjust), size the gene ratio. Pick a method on the left to view it.</p>%s</div>'),
      build_select_group(ora_items, prefix = prefix_html, ctrl_title = "Method", view_header = "")
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

