# 役割: CSS/JS・タブ/ページ組み立て・Overview/Downloads・orchestration・xenograft。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。

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

# DEG 数 (contrast -> 合計) を DEG_count.csv から取得。無ければ stats から算出する。
# (通常レポート / xenograft 画分の両経路で共通)
resolve_deg_counts <- function(desc, stats, fdr) {
  dc <- pic_load_deg_counts(desc$deseq2_dir, desc$project)
  if (is.null(dc)) {
    dc <- list()
    for (a in unique(stats$aspect)) dc[[a]] <- sum(stats$aspect == a & !is.na(stats$padj) & stats$padj < fdr)
  }
  dc
}

# msum から sample->group 対応と group 配色を得る。列が無ければ両方 NULL。
# (通常レポート / xenograft 画分の両経路で共通)
group_map_pal <- function(msum, group_order) {
  if (is.null(msum) || !all(c("sample", "group") %in% colnames(msum))) return(list(map = NULL, pal = NULL))
  gmap <- stats::setNames(as.character(msum$group), as.character(msum$sample))
  go <- if (!is.null(group_order)) pic_reorder_vec(unique(as.character(msum$group)), group_order)
        else unique(as.character(msum$group))
  list(map = gmap, pal = group_palette(go))
}

# 解析タブ (QC/Aggregation/Correlation/PCA/DEG/Expression/Enrichment) を組み立てて
# data_tabs (build_tabs 入力) を返す。id_prefix を付けると section / sub / plot の
# ID が prefix 付きになり、同一ページに複数 (xenograft の graft/host) を共存できる。
# base_dir: mapping_sum / aggregate が置かれるディレクトリ (通常は out_dir、画分では frac_dir)。
build_analysis_data_tabs <- function(reg, desc, base_dir, out_dir, id_prefix, fdr,
                                     msum, group_map, group_pal, sample_order, group_order,
                                     stats, deg_counts, tmp_dir) {
  sid <- function(k) paste0(id_prefix, k)
  project <- desc$project

  # 1. Mapping QC
  msum_file <- { mf <- pic_list_summary(base_dir, "^mapping_sum__.*\\.tsv$"); if (length(mf) > 0) mf[[1]] else "" }
  sec_qc <- if (!is.null(msum)) section_mapping_qc(msum, sid("qc"), "1. Mapping QC", report_rel_path(msum_file, out_dir), group_map, group_pal) else ""

  # 2. Aggregation (TSS-TES)
  agg_html <- build_aggregate_html(reg, base_dir, id_prefix, out_dir, sample_order, group_pal)
  sec_agg <- if (nzchar(agg_html)) paste0(sprintf('<section id="%s"><h2>2. Aggregation</h2>', sid("aggregate")), agg_html, '</section>') else ""

  # 3. Sample correlation
  cor_html <- build_correlation_html(reg, desc$deseq2_dir, project, group_map, group_pal, sample_order, id_prefix, out_dir)
  sec_cor <- if (nzchar(cor_html)) paste0(sprintf('<section id="%s"><h2>3. Sample Correlation</h2>', sid("cor")), cor_html, '</section>') else ""

  # 4. PCA
  pca_html <- build_pca_plots(reg, desc$deseq2_dir, project, group_pal, id_prefix, out_dir)
  sec_pca <- paste0(sprintf('<section id="%s"><h2>4. PCA</h2>', sid("pca")),
                    if (!is.null(pca_html)) pca_html else "<p>No PCA data.</p>",
                    '</section>')

  # 5. DEG (heatmap + MA/volcano)
  hm_html <- build_heatmap_html(desc$deseq2_dir, project, group_map, group_pal, out_dir, sample_order, id_prefix)
  contrast_html <- build_contrast_plots(reg, stats, deg_counts, fdr, id_prefix, report_rel_path(desc$stats_csv, out_dir), group_order, group_pal)
  sec_deg <- paste0(sprintf('<section id="%s"><h2>5. Differential Expression (FDR = %s)</h2>', sid("deg"), format(fdr, trim = TRUE)),
                    if (!is.null(hm_html)) hm_html else "",
                    contrast_html,
                    '</section>')

  # 6. Gene expression (normalized counts, interactive)
  sec_expr <- build_expression_section(reg, desc$deseq2_dir, project, group_map, group_pal, stats, sid("expr"), sid("expr"), "6. Gene Expression", out_dir, sample_order, group_order)

  # 7. Enrichment (section_enrichment は id="enrich" を埋め込むので prefix 付けは後処理)
  sec_enrich <- section_enrichment(desc$enrich_dir, project, tmp_dir, deg_counts, desc$deseq2_dir, group_pal, "7. Enrichment", out_dir, id_prefix, reg, group_order)
  if (is.null(sec_enrich)) sec_enrich <- ""
  sec_enrich <- sub('<section id="enrich">', sprintf('<section id="%s">', sid("enrich")), sec_enrich, fixed = TRUE)

  grp_panel <- group_toggle_panel(group_pal)
  list(
    list(html = sec_qc, label = "QC", shared_controls = grp_panel, subs = list(
      list(label = "Read distribution", id = sid("qc-readdist"), marker = ""),
      list(label = "Sequencing depth",  id = sid("qc-depth"),    marker = '<div class="pic-sub-hd"><h3>Sequencing depth'))),
    list(html = sec_agg,    label = "Aggregation", controls = grp_panel),
    list(html = sec_cor,    label = "Correlation", controls = grp_panel),
    list(html = sec_pca,    label = "PCA",         controls = grp_panel),
    list(html = sec_deg, label = "DEG", subs = list(
      list(label = "Heatmap",     id = sid("deg-heatmap"),   marker = ""),
      list(label = "M-A / Volcano", id = sid("deg-mavolcano"), marker = '<div class="pic-degsel'))),
    list(html = sec_expr,   label = "Expression"),
    list(html = sec_enrich, label = "Enrichment", subs = list(
      list(label = "GSEA", id = sid("enrich-gsea"), marker = ""),
      list(label = "ORA",  id = sid("enrich-ora"),  marker = '<div class="pic-enrich-ora">')))
  )
}

build_report_for_project <- function(desc, out_dir, msum, asset_dir) {
  options(pic.report.asset_dir = asset_dir)
  fdr <- pic_plot_spec()$defaults$fdr
  project <- desc$project
  # project = <run>_<genome>。deseq2/ は平坦化済みなので run は mapping_sum から、
  # genome は sample_sheet の genome 列 (無ければ project から run 接頭辞を除去) で解決する。
  run <- pic_report_run_name(out_dir)
  if (!nzchar(run)) run <- project
  genome <- pic_report_genome_of(out_dir, project)
  reg <- pic_report_registry()
  # src_note / register_file が参照する埋め込みレジストリと基準ディレクトリ
  options(pic.report.reg = reg, pic.report.projdir = out_dir)

  stats <- suppressMessages(readr::read_csv(desc$stats_csv, show_col_types = FALSE, progress = FALSE))
  stats <- as.data.frame(stats, check.names = FALSE)

  # DEG 数 (contrast -> 合計) を DEG_count.csv から取得 (なければ stats から)
  deg_counts <- resolve_deg_counts(desc, stats, fdr)

  tmp_dir <- file.path(tempdir(), paste0("picreport_", project))
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  # sample_sheet の順を正順とする (なければ mapping_sum の順)。
  sheet <- pic_read_sample_sheet(out_dir)
  # sample_sheet に batch 列があれば、ラベル / マーカーの batch 表現を有効にする。
  pic_set_batch_info(if (!is.null(sheet)) sheet$batch else NULL,
                     if (!is.null(sheet)) sheet$batch_levels else NULL)
  sample_order <- if (!is.null(sheet)) sheet$samples else NULL
  group_order  <- if (!is.null(sheet) && length(sheet$groups) > 0) sheet$groups else NULL
  if (!is.null(msum)) msum <- pic_reorder_rows(msum, "sample", sample_order)

  # サンプル -> グループの対応とグループ配色 (PCA / heatmap で共有)
  mp <- group_map_pal(msum, group_order)
  group_map <- mp$map
  group_pal <- mp$pal

  data_tabs <- build_analysis_data_tabs(reg, desc, out_dir, out_dir, "", fdr,
                                        msum, group_map, group_pal, sample_order, group_order,
                                        stats, deg_counts, tmp_dir)
  overview_sec  <- build_overview_section(run, genome, params = pic_read_analysis_params(out_dir))
  downloads_sec <- build_downloads_section(data_tabs, run, genome, group_pal)
  tb <- build_tabs(c(list(list(html = overview_sec, label = "Overview")),
                     data_tabs,
                     list(list(html = downloads_sec, label = "Downloads"))))
  out_html <- file.path(out_dir, sprintf("report_%s.html", project))
  render_report_page(run, tb$bar, tb$panels, reg, asset_dir, out_html)
  unlink(tmp_dir, recursive = TRUE)
  out_html
}

# 文字列から source (pic-src) を抽出し、除去後の HTML と結合 src を返す。
pic_pull_src <- function(s) {
  srcs <- unique(unlist(regmatches(s, gregexpr('<span class="pic-src">.*?</a></span>', s, perl = TRUE))))
  list(src = paste(srcs, collapse = ""),
       html = gsub('<span class="pic-src">.*?</a></span>', '', s, perl = TRUE))
}

# 文字列から説明文 (<p class="pic-note">) を抽出し、除去後の HTML と結合 desc を返す。
pic_pull_notes <- function(s) {
  notes <- unlist(regmatches(s, gregexpr('<p class="pic-note">.*?</p>', s, perl = TRUE)))
  desc <- if (length(notes) > 0) paste(sub('<p class="pic-note">', '<p>', notes), collapse = "") else ""
  list(desc = desc, html = gsub('<p class="pic-note">.*?</p>', '', s, perl = TRUE))
}

# <section id><h2>title</h2>inner</section> から (id, title, inner) を取り出す。
# 説明文・source は inner に残し pic_tab_panel 側で処理する。
pic_extract_section <- function(sec_html) {
  id <- regmatches(sec_html, regexpr('(?<=<section id=")[^"]+', sec_html, perl = TRUE))
  title <- regmatches(sec_html, regexpr('(?<=<h2>).*?(?=</h2>)', sec_html, perl = TRUE))
  inner <- sub('^<section[^>]*><h2>.*?</h2>', '', sec_html)
  inner <- sub('</section>[[:space:]]*$', '', inner)
  list(id = id, title = title, inner = inner)
}

info_badge <- function(desc) {
  if (!nzchar(desc)) return("")
  sprintf('<span class="pic-info" tabindex="0">i<span class="pic-info-pop">%s</span></span>', desc)
}

# inner を markers (各 sub の開始文字列。先頭は "") で分割してチャンクを返す。
pic_split_inner <- function(inner, markers) {
  pos <- vapply(markers, function(mk) {
    if (!nzchar(mk)) return(1L)
    p <- regexpr(mk, inner, fixed = TRUE)[[1]]; if (p < 0) NA_integer_ else as.integer(p)
  }, integer(1))
  n <- length(markers)
  ends <- c(pos[-1], nchar(inner) + 1L)
  vapply(seq_len(n), function(i) {
    s <- pos[[i]]; e <- ends[[i]]
    if (is.na(s)) return("")
    if (is.na(e)) e <- nchar(inner) + 1L
    substr(inner, s, e - 1L)
  }, character(1))
}

# 1 タブを組み立てる。subs 指定時はセクション内サブタブ、無ければ通常 (左=controls / 右=view)。
# subs: list(list(label=, id=, marker=), ...)。shared_controls 指定時は左パネルを共有し
# サブタブは右ビューのみ切替。無指定時はサブパネルを丸ごと切替。
pic_tab_panel <- function(sec_html, active = FALSE, controls = "", subs = NULL, shared_controls = "") {
  x <- pic_extract_section(sec_html)
  if (is.null(subs)) {
    n <- pic_pull_notes(x$inner)
    info <- info_badge(n$desc)
    p <- pic_pull_src(n$html)
    # ヘッダ右のアクション (例: correlation の PNG ボタン) を source の隣へ移す
    acts <- unique(unlist(regmatches(p$html, gregexpr('<span class="pic-headact">.*?</span>', p$html, perl = TRUE))))
    inner2 <- gsub('<span class="pic-headact">.*?</span>', '', p$html, perl = TRUE)
    acts_html <- gsub('</?span[^>]*>', '', paste(acts, collapse = ""))  # ラッパ span を除去
    head_right <- if (nzchar(p$src) || nzchar(acts_html))
      sprintf('<div class="pic-head-right">%s%s</div>', p$src, acts_html) else ""
    pane_cls <- if (nzchar(controls)) "pic-2pane" else "pic-2pane pic-1pane"
    ctrl <- if (nzchar(controls)) sprintf('<aside class="pic-ctrl">%s</aside>', controls) else ""
    return(sprintf('<section class="pic-tab%s" id="%s"><div class="pic-tab-head"><h2>%s</h2>%s%s</div><div class="%s">%s<div class="pic-view">%s</div></div></section>',
                   if (active) " active" else "", x$id, x$title, info, head_right, pane_cls, ctrl, inner2))
  }
  # サブタブは各サブパネルの左コントロール内 (View 見出し + 縦並びボタン) で選択する。
  chunks <- pic_split_inner(x$inner, vapply(subs, function(s) s$marker, character(1)))
  # 各サブタブの説明文を抽出し、どの View にかかるかを明記して info バッジに集約
  desc_parts <- character(0)
  for (i in seq_along(subs)) {
    pn <- pic_pull_notes(chunks[[i]])
    chunks[[i]] <- pn$html
    if (nzchar(pn$desc))
      desc_parts <- c(desc_parts, sprintf('<p class="pic-pop-sub">%s</p>%s', html_escape(subs[[i]]$label), pn$desc))
  }
  info <- info_badge(paste(desc_parts, collapse = ""))
  btns <- vapply(seq_along(subs), function(i) sprintf(
    '<button class="pic-subtabbtn%s" type="button" data-sub="%s">%s</button>',
    if (i == 1L) " active" else "", subs[[i]]$id, html_escape(subs[[i]]$label)), character(1))
  subnav <- sprintf('<nav class="pic-subtabs"><h4>View</h4>%s</nav>', paste(btns, collapse = ""))
  if (nzchar(shared_controls)) {
    # QC など: 左パネル (.pic-ctrl) にサブタブ + 共有コントロール、右にサブパネル
    panels <- vapply(seq_along(subs), function(i) sprintf(
      '<div class="pic-subpanel%s" id="%s">%s</div>', if (i == 1L) " active" else "", subs[[i]]$id, chunks[[i]]), character(1))
    body <- sprintf('<div class="pic-2pane"><aside class="pic-ctrl">%s%s</aside><div class="pic-view">%s</div></div>',
                    subnav, shared_controls, paste(panels, collapse = ""))
  } else {
    # DEG/Enrichment: サブパネル内のコントロール先頭 (<!--SUBNAV-->) にサブタブを差し込む
    panels <- vapply(seq_along(subs), function(i) {
      ck <- chunks[[i]]
      if (grepl("<!--SUBNAV-->", ck, fixed = TRUE)) {
        sp <- strsplit(ck, "<!--SUBNAV-->", fixed = TRUE)[[1]]
        ck <- paste0(sp[[1]], subnav, if (length(sp) >= 2) paste(sp[-1], collapse = "<!--SUBNAV-->") else "")
      } else {
        ck <- paste0(subnav, ck)
      }
      sprintf('<div class="pic-subpanel%s" id="%s">%s</div>', if (i == 1L) " active" else "", subs[[i]]$id, ck)
    }, character(1))
    body <- sprintf('<div class="pic-subbody">%s</div>', paste(panels, collapse = ""))
  }
  # ヘッダ: タイトル -> (info)。サブタブは左コントロールへ。
  sprintf('<section class="pic-tab%s" id="%s"><div class="pic-tab-head"><h2>%s</h2>%s</div>%s</section>',
          if (active) " active" else "", x$id, x$title, info, body)
}

# タブ html から source ファイルのパス一覧を取得する。
pic_tab_src_paths <- function(tab_html) {
  paths <- unique(regmatches(tab_html, gregexpr('(?<=<a class="pic-dlcsv" href=")[^"]+', tab_html, perl = TRUE))[[1]])
  if (length(paths) <= 1) return(paths)
  # 同一ディレクトリに多数 (>=3) のファイルがある場合はフォルダにまとめる (GSEA/ORA の大量 CSV 対策)
  collapse_once <- function(ps) {
    d <- dirname(ps)
    cnt <- table(d)
    unique(vapply(seq_along(ps), function(i)
      if (d[[i]] != "." && cnt[[d[[i]]]] >= 3L) d[[i]] else ps[[i]], character(1)))
  }
  repeat {
    nxt <- collapse_once(paths)
    if (identical(sort(nxt), sort(paths))) break
    paths <- nxt
  }
  paths
}

# Overview セクション。解析パイプラインを Materials & Methods 風のフローチャート
# (各ステップの使用ツール + パラメータ + 再現用コード) で提示する。
# summary/analysis_params.tsv (mapping / deseq2 が書く) を key -> value の list で返す。
# 無い場合 (旧バージョンの出力) は空 list を返し、従来の既定表記にフォールバックする。
pic_read_analysis_params <- function(out_dir) {
  f <- file.path(out_dir, "summary", "analysis_params.tsv")
  if (!file.exists(f)) return(list())
  tbl <- tryCatch(utils::read.delim(f, stringsAsFactors = FALSE, colClasses = "character"),
                  error = function(e) NULL)
  if (is.null(tbl) || !all(c("key", "value") %in% names(tbl))) return(list())
  stats::setNames(as.list(tbl$value), tbl$key)
}

pic_param <- function(params, key, default = "") {
  v <- params[[key]]
  if (is.null(v) || is.na(v)) default else as.character(v)
}

build_overview_section <- function(run, genome, pre_steps = list(), extra_tools = character(0),
                                   params = list()) {
  code_block <- function(txt) sprintf('<pre class="pic-code">%s</pre>', txt)
  # ツールバッジ: name|version を 1 つの monospace バッジ "name version" にまとめる。
  tools_html <- function(ts) paste(vapply(ts, function(t) {
    p <- strsplit(t, "|", fixed = TRUE)[[1]]
    label <- if (length(p) >= 2) paste0(p[[1]], " ", p[[2]]) else p[[1]]
    sprintf('<span class="pic-tool">%s</span>', html_escape(label))
  }, character(1)), collapse = "")
  step <- function(n, name, desc, code)  # ツールはステップに書かず冒頭にまとめる
    sprintf(paste0('<div class="pic-flow-step"><div class="pic-flow-head"><span class="pic-flow-n">%s</span>',
                   '<h3>%s</h3></div><p>%s</p>%s</div>'),
            n, name, desc, code)
  arrow <- '<div class="pic-flow-arrow">&#9660;</div>'
  # 全ステップで使用するツール + バージョンを 1 箇所にまとめて表示 (extra_tools は先頭に)
  tools_block <- sprintf('<div class="pic-flow-tools"><span class="pic-flow-tools-lbl">Tools</span>%s</div>',
    tools_html(c(extra_tools, "Trim Galore|0.6.10", "HISAT2|2.2.1", "samtools|1.23.1", "featureCounts|2.1.1",
                 "UMI-tools|1.1.4", "DESeq2|1.46.0", "R|4.4.3")))

  # 実際に使った hisat2 オプションを再現する (未記録の旧出力では既定設定として表示)。
  hisat2_opts <- ""
  if (identical(pic_param(params, "hisat2_very_sensitive"), "1")) {
    hisat2_opts <- paste0(hisat2_opts, " --very-sensitive")
  }
  score_min <- pic_param(params, "hisat2_score_min")
  if (nzchar(score_min)) hisat2_opts <- paste0(hisat2_opts, " --score-min ", score_min)
  align_note <- if (nzchar(hisat2_opts)) {
    sprintf(' Alignment sensitivity was raised with <code>%s</code>.', html_escape(trimws(hisat2_opts)))
  } else {
    ' HISAT2 defaults (<code>--score-min L,0,-0.2</code>) were used.'
  }

  # 実際に採用したサイズ因子の方法を反映する。UMI カウントが極端に疎な場合、
  # poscounts は深度差を検出できないためライブラリサイズ正規化にフォールバックする。
  # genome 別キーを優先し、旧出力 (genome 非依存キー) にはフォールバックする。
  sf_method <- pic_param(params, paste0("size_factor_method__", genome),
                         pic_param(params, "size_factor_method", "poscounts"))
  if (identical(sf_method, "libsize")) {
    sf_label <- "library-size"
    sf_code <- paste0(
'lib &lt;- colSums(counts(dds))\n',
'sizeFactors(dds) &lt;- lib / exp(mean(log(lib)))\n',
'dds &lt;- DESeq(dds)\n')
  } else {
    sf_label <- "poscounts"
    sf_code <- 'dds &lt;- DESeq(dds, sfType = "poscounts")\n'
  }

  step_defs <- c(pre_steps, list(
    list(name = "Adapter Trimming &amp; Alignment",
      desc = paste0('Trim the 3&#39; PIC adapter (<code>Trim Galore</code>), then align single-end to the ',
                    '<code>HISAT2</code> index.', align_note),
      code = code_block(sprintf(paste0(
'trim_galore -a GATCGTCGGACT -o trim/ demux/${sample}.fastq.gz\n',
'hisat2 -x hisat2_index/%s%s -U trim/${sample}_trimmed.fq.gz -S map/${sample}.sam\n',
'samtools sort -o map/${sample}.bam map/${sample}.sam\n',
'samtools index map/${sample}.bam'), html_escape(genome), html_escape(hisat2_opts)))),

    list(name = "Read-to-Gene Assignment",
      desc = '<code>featureCounts</code> on the forward/sense strand (<code>-s 1</code>); the gene id is stored in each read&#39;s <code>XT</code> tag.',
      code = code_block(sprintf(paste0(
'featureCounts -s 1 -t exon -g gene_id -a %s.gtf -R BAM \\\n',
'    -o count/${sample}.fc.tsv map/${sample}.bam\n',
'samtools sort -o map/${sample}.assigned.bam map/${sample}.bam.featureCounts.bam\n',
'samtools index map/${sample}.assigned.bam'), html_escape(genome)))),

    list(name = "UMI Counting",
      desc = 'Collapse UMIs per gene (<code>XT</code>) per barcode with <code>UMI-tools</code> (exact-<code>unique</code> method).',
      code = code_block(paste0(
'umi_tools count --method=unique --per-gene --gene-tag=XT --per-cell \\\n',
'    -I map/${sample}.assigned.bam -S count/${sample}.umi.tsv'))),

    list(name = "Differential Expression",
      desc = paste0('Join the per-sample UMI counts into a genes&times;samples matrix and run <code>DESeq2</code>: ',
                    '<code>~group</code> design, <b>', sf_label, '</b> size factors, all pairwise group contrasts, ',
                    'DEGs at <code>padj&nbsp;&lt;&nbsp;0.1</code>.'),
      code = code_block(paste0(
'<span class="c"># R</span>\n',
'grp &lt;- c("Cntl_Nega","Cntl_Nega","Cntl_Nega", "Cntl_Posi","Cntl_Posi","Cntl_Posi", ...)  <span class="c"># one per sample</span>\n',
'coldata &lt;- data.frame(group = factor(grp), row.names = colnames(umi_counts))\n',
'dds &lt;- DESeqDataSetFromMatrix(umi_counts, coldata, design = ~ group)\n',
sf_code,
'<span class="c"># for each pairwise group contrast (A vs B):</span>\n',
'res &lt;- results(dds, contrast = c("group", "A", "B"),\n',
'               independentFiltering = FALSE, cooksCutoff = FALSE)'))))
  )
  steps <- vapply(seq_along(step_defs), function(i) {
    d <- step_defs[[i]]; step(as.character(i), d$name, d$desc, d$code)
  }, character(1))

  overview <- sprintf(paste0('<p class="pic-ov-intro">This self-contained report presents the <b>PIC</b> ',
    '(photo-isolation chemistry) 3&prime;-biased RNA-seq run <b>%s</b>, mapped to <b>%s</b>.</p>'),
    html_escape(run), html_escape(genome))
  pipe_intro <- '<p class="pic-flow-sub">Starting from the demultiplexed per-sample FASTQ. Run each command per sample.</p>'
  sprintf(paste0('<section id="overview"><h2>Overview</h2>%s',
                 '<h3 class="pic-flow-title">Analysis Pipeline</h3>%s%s<div class="pic-flow">%s</div></section>'),
          overview, pipe_intro, tools_block, paste(steps, collapse = arrow))
}

# Downloads セクション。プロジェクトの全 csv/tsv をカテゴリ別 (折りたたみ) に列挙。
# 各ファイルは HTML 埋め込み・クリックで DL。
build_downloads_section <- function(tabs, run, genome, group_pal = NULL, geno_map = NULL) {
  proj <- getOption("pic.report.projdir")
  if (!is.null(proj)) {
    all <- list.files(proj, pattern = "\\.(csv|tsv)$", recursive = TRUE)
    all <- all[!grepl("(^|/)enrich_old_backup/", all)]
    for (rel in all) register_file(rel)
  }
  reg <- getOption("pic.report.reg")
  files <- if (!is.null(reg) && !is.null(reg$files)) reg$files else list()
  proj_name <- paste0(genome, "_", run)

  dl_btn <- function(i, label) sprintf('<a class="pic-dlcsv pic-ov-dl" role="button" tabindex="0" data-file="%s">%s</a>', i, label)
  ext_label <- function(name) { ext <- toupper(tools::file_ext(name)); if (nzchar(ext)) sprintf("Download %s", ext) else "Download" }
  file_li <- function(i) {
    f <- files[[i]]
    sprintf('<li>%s<span class="pic-ov-desc">%s</span></li>',
            dl_btn(i, ext_label(f$name)), html_escape(f$desc))
  }

  # 2 欄グリッド (行優先) の並び: 左=Mapping/Aggregation/Count Table/GSEA、右=PCA/Correlation/Diff/ORA
  order_cats <- c("Mapping", "PCA", "Aggregation", "Sample Correlation",
                  "Count Table", "Differential Expression",
                  "Enrichment — GSEA", "Enrichment — ORA")
  # xenograft では "… — graft/host" 接尾辞が付く。base 部で並べ、画分は graft->host。
  ct_base <- function(ct) sub(" — (graft|host)$", "", ct)
  ct_frac <- function(ct) if (grepl(" — graft$", ct)) 1L else if (grepl(" — host$", ct)) 2L else 0L
  meth_of <- function(i) { segs <- strsplit(files[[i]]$path, "/", fixed = TRUE)[[1]]
    gi <- which(segs %in% c("GSEA", "ORA")); if (length(gi) && gi[[1]] < length(segs)) segs[[gi[[1]] + 1L]] else "?" }
  cats <- unique(vapply(files, function(f) f$cat, character(1)))
  cats <- cats[order(vapply(cats, ct_frac, integer(1)),
                     vapply(cats, function(ct) { r <- match(ct_base(ct), order_cats); if (is.na(r)) length(order_cats) + 1L else r }, integer(1)),
                     cats)]

  secs <- character(0); sec_cts <- character(0)
  for (ct in cats) {
    ids <- names(files)[vapply(files, function(f) identical(f$cat, ct), logical(1))]
    if (length(ids) == 0) next
    tab <- files[[ids[[1]]]]$tab
    open_this <- !(ct_base(ct) %in% c("Enrichment — GSEA", "Enrichment — ORA"))
    if (identical(ct_base(ct), "Enrichment — GSEA")) {
      meth_order <- strsplit(pic_plot_spec()$defaults$enrich_methods_csv, ",", fixed = TRUE)[[1]]
      order_meths <- function(ms) c(meth_order[meth_order %in% ms], sort(setdiff(ms, meth_order)))
      meth <- vapply(ids, meth_of, character(1))
      subs <- character(0)
      for (mth in order_meths(unique(meth))) {
        mids <- ids[meth == mth]
        clis <- vapply(mids, function(i) {
          cf <- sub(paste0("^GSEA_", mth, "_"), "", files[[i]]$name)
          cf <- sub(paste0("_", proj_name, ".csv"), "", cf, fixed = TRUE)
          sprintf('<li>%s<span class="pic-ov-fname">%s</span></li>', dl_btn(i, "Download CSV"), color_vs_contrast(cf, group_pal))
        }, character(1))
        subs <- c(subs, sprintf('<details class="pic-ov-sub"><summary>%s <span class="pic-ov-n">(%d)</span></summary><ul class="pic-ov-list">%s</ul></details>',
                                html_escape(mth), length(mids), paste(clis, collapse = "")))
      }
      body <- paste(subs, collapse = "")
    } else if (identical(ct_base(ct), "Enrichment — ORA")) {
      # ORA は method ごとに 1 本の結合 CSV。method 順に並べ、名前 + Download を列挙
      meth_order <- strsplit(pic_plot_spec()$defaults$enrich_methods_csv, ",", fixed = TRUE)[[1]]
      order_meths <- function(ms) c(meth_order[meth_order %in% ms], sort(setdiff(ms, meth_order)))
      meth <- vapply(ids, meth_of, character(1))
      lis <- unlist(lapply(order_meths(unique(meth)), function(mth) {
        vapply(ids[meth == mth], function(i)
          sprintf('<li>%s<span class="pic-ov-desc">%s</span></li>',
                  dl_btn(i, "Download CSV"), html_escape(mth)), character(1))
      }))
      body <- sprintf('<ul class="pic-ov-list">%s</ul>', paste(lis, collapse = ""))
    } else {
      ids <- ids[order(vapply(ids, function(i) files[[i]]$name, character(1)))]
      ids <- ids[!duplicated(vapply(ids, function(i) files[[i]]$name, character(1)))]  # 同名は 1 つに
      body <- sprintf('<ul class="pic-ov-list">%s</ul>', paste(vapply(ids, file_li, character(1)), collapse = ""))
    }
    jump <- if (!is.null(tab) && nzchar(tab)) sprintf('<a class="pic-ov-jump" data-target="%s" tabindex="0">open tab &rarr;</a>', tab) else ""
    # xenograft は画分名にゲノムを併記 (例 "Mapping — host (mm10)")
    ct_disp <- ct
    frac_of_ct <- if (grepl(" — graft$", ct)) "graft" else if (grepl(" — host$", ct)) "host" else ""
    if (nzchar(frac_of_ct) && !is.null(geno_map) && frac_of_ct %in% names(geno_map) && nzchar(geno_map[[frac_of_ct]]))
      ct_disp <- sprintf("%s (%s)", ct, geno_map[[frac_of_ct]])
    # GSEA/ORA 以外は展開、GSEA/ORA は畳む
    sec <- sprintf('<details class="pic-ov-sec"%s><summary><span class="pic-ov-cat">%s</span>%s</summary>%s</details>',
                   if (open_this) " open" else "", html_escape(ct_disp), jump, body)
    secs <- c(secs, sec); sec_cts <- c(sec_cts, ct)
  }
  fr <- vapply(sec_cts, ct_frac, integer(1))
  if (!any(fr > 0L)) {
    # 通常レポート: 2 欄グリッド
    return(sprintf('<section id="downloads"><h2>Downloads</h2><div class="pic-ov-groups pic-ov-2col">%s</div></section>',
                   paste(secs, collapse = "")))
  }
  # xenograft: Genome split を上部 (全幅) に、graft/host は 2 欄でゲノム切替
  split_i <- which(sec_cts == "Genome split")
  graft_i <- which(fr == 1L); host_i <- which(fr == 2L)
  top <- if (length(split_i)) sprintf('<div class="pic-dl-top">%s</div>', paste(secs[split_i], collapse = "")) else ""
  grid <- function(sel) sprintf('<div class="pic-ov-groups pic-ov-2col">%s</div>', paste(secs[sel], collapse = ""))
  graft_block <- sprintf('<div class="pic-geno-dl" data-geno="graft">%s</div>', grid(graft_i))
  host_block  <- sprintf('<div class="pic-geno-dl" data-geno="host" hidden>%s</div>', grid(host_i))
  sprintf('<section id="downloads"><h2>Downloads</h2>%s%s%s</section>', top, graft_block, host_block)
}

# タブバー + パネルを組み立てる。tabs: list(list(html=, label=, controls=, subs=, shared_controls=))。
# 各 tab は任意で target (ボタンの data-target。未指定なら section id)、active (既定は先頭)、
# prebuilt (TRUE なら html は pic_tab_panel 済みとして 1 ボタン複数セクションに使う) を持てる。
# extra_nav: タブボタンの後に付ける追加ナビ HTML。
build_tabs <- function(tabs, extra_nav = "") {
  tabs <- Filter(function(t) nzchar(t$html), tabs)
  if (length(tabs) == 0) return(list(bar = "", panels = ""))
  btns <- character(0); panels <- character(0)
  for (i in seq_along(tabs)) {
    t <- tabs[[i]]
    active <- if (!is.null(t$active)) isTRUE(t$active) else (i == 1L)
    target <- if (!is.null(t$target)) t$target else regmatches(t$html, regexpr('(?<=<section id=")[^"]+', t$html, perl = TRUE))
    btns <- c(btns, sprintf('<button class="pic-tabbtn%s" type="button" data-target="%s">%s</button>',
                            if (active) " active" else "", target, html_escape(t$label)))
    panels <- c(panels, if (isTRUE(t$prebuilt)) t$html else pic_tab_panel(t$html, active,
                                      controls = if (!is.null(t$controls)) t$controls else "",
                                      subs = t$subs,
                                      shared_controls = if (!is.null(t$shared_controls)) t$shared_controls else ""))
  }
  list(bar = paste0(paste(btns, collapse = ""), extra_nav),
       panels = paste(panels, collapse = "\n"))
}

# HTML ページを組み立てて書き出す (plotly を内包)。
render_report_page <- function(title, nav_html, body_html, reg, asset_dir, out_html, extra_nav = "") {
  plots_json <- jsonlite::toJSON(reg$plots, auto_unbox = TRUE, null = "null", na = "null", digits = 6)
  expr_json <- jsonlite::toJSON(if (length(reg$expr) > 0) reg$expr else stats::setNames(list(), character(0)),
                                auto_unbox = TRUE, null = "null", na = "null", digits = 6)
  # 埋め込みファイル (gzip+base64) — DL 用に {id:{n:name, gz:base64}} だけを出力
  files_min <- if (!is.null(reg$files)) lapply(reg$files, function(f) list(n = f$name, gz = f$gz)) else list()
  files_json <- jsonlite::toJSON(if (length(files_min) > 0) files_min else stats::setNames(list(), character(0)),
                                 auto_unbox = TRUE, null = "null", na = "null")
  plotly_js <- paste(readLines(file.path(asset_dir, "plotly.min.js"), warn = FALSE), collapse = "\n")
  page <- paste0(
    '<!DOCTYPE html><html lang="ja"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>pic report — ', html_escape(title), '</title>',
    '<style>', report_css(), '</style>',
    '</head><body>',
    # ナビバー左端にレポート名 (ただの文字)、区切り線、右側にタブを両端揃え
    '<nav class="pic-tabs"><span class="pic-brand">', html_escape(title), ' Report</span>',
    '<div class="pic-tabbtns">', nav_html, '</div>', extra_nav, '</nav>',
    '<main class="pic-main">', body_html, '</main>',
    '<script>', plotly_js, '</script>',
    '<script>var PIC_PLOTS=', plots_json, ';var PIC_EXPR=', expr_json, ';var PIC_FILES=', files_json, ';</script>',
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

# 分類サマリ TSV から「Genome split」タブ (2 ゲノムへの振り分け割合) を作る。
# Read distribution と同じ 2 ペイン設計 (左に group トグル、行は data-group + 色付きラベル)。
# geno_map: c(graft="hg38", host="mm10") — 凡例ラベル用。group_map/group_pal/sample_order で
# サンプルの色・並び・group トグル対応。
section_genome_split <- function(summary_df, src = "", geno_map = NULL,
                                 group_map = NULL, group_pal = NULL, sample_order = NULL) {
  cats <- names(PIC_XENO_COLORS)
  have <- all(cats %in% colnames(summary_df))
  cat_label <- function(cc) if (!is.null(geno_map) && cc %in% names(geno_map) && nzchar(geno_map[[cc]]))
    sprintf("%s (%s)", cc, geno_map[[cc]]) else cc
  scol <- function(s) sample_color(s, group_map, group_pal)
  grpof <- function(s) if (!is.null(group_map) && s %in% names(group_map)) group_map[[s]] else s
  note <- paste0('<p class="pic-note">Each library was split between the two reference genomes with ',
                 '<b>xengsort</b>. Each bar is one sample (100&thinsp;% stacked): <b>graft</b> / <b>host</b> are reads ',
                 'assigned uniquely to each genome; <b>both</b> map to either, <b>ambiguous</b> are undecided, ',
                 '<b>neither</b> map to no genome. A clean split is dominated by graft + host. Hover a segment for counts.</p>')
  parts <- c('<section id="genomesplit"><h2>Genome split</h2>', note, src_note(src))
  if (!have) {
    parts <- c(parts, '<p>classification summary not found.</p></section>')
    return(paste(parts, collapse = "\n"))
  }
  samples <- as.character(summary_df$sample)
  ord <- if (!is.null(sample_order)) order(match(samples, sample_order)) else seq_along(samples)
  ruler <- '<div class="pic-bar-ruler"><span>0%</span><span>20%</span><span>40%</span><span>60%</span><span>80%</span><span>100%</span></div>'
  rows <- character(0); grp_prev <- NULL
  for (i in ord) {
    sample <- samples[[i]]
    g_cur <- grpof(sample)
    sep_cls <- if (!is.null(grp_prev) && !identical(g_cur, grp_prev)) " pic-grp-start" else ""
    grp_prev <- g_cur
    vals <- vapply(cats, function(cc) suppressWarnings(as.numeric(summary_df[[cc]][[i]])), numeric(1))
    vals[!is.finite(vals)] <- 0
    tot <- sum(vals); if (tot <= 0) tot <- 1
    segs <- character(0)
    for (cc in cats) {
      pct <- 100 * vals[[cc]] / tot
      if (pct <= 0) next
      tip <- sprintf("%s: %s (%.1f%%)", cat_label(cc), fmt_int(vals[[cc]]), pct)
      segs <- c(segs, sprintf('<div class="pic-seg" style="width:%.4f%%;background:%s" title="%s — %s"></div>',
                              pct, PIC_XENO_COLORS[[cc]], html_escape(sample), html_escape(tip)))
    }
    rows <- c(rows, sprintf('<div class="pic-bar-row%s" data-group="%s"><div class="pic-bar-label" style="color:%s;font-weight:600">%s</div><div class="pic-bar-track">%s</div></div>',
                            sep_cls, html_escape(g_cur), scol(sample), html_escape(sample), paste(segs, collapse = "")))
  }
  legend <- vapply(cats, function(cc) sprintf('<span class="pic-legend-item"><span class="pic-swatch" style="background:%s"></span>%s</span>',
                                              PIC_XENO_COLORS[[cc]], html_escape(cat_label(cc))), character(1))
  parts <- c(parts,
    sprintf('<div class="pic-legend" style="margin-left:0;margin-bottom:6px">%s</div>', paste(legend, collapse = "")),
    '<div class="pic-bars">', ruler, paste(rows, collapse = "\n"), '</div>',
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

# 1 画分 (graft/host) の解析タブ (data_tabs) を id_prefix "<frac>__" 付きで構築する。
# 戻り値: list(key, genome, data_tabs) または NULL (DESeq2 出力なし)。
build_fraction_analysis <- function(reg, frac_key, out_dir, frac_dir, tmp_dir) {
  projs <- pic_report_discover_projects(frac_dir)
  if (length(projs) == 0) return(NULL)
  desc <- projs[[1]]
  fdr <- pic_plot_spec()$defaults$fdr
  id_prefix <- paste0(frac_key, "__")
  genome <- fraction_genome(out_dir, frac_key, desc)

  stats <- suppressMessages(readr::read_csv(desc$stats_csv, show_col_types = FALSE, progress = FALSE))
  stats <- as.data.frame(stats, check.names = FALSE)
  deg_counts <- resolve_deg_counts(desc, stats, fdr)

  sheet <- pic_read_sample_sheet(out_dir, frac_key)
  sample_order <- if (!is.null(sheet)) sheet$samples else NULL
  group_order  <- if (!is.null(sheet) && length(sheet$groups) > 0) sheet$groups else NULL

  msum <- NULL; group_map <- NULL; group_pal <- NULL
  msum_files <- pic_list_summary(frac_dir, "^mapping_sum__.*\\.tsv$")
  if (length(msum_files) > 0) {
    msum <- pic_reorder_rows(read_mapping_sum(msum_files[[1]]), "sample", sample_order)
    mp <- group_map_pal(msum, group_order); group_map <- mp$map; group_pal <- mp$pal
  }

  data_tabs <- build_analysis_data_tabs(reg, desc, frac_dir, out_dir, id_prefix, fdr,
                                        msum, group_map, group_pal, sample_order, group_order,
                                        stats, deg_counts, tmp_dir)
  list(key = frac_key, genome = genome, data_tabs = data_tabs,
       group_pal = group_pal, group_map = group_map, sample_order = sample_order)
}

# out_dir 配下に xenograft 分類結果 (classify summary + graft/host) があるか
pic_is_xenograft_out <- function(out_dir) {
  s <- list.files(out_dir, pattern = "^xenograft_classify_summary__.*\\.tsv$", full.names = TRUE)
  length(s) > 0 && (dir.exists(file.path(out_dir, "graft")) || dir.exists(file.path(out_dir, "host")))
}

# xenograft レポート用の Overview。通常 Overview (Analysis Pipeline 含む) を再利用し、
# 冒頭の intro のみ xenograft 用 (2 ゲノム + xengsort split + ゲノム切替の案内) に差し替える。
build_overview_section_xeno <- function(run, fr_list, params = list()) {
  base_genome <- if (length(fr_list) > 0) fr_list[[1]]$genome else ""
  # index 名は hisat2 と同様に模式的に <graft>_on_<host> で示す
  graft_g <- if (!is.null(fr_list[["graft"]])) fr_list[["graft"]]$genome else "graft"
  host_g  <- if (!is.null(fr_list[["host"]]))  fr_list[["host"]]$genome  else "host"
  idx_name <- sprintf("%s_on_%s", graft_g, host_g)
  # パイプライン先頭に xengsort によるゲノム振り分けステップを追加 (実コマンド)
  xengsort_step <- list(
    name = "Genome split",
    desc = 'Classify each demultiplexed read to the <b>graft</b> or <b>host</b> genome with <code>xengsort</code>; each fraction is then mapped to its own reference by the steps below (run per genome).',
    code = paste0('<pre class="pic-code">xengsort classify --index ', html_escape(idx_name), ' --fastq demux/${sample}.fastq.gz \\\n',
                  '    -o xengsort/${sample} --mode count\n',
                  '<span class="c"># writes xengsort/${sample}-{graft,host,both,neither,ambiguous}.fq.gz</span></pre>'))
  base <- build_overview_section(run, base_genome, pre_steps = list(xengsort_step),
                                 extra_tools = "xengsort|2.2.0", params = params)  # <section><h2>Overview</h2><intro><h3>Analysis Pipeline</h3>...
  xeno_intro <- sprintf(paste0('<p class="pic-ov-intro">This self-contained report presents the <b>PIC</b> ',
    '(photo-isolation chemistry) 3&prime;-biased RNA-seq run <b>%s</b>, a <b>xenograft</b> library.</p>'),
    html_escape(run))
  # <h2>Overview</h2> と <h3 class="pic-flow-title"> の間 (= 元の intro) を差し替え
  marker <- '<h3 class="pic-flow-title">'
  segs <- strsplit(base, marker, fixed = TRUE)[[1]]
  paste0('<section id="overview"><h2>Overview</h2>', xeno_intro, marker, paste(segs[-1], collapse = marker))
}

build_xenograft_report <- function(out_dir, asset_dir) {
  options(pic.report.asset_dir = asset_dir)
  reg <- pic_report_registry()
  # Downloads の埋め込み・相対パス基準
  options(pic.report.reg = reg, pic.report.projdir = out_dir)

  summary_files <- list.files(out_dir, pattern = "^xenograft_classify_summary__.*\\.tsv$", full.names = TRUE)
  run <- sub("\\.tsv$", "", sub("^xenograft_classify_summary__", "", basename(summary_files[[1]])))
  summary_df <- suppressMessages(readr::read_tsv(summary_files[[1]], show_col_types = FALSE, progress = FALSE))
  summary_df <- as.data.frame(summary_df, check.names = FALSE)

  tmp_dir <- file.path(tempdir(), paste0("picxenoreport_", run))

  # 各画分の解析タブを構築 (graft を先に)
  fr_list <- list()
  for (frac in c("graft", "host")) {
    if (!dir.exists(file.path(out_dir, frac))) next
    tf <- file.path(tmp_dir, frac); dir.create(tf, recursive = TRUE, showWarnings = FALSE)
    fa <- build_fraction_analysis(reg, frac, out_dir, file.path(out_dir, frac), tf)
    if (!is.null(fa)) fr_list[[frac]] <- fa
  }
  first_geno <- if (length(fr_list) > 0) names(fr_list)[[1]] else "graft"
  geno_map <- vapply(fr_list, function(fa) fa$genome, character(1))  # c(graft="hg38", host="mm10")

  keys <- c("qc", "aggregate", "cor", "pca", "deg", "expr", "enrich")

  # --- 通常レポートと同じ build_tabs 経路でタブを組み立てる ---
  # 各 tab は prebuilt (pic_tab_panel 済み) + target (論理キー) で 1 ボタンに対応。
  # 解析タブは graft/host の 2 セクションを 1 タブに束ね、JS のゲノム切替で出し分ける。
  tabs <- list(list(target = "overview", label = "Overview", active = TRUE, prebuilt = TRUE,
                    html = pic_tab_panel(build_overview_section_xeno(run, fr_list, pic_read_analysis_params(out_dir)), active = TRUE)))

  # Genome split (Read distribution と同じ 2 ペイン: 左に group トグル)
  gsplit_pal <- if (length(fr_list) > 0) fr_list[[1]]$group_pal else NULL
  gsplit_map <- if (length(fr_list) > 0) fr_list[[1]]$group_map else NULL
  gsplit_ord <- if (length(fr_list) > 0) fr_list[[1]]$sample_order else NULL
  gs <- section_genome_split(summary_df, report_rel_path(summary_files[[1]], out_dir), geno_map,
                             gsplit_map, gsplit_pal, gsplit_ord)
  tabs <- c(tabs, list(list(target = "genomesplit", label = "Genome split", prebuilt = TRUE,
                            html = pic_tab_panel(gs, active = FALSE, controls = group_toggle_panel(gsplit_pal)))))

  # 解析タブ: 論理キーごとに両画分の <section id="<frac>__<key>"> を 1 タブに束ねる
  for (ki in seq_along(keys)) {
    k <- keys[[ki]]; sec_html <- character(0); label <- NULL
    for (frac in names(fr_list)) {
      t <- fr_list[[frac]]$data_tabs[[ki]]
      if (is.null(t) || !nzchar(t$html)) next
      if (is.null(label)) label <- t$label
      sec_html <- c(sec_html, pic_tab_panel(t$html, active = FALSE,
                                            controls = if (!is.null(t$controls)) t$controls else "",
                                            subs = t$subs,
                                            shared_controls = if (!is.null(t$shared_controls)) t$shared_controls else ""))
    }
    if (length(sec_html) > 0)
      tabs <- c(tabs, list(list(target = k, label = label, prebuilt = TRUE, html = paste(sec_html, collapse = ""))))
  }

  # Downloads (out_dir を再帰スキャンして graft/host 両方を収録)
  merged_tabs <- if (length(fr_list) > 0) fr_list[[1]]$data_tabs else list()
  dl <- build_downloads_section(merged_tabs, run, if (length(fr_list) > 0) fr_list[[1]]$genome else "", NULL, geno_map)
  tabs <- c(tabs, list(list(target = "downloads", label = "Downloads", prebuilt = TRUE,
                            html = pic_tab_panel(dl, active = FALSE))))

  # ゲノム切替コントロール (navbar 右端 = render の extra_nav)
  geno_btns <- vapply(names(fr_list), function(frac) {
    fa <- fr_list[[frac]]
    lab <- sprintf("%s%s", frac, if (nzchar(fa$genome)) sprintf(" (%s)", fa$genome) else "")
    sprintf('<button class="pic-genobtn%s" type="button" data-geno="%s">%s</button>',
            if (frac == first_geno) " active" else "", html_escape(frac), html_escape(lab))
  }, character(1))
  geno_switch <- if (length(geno_btns) > 0) sprintf('<div class="pic-genoswitch">%s</div>', paste(geno_btns, collapse = "")) else ""

  tb <- build_tabs(tabs)
  out_html <- file.path(out_dir, sprintf("report_%s.html", run))
  render_report_page(run, tb$bar, tb$panels, reg, asset_dir, out_html, extra_nav = geno_switch)
  unlink(tmp_dir, recursive = TRUE)
  out_html
}

pic_load_deg_counts <- function(deseq2_dir, project) {
  f <- file.path(deg_dir_of(deseq2_dir, project), sprintf("DEG_count_%s.csv", project))
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

