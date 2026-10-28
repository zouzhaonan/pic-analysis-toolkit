# 役割: QC タブの「RNA integrity」サブタブ (TapeStation の electropherogram と、各波形の隣の RINe)。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が
#       report_build.R (ローダ) 経由で source する。単体では動作しない。
#
# データ: TapeStation の export (<stem>_sampleTable.csv / _compactPeakTable.csv /
#         _Electropherogram.csv)。RIN データセットのフォルダ直下で、名前 (フォルダ名
#         またはファイル名) がプロジェクト ID (= run 名) で始まるものを使う。
#         見つからなければサブタブを出さない。
# サンプル名は sampleTable の Sample Description に従う (sample_sheet とは突き合わせない)。
# ladder は electropherogram の横軸 (データ点 -> nt) の較正だけに使い、描画しない。
# テープが普通の RNA ScreenTape (assay 名 "RNA"。HS RNA は "HSRNA") のときは、
# ポジティブコントロール (pic_rin_pc_dir() の泳動) を毎回、本番サンプルの下の別パネルに出す
# (y 軸は別に autoscale)。HS RNA には入れない。
# 1 run に複数 genome があるときは、レポートの genome のサンプルに当たる測定だけ出す
# (pic_rin_filter_genome)。

PIC_RIN_DEFAULT_DIR <- "~/OkiLab Dropbox/OkiLab_Data/TapeStation/RIN_raw_data"

# RIN データセットのフォルダ。優先: pic report --rin-dir → 環境変数 PIC_RIN_DIR → 既定パス。
pic_rin_dir <- function() {
  d <- getOption("pic.report.rin_dir", default = NULL)
  if (is.null(d)) d <- Sys.getenv("PIC_RIN_DIR", unset = "")
  if (!nzchar(d)) d <- PIC_RIN_DEFAULT_DIR
  d <- path.expand(d)
  d
}

PIC_RIN_PC_COLOUR <- "#555555"

# ポジティブコントロール (RNA ScreenTape 用) のフォルダ。
# 優先: pic report --rin-pc-dir → 環境変数 PIC_RIN_PC_DIR → <rin_dir>/RNA_positive_control
#       (あれば) → <既定パス>/RNA_positive_control。
pic_rin_pc_dir <- function(rin_dir = pic_rin_dir()) {
  d <- getOption("pic.report.rin_pc_dir", default = NULL)
  if (is.null(d)) d <- Sys.getenv("PIC_RIN_PC_DIR", unset = "")
  if (nzchar(d)) return(path.expand(d))
  d <- file.path(rin_dir, "RNA_positive_control")
  if (dir.exists(d)) return(d)
  file.path(path.expand(PIC_RIN_DEFAULT_DIR), "RNA_positive_control")
}

# CSV の一覧を TapeStation の 1 回の泳動 (stem) ごとに 3 種のファイルへまとめる。
pic_rin_group_files <- function(files) {
  if (length(files) == 0) return(list())
  kinds <- c(sample = "sampleTable", peak = "compactPeakTable", ep = "Electropherogram")
  re <- sprintf("_(%s)\\.csv$", paste(kinds, collapse = "|"))
  files <- files[grepl(re, basename(files)) & file.info(files)$size > 0]
  stems <- sub(re, "", files)
  out <- list()
  for (st in sort(unique(stems))) {
    fs <- files[stems == st]
    ds <- list(stem = basename(st))
    for (k in names(kinds)) {
      f <- fs[grepl(sprintf("_%s\\.csv$", kinds[[k]]), basename(fs))]
      ds[[k]] <- if (length(f) > 0) f[[1]] else NULL
    }
    if (!is.null(ds$sample)) out[[length(out) + 1]] <- ds
  }
  out
}

# ポジティブコントロールの泳動 (pc_dir 以下の CSV すべて)。
pic_rin_find_pc <- function(pc_dir) {
  if (!dir.exists(pc_dir)) return(list())
  pic_rin_group_files(list.files(pc_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE))
}

# テープの種類 (assay 名)。compactPeakTable の FileName ("2026-09-02 - 15.27.32.RNA") の
# 最後の "." の後ろ、無ければ stem の最後の "-" の後ろ。"RNA" / "HSRNA" など (大文字)。
pic_rin_assay <- function(ds) {
  pk <- pic_rin_read_csv(ds$peak)
  if (!is.null(pk) && nrow(pk) > 0) {
    c_fn <- pic_rin_col(pk, "^FileName$")
    if (!is.null(c_fn) && nzchar(pk[[c_fn]][[1]])) return(toupper(sub("^.*\\.", "", pk[[c_fn]][[1]])))
  }
  toupper(sub("^.*-", "", ds$stem))
}

# rin_dir 直下で名前が project_id で始まるもの (直後が英数字でないこと: pic19 は pic192 に
# 当たらない) を集め、TapeStation の 1 回の泳動 (stem) ごとに 3 種のファイルをまとめる。
# 複数当たれば全部使う (泳動ごとに 1 ブロック)。sampleTable が無い / 空の stem は捨てる。
pic_rin_find_datasets <- function(rin_dir, project_id) {
  if (is.null(project_id) || !nzchar(project_id) || !dir.exists(rin_dir)) return(list())
  hits <- list.files(rin_dir, full.names = TRUE)
  b <- tolower(basename(hits)); id <- tolower(project_id)
  nxt <- substr(b, nchar(id) + 1L, nchar(id) + 1L)
  hits <- hits[startsWith(b, id) & !grepl("[a-z0-9]", nxt)]
  files <- unlist(lapply(hits, function(h)
    if (dir.exists(h)) list.files(h, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE) else h))
  pic_rin_group_files(files)
}

# TapeStation の CSV は latin1 (µ が 0xB5)。全列を文字列で読み、読めなければ NULL。
pic_rin_read_csv <- function(path, header = TRUE) {
  if (is.null(path)) return(NULL)
  tryCatch(utils::read.csv(path, header = header, check.names = FALSE, stringsAsFactors = FALSE,
                           colClasses = "character", fileEncoding = "latin1", strip.white = TRUE),
           error = function(e) NULL)
}

# 列名を正規表現で引く (単位表記の揺れや µ の文字化けに左右されないように)。
pic_rin_col <- function(d, re) {
  hit <- grep(re, colnames(d), ignore.case = TRUE, value = TRUE)
  if (length(hit) > 0) hit[[1]] else NULL
}

pic_rin_num <- function(x) suppressWarnings(as.numeric(x))

# ladder の行: Observations が Ladder か、Well が EL<n> (Electronic Ladder)。
pic_rin_is_ladder <- function(well, obs) {
  grepl("ladder", obs, ignore.case = TRUE) | grepl("^EL[0-9]+$", well)
}

# 局所最大のうち高い順に n 個 (近接 (min_gap 点以内) の重複は除く) を位置順で返す。
pic_rin_top_peaks <- function(y, n, min_gap = 5L) {
  m <- length(y)
  if (m < 3 || n < 1) return(integer(0))
  i <- 2:(m - 1)
  cand <- i[y[i] > y[i - 1] & y[i] >= y[i + 1]]
  cand <- cand[order(-y[cand])]
  keep <- integer(0)
  for (c in cand) {
    if (all(abs(keep - c) > min_gap)) keep <- c(keep, c)
    if (length(keep) == n) break
  }
  sort(keep)
}

# 区分線形 (端は端の区間の傾きで外挿)。
pic_rin_interp <- function(x, xs, ys) {
  out <- stats::approx(xs, ys, xout = x, rule = 2)$y
  k <- length(xs)
  lo <- x < xs[[1]]; hi <- x > xs[[k]]
  out[lo] <- ys[[1]] + (x[lo] - xs[[1]]) * (ys[[2]] - ys[[1]]) / (xs[[2]] - xs[[1]])
  out[hi] <- ys[[k]] + (x[hi] - xs[[k]]) * (ys[[k]] - ys[[k - 1]]) / (xs[[k]] - xs[[k - 1]])
  out
}

# 1 泳動分を読み、samples (ladder 除く) と、較正済みの trace / ピーク注釈を返す。
pic_rin_load <- function(ds) {
  st <- pic_rin_read_csv(ds$sample)
  if (is.null(st) || nrow(st) == 0) return(NULL)
  c_well <- pic_rin_col(st, "^Well$"); c_desc <- pic_rin_col(st, "Sample Description")
  c_obs <- pic_rin_col(st, "^Observations$")
  if (is.null(c_well)) return(NULL)
  well <- st[[c_well]]
  obs <- if (!is.null(c_obs)) st[[c_obs]] else rep("", nrow(st))
  lad <- pic_rin_is_ladder(well, obs)
  c_rin <- pic_rin_col(st, "^RINe?$")
  samples <- data.frame(
    well = well, sample = if (!is.null(c_desc)) st[[c_desc]] else well,
    rine = if (!is.null(c_rin)) pic_rin_num(st[[c_rin]]) else NA_real_,
    stringsAsFactors = FALSE)[!lad, , drop = FALSE]
  samples$sample[!nzchar(samples$sample)] <- samples$well[!nzchar(samples$sample)]
  if (nrow(samples) == 0) return(NULL)
  ladder_wells <- well[lad]

  # compactPeakTable: well -> (size, 注釈名)
  pk <- pic_rin_read_csv(ds$peak)
  peaks <- NULL
  if (!is.null(pk) && nrow(pk) > 0) {
    c_pw <- pic_rin_col(pk, "^Well$"); c_sz <- pic_rin_col(pk, "^Size"); c_po <- pic_rin_col(pk, "^Observations$")
    if (!is.null(c_pw) && !is.null(c_sz))
      peaks <- data.frame(well = pk[[c_pw]], size = pic_rin_num(pk[[c_sz]]),
                          label = if (!is.null(c_po)) pk[[c_po]] else "", stringsAsFactors = FALSE)
  }

  # electropherogram: 1 列 = 1 well ("A1: name")、行 = データ点。横軸の値は持たない。
  traces <- list()
  lsz <- numeric(0)
  ep <- pic_rin_read_csv(ds$ep, header = FALSE)
  if (!is.null(ep) && nrow(ep) > 1 && !is.null(peaks)) {
    ep_well <- trimws(sub(":.*$", "", as.character(unlist(ep[1, ]))))
    ep <- ep[-1, , drop = FALSE]
    lw <- intersect(ladder_wells, ep_well)
    lsz <- if (length(lw) > 0) sort(peaks$size[peaks$well == lw[[1]] & is.finite(peaks$size)]) else numeric(0)
    if (length(lsz) >= 2) {
      ly <- pic_rin_num(ep[[match(lw[[1]], ep_well)]])
      lpos <- pic_rin_top_peaks(ly, length(lsz))
      if (length(lpos) == length(lsz)) {
        lm_pos <- lpos[[1]]
        for (i in seq_len(nrow(samples))) {
          j <- match(samples$well[[i]], ep_well)
          if (is.na(j)) next
          y <- pic_rin_num(ep[[j]])
          # 泳動ごとの位置ずれは lower marker (先頭ピーク) で合わせる
          win <- max(1L, lm_pos - 15L):min(length(y), lm_pos + 15L)
          shift <- win[[which.max(y[win])]] - lm_pos
          size <- 10^pic_rin_interp(seq_along(y) - shift, lpos, log10(lsz))
          keep <- is.finite(y) & size >= lsz[[1]] * 0.8 & size <= lsz[[length(lsz)]] * 1.25
          sp <- peaks[peaks$well == samples$well[[i]] & is.finite(peaks$size) & nzchar(peaks$label), , drop = FALSE]
          # 注釈の高さ: ピーク位置付近 (±3 点) の最大値
          sp$y <- vapply(sp$size, function(s) {
            k <- which.min(abs(log10(size) - log10(s)))
            max(y[max(1L, k - 3L):min(length(y), k + 3L)], na.rm = TRUE)
          }, numeric(1))
          traces[[samples$well[[i]]]] <- list(size = size[keep], y = y[keep], peaks = sp)
        }
      }
    }
    if (length(traces) == 0) message(sprintf("[WARN] RIN: electropherogram の ladder 較正に失敗しました: %s", ds$stem))
  }
  list(stem = ds$stem, samples = samples, traces = traces, ladder_sizes = lsz)
}

# ポジティブコントロール (pc_dir の泳動すべて) を 1 つの泳動のように読む (well は "PC<n>:<well>")。
# 名前は "Positive control" (PC が複数あれば " (<sample>)" を付ける)。無ければ NULL。
pic_rin_load_pc <- function(pc_dir) {
  pls <- Filter(Negate(is.null), lapply(pic_rin_find_pc(pc_dir), pic_rin_load))
  if (length(pls) == 0) return(NULL)
  n_pc <- sum(vapply(pls, function(p) nrow(p$samples), integer(1)))
  out <- list(stem = "positive control", samples = NULL, traces = list(), ladder_sizes = pls[[1]]$ladder_sizes)
  for (j in seq_along(pls)) {
    p <- pls[[j]]; s <- p$samples
    new_well <- sprintf("PC%d:%s", j, s$well)
    for (i in seq_len(nrow(s))) {
      tr <- p$traces[[s$well[[i]]]]
      if (!is.null(tr)) { tr$peaks$well <- rep(new_well[[i]], nrow(tr$peaks)); out$traces[[new_well[[i]]]] <- tr }
    }
    s$sample <- if (n_pc > 1) sprintf("Positive control (%s)", s$sample) else "Positive control"
    s$well <- new_well
    out$samples <- rbind(out$samples, s)
  }
  out
}

# 1 run に複数 genome があるとき、このレポート (project) の genome のサンプルに当たる測定だけ残す。
# 照合: summary/deftable_<run>_<genome>.tsv の sample / group 名と TapeStation の Sample Description
# (大小文字を無視)。一致、または名前を区切り (_ - . 空白) で割った語の並びの先頭・途中に現れれば当たり
# (例: "tora" は "tora_posi" に当たる)。どの genome にも当たらない測定は全レポートに残す。
# genome が 1 つ (deftable が 1 つ以下) なら絞らない。
pic_rin_genome_names <- function(base_dir, run) {
  if (is.null(base_dir) || !nzchar(run)) return(list())
  dfs <- pic_list_summary(base_dir, "^deftable_.*\\.tsv$")
  dfs <- dfs[startsWith(basename(dfs), paste0("deftable_", run, "_"))]
  out <- list()
  for (f in dfs) {
    d <- tryCatch(utils::read.delim(f, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character"),
                  error = function(e) NULL)
    if (is.null(d)) next
    nm <- unique(c(d$sample, d$group)); nm <- nm[!is.na(nm) & nzchar(nm)]
    out[[sub("\\.tsv$", "", sub("^deftable_", "", basename(f)))]] <- tolower(nm)
  }
  out
}

pic_rin_name_hits <- function(ts, names) {
  t <- tolower(trimws(ts))
  if (!nzchar(t)) return(FALSE)
  tt <- strsplit(t, "[_. -]+")[[1]]
  any(vapply(names, function(n) {
    if (n == t) return(TRUE)
    tok <- strsplit(n, "[_. -]+")[[1]]
    k <- length(tt)
    if (k > length(tok)) return(FALSE)
    any(vapply(seq_len(length(tok) - k + 1L), function(i) all(tok[i:(i + k - 1L)] == tt), logical(1)))
  }, logical(1)))
}

pic_rin_filter_genome <- function(r, project, gnames) {
  if (is.null(project) || length(gnames) < 2 || !(project %in% names(gnames))) return(r)
  keep <- vapply(r$samples$sample, function(ts) {
    hit <- names(gnames)[vapply(gnames, function(n) pic_rin_name_hits(ts, n), logical(1))]
    length(hit) == 0 || project %in% hit
  }, logical(1))
  r$samples <- r$samples[keep, , drop = FALSE]
  r$traces <- r$traces[intersect(names(r$traces), r$samples$well)]
  r
}

# RINe ラベルの縦位置: 望む位置 (want) をなるべく保ち、隣同士を gap 以上離す。
# 下から押し上げ、上限 hi を超えたら上から押し下げる (順序は want の順のまま)。
pic_rin_spread <- function(want, gap, lo, hi) {
  n <- length(want)
  if (n == 0) return(numeric(0))
  o <- order(want); y <- want[o]
  y[[1]] <- max(y[[1]], lo)
  if (n > 1) for (i in 2:n) y[[i]] <- max(y[[i]], y[[i - 1]] + gap)
  if (y[[n]] > hi) {
    y[[n]] <- hi
    if (n > 1) for (i in (n - 1):1) y[[i]] <- min(y[[i]], y[[i + 1]] - gap)
  }
  out <- numeric(n); out[o] <- y
  out
}

# 波形のどこから RINe ラベルへ線を引くか: 28S ピークの右肩 (28S の注釈文字を避ける)。
# 28S が無ければ (分解した RNA など) 200 nt 以上で最も高い点の右肩。
pic_rin_anchor <- function(tr) {
  p28 <- tr$peaks$size[grepl("^28S", tr$peaks$label)]
  big <- which(tr$size >= 200)
  if (length(big) == 0) big <- seq_along(tr$size)
  x0 <- if (length(p28) > 0) p28[[1]] else tr$size[[big[[which.max(tr$y[big])]]]]
  k <- which.min(abs(log10(tr$size) - log10(x0 * 1.2)))
  c(x = tr$size[[k]], y = tr$y[[k]])
}

# 1 パネル分の plotly 図 (electropherogram + RINe ラベル)。x_lim は log10 前の [min, max] (パネル間で揃える)。
# y 軸はパネルごとに autoscale。図を作れなければ ""。
# colour を渡すと全サンプルをその色で描く (PC 用)。legend = FALSE で凡例を出さない。
pic_rin_panel <- function(reg, id, r, x_lim, colour = NULL, legend = TRUE) {
  if (is.null(r) || length(r$traces) == 0) return("")
  s <- r$samples[r$samples$well %in% names(r$traces), , drop = FALSE]
  if (nrow(s) == 0) return("")
  pal <- if (is.null(colour)) group_palette(s$sample) else stats::setNames(rep(colour, nrow(s)), s$sample)
  traces <- list()
  for (i in seq_len(nrow(s))) {
    tr <- r$traces[[s$well[[i]]]]
    nm <- s$sample[[i]]; col <- unname(pal[[nm]])
    traces[[length(traces) + 1]] <- list(
      x = as.list(round(tr$size, 1)), y = as.list(round(tr$y, 3)),
      name = sprintf('<span style="color:%s">%s</span>', col, html_escape(nm)), legendgroup = nm,
      showlegend = legend, mode = "lines", type = "scatter", line = list(color = col, width = 1.4),
      hovertemplate = sprintf("%s<br>%%{x:.0f} nt<br>%%{y:.1f} FU<extra></extra>", html_escape(nm)))
    if (nrow(tr$peaks) > 0) {
      traces[[length(traces) + 1]] <- list(
        x = as.list(tr$peaks$size), y = as.list(tr$peaks$y), text = as.list(html_escape(tr$peaks$label)),
        name = html_escape(nm), legendgroup = nm, showlegend = FALSE,
        mode = "markers+text", type = "scatter", textposition = "top center", cliponaxis = FALSE,
        textfont = list(size = 11, color = col),
        marker = list(color = col, size = 7, symbol = "triangle-down"),
        customdata = as.list(tr$peaks$size),
        hovertemplate = sprintf("%s<br>%%{text}: %%{customdata} nt<extra></extra>", html_escape(nm)))
    }
  }
  # RINe ラベル: 波形の右に 1 列に並べ (重ならないよう縦にずらす)、細い点線で波形とつなぐ。
  # 凡例のクリックで波形と一緒に隠れるよう、同じ legendgroup に入れる。
  trs <- r$traces[s$well]
  lab <- s[is.finite(s$rine), , drop = FALSE]
  ymax <- max(vapply(trs, function(t) max(c(t$y, t$peaks$y), na.rm = TRUE), numeric(1)))
  n_leg <- nrow(s)
  # 縦の大きさ: 凡例の行 (5 サンプル/行ほど) を除いた描画域に、ラベル 1 つ 17 px が収まるように
  leg_px <- if (legend) 22 * ceiling(n_leg / 5) else 0
  height <- max(420, ceiling(76 + leg_px + nrow(lab) * 17 / 0.9))
  xcol <- x_lim[[2]] * 1.08
  xrange <- list(log10(x_lim[[1]]), log10(xcol * 1.5))
  if (nrow(lab) > 0 && is.finite(ymax) && ymax > 0) {
    area <- height - 76 - leg_px
    anc <- t(vapply(lab$well, function(w) pic_rin_anchor(r$traces[[w]]), numeric(2)))
    ly <- pic_rin_spread(anc[, "y"], gap = ymax * 17 / area, lo = ymax * 0.03, hi = ymax * 0.97)
    for (i in seq_len(nrow(lab))) {
      nm <- lab$sample[[i]]; col <- unname(pal[[nm]])
      txt <- formatC(lab$rine[[i]], format = "f", digits = 1)
      traces[[length(traces) + 1]] <- list(
        x = list(anc[i, "x"], xcol), y = list(anc[i, "y"], ly[[i]]),
        text = list("", sprintf("<b>%s</b>", txt)), name = html_escape(nm), legendgroup = nm, showlegend = FALSE,
        mode = "lines+text", type = "scatter", textposition = "middle right", cliponaxis = FALSE,
        textfont = list(size = 12, color = col), line = list(color = col, width = 0.8, dash = "dot"),
        hovertemplate = sprintf("%s<br>RINe %s<extra></extra>", html_escape(nm), txt))
    }
  }
  layout <- list(
    xaxis = list(title = "size [nt]", type = "log", tickmode = "array", tickvals = as.list(r$ladder_sizes),
                 ticktext = as.list(format(r$ladder_sizes, trim = TRUE)), zeroline = FALSE, range = xrange),
    yaxis = list(title = "sample intensity [FU]", autorange = TRUE),
    hovermode = "closest", margin = list(l = 60, r = 20, t = if (legend) 30 else 16, b = 46),
    showlegend = legend,
    legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0))
  register_plot(reg, id, list(data = traces, layout = layout))
  sprintf('<div id="%s" class="pic-plot" style="height:%dpx"></div>', id, height)
}

# QC タブ「RNA integrity」の中身 (説明 + 本番サンプルのパネル + その下に PC のパネル)。
# データが無ければ ""。project / base_dir を渡すと、複数 genome の run で genome ごとに絞る。
build_rin_html <- function(reg, project_id, id_prefix = "", rin_dir = pic_rin_dir(),
                           project = NULL, base_dir = NULL) {
  dss <- pic_rin_find_datasets(rin_dir, project_id)
  if (length(dss) == 0) return("")
  gnames <- pic_rin_genome_names(base_dir, project_id)
  runs <- list(); rna_tape <- FALSE
  for (k in seq_along(dss)) {
    r <- pic_rin_load(dss[[k]])
    if (is.null(r)) next
    r <- pic_rin_filter_genome(r, project, gnames)
    if (nrow(r$samples) == 0 || length(r$traces) == 0) next
    if (identical(pic_rin_assay(dss[[k]]), "RNA")) rna_tape <- TRUE
    runs[[length(runs) + 1]] <- r
  }
  if (length(runs) == 0) return("")
  pc <- NULL
  if (rna_tape) {
    pc_dir <- pic_rin_pc_dir(rin_dir)
    pc <- pic_rin_load_pc(pc_dir)
    if (is.null(pc) || length(pc$traces) == 0) {
      message(sprintf("[WARN] RIN: RNA ScreenTape ですが、ポジティブコントロールが見つかりません: %s", pc_dir))
      pc <- NULL
    }
  }
  # x 軸はパネル間で揃える
  all_tr <- unlist(lapply(c(runs, if (!is.null(pc)) list(pc)), function(r) r$traces), recursive = FALSE)
  x_lim <- c(min(vapply(all_tr, function(t) min(t$size), numeric(1))),
             max(vapply(all_tr, function(t) max(t$size), numeric(1))))
  blocks <- character(0)
  for (k in seq_along(runs)) {
    plot <- pic_rin_panel(reg, sprintf("%srin_%d", id_prefix, k), runs[[k]], x_lim)
    if (!nzchar(plot)) next
    head <- if (length(runs) > 1) sprintf('<h4>%s</h4>', html_escape(runs[[k]]$stem)) else ""
    blocks <- c(blocks, sprintf('<div class="pic-plot-cell">%s%s</div>', head, plot))
  }
  if (length(blocks) == 0) return("")
  pc_used <- FALSE
  if (!is.null(pc)) {
    # PC はパネルの見出しで分かるので凡例は出さず、本番と区別できる暗いグレーで描く
    plot <- pic_rin_panel(reg, sprintf("%srin_pc", id_prefix), pc, x_lim, colour = PIC_RIN_PC_COLOUR, legend = FALSE)
    if (nzchar(plot)) {
      pc_used <- TRUE
      blocks <- c(blocks, sprintf('<div class="pic-plot-cell"><h4>Positive control</h4>%s</div>', plot))
    }
  }
  paste0('<p class="pic-note">RNA integrity of the input RNA measured on the TapeStation. The electropherogram shows each sample&rsquo;s fluorescence along fragment size; the lower marker and the ',
         '18S / 28S ribosomal RNA peaks are marked. The number next to each trace, in the same colour, is its <b>RINe</b> ',
         '(1 = degraded, 10 = intact). Click a sample in the legend to show or hide it. ',
         'Sample names follow the TapeStation sample table.',
         if (pc_used) ' The lower panel is the <b>positive control</b>, a reference RNA of known good integrity run on a separate RNA ScreenTape, for comparison (its intensity axis is scaled separately).' else '',
         '</p><div class="pic-rin-panels" style="display:grid;grid-template-columns:1fr;gap:16px">',
         paste(blocks, collapse = ""), '</div>')
}
