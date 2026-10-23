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

PIC_RIN_DEFAULT_DIR <- "/Users/zou/OkiLab Dropbox/OkiLab_Data/TapeStation/RIN_raw_data"

# RIN データセットのフォルダ。優先: pic report --rin-dir → 環境変数 PIC_RIN_DIR → 既定パス。
pic_rin_dir <- function() {
  d <- getOption("pic.report.rin_dir", default = NULL)
  if (is.null(d)) d <- Sys.getenv("PIC_RIN_DIR", unset = "")
  if (!nzchar(d)) d <- PIC_RIN_DEFAULT_DIR
  d
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

# QC タブ「RNA integrity」の中身 (説明 + electropherogram + RINe ラベル)。データが無ければ ""。
build_rin_html <- function(reg, project_id, id_prefix = "", rin_dir = pic_rin_dir()) {
  dss <- pic_rin_find_datasets(rin_dir, project_id)
  if (length(dss) == 0) return("")
  blocks <- character(0)
  for (k in seq_along(dss)) {
    r <- pic_rin_load(dss[[k]])
    if (is.null(r)) next
    s <- r$samples
    pal <- group_palette(s$sample)
    plot <- ""
    if (length(r$traces) > 0) {
      traces <- list()
      for (i in seq_len(nrow(s))) {
        tr <- r$traces[[s$well[[i]]]]
        if (is.null(tr)) next
        nm <- s$sample[[i]]; col <- unname(pal[[nm]])
        traces[[length(traces) + 1]] <- list(
          x = as.list(round(tr$size, 1)), y = as.list(round(tr$y, 3)),
          name = sprintf('<span style="color:%s">%s</span>', col, html_escape(nm)), legendgroup = nm,
          mode = "lines", type = "scatter", line = list(color = col, width = 1.4),
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
      drawn <- s[s$well %in% names(r$traces), , drop = FALSE]
      lab <- drawn[is.finite(drawn$rine), , drop = FALSE]
      ymax <- max(vapply(r$traces, function(t) max(c(t$y, t$peaks$y), na.rm = TRUE), numeric(1)))
      xmax <- max(vapply(r$traces, function(t) max(t$size), numeric(1)))
      xmin <- min(vapply(r$traces, function(t) min(t$size), numeric(1)))
      n_leg <- nrow(drawn)
      # 縦の大きさ: 凡例の行 (5 サンプル/行ほど) を除いた描画域に、ラベル 1 つ 17 px が収まるように
      leg_px <- 22 * ceiling(n_leg / 5)
      height <- max(420, ceiling(76 + leg_px + nrow(lab) * 17 / 0.9))
      xrange <- NULL
      if (nrow(lab) > 0 && is.finite(ymax) && ymax > 0) {
        area <- height - 76 - leg_px
        anc <- t(vapply(lab$well, function(w) pic_rin_anchor(r$traces[[w]]), numeric(2)))
        ly <- pic_rin_spread(anc[, "y"], gap = ymax * 17 / area, lo = ymax * 0.03, hi = ymax * 0.97)
        xcol <- xmax * 1.08
        xrange <- list(log10(xmin), log10(xcol * 1.5))
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
      id <- sprintf("%srin_%d", id_prefix, k)
      ticks <- as.list(r$ladder_sizes)
      xaxis <- list(title = "size [nt]", type = "log", tickmode = "array", tickvals = ticks,
                    ticktext = as.list(format(r$ladder_sizes, trim = TRUE)), zeroline = FALSE)
      if (!is.null(xrange)) xaxis$range <- xrange
      layout <- list(
        xaxis = xaxis,
        yaxis = list(title = "sample intensity [FU]", autorange = TRUE),
        hovermode = "closest", margin = list(l = 60, r = 20, t = 30, b = 46),
        legend = list(orientation = "h", yanchor = "bottom", y = 1.02, x = 0))
      register_plot(reg, id, list(data = traces, layout = layout))
      plot <- sprintf('<div id="%s" class="pic-plot" style="height:%dpx"></div>', id, height)
    }
    if (!nzchar(plot)) next
    head <- if (length(dss) > 1) sprintf('<h4>%s</h4>', html_escape(r$stem)) else ""
    blocks <- c(blocks, sprintf('<div class="pic-plot-cell">%s%s</div>', head, plot))
  }
  if (length(blocks) == 0) return("")
  paste0('<p class="pic-note">RNA integrity of the input RNA measured on the TapeStation. The electropherogram shows each sample&rsquo;s fluorescence along fragment size; the lower marker and the ',
         '18S / 28S ribosomal RNA peaks are marked. The number next to each trace, in the same colour, is its <b>RINe</b> ',
         '(1 = degraded, 10 = intact). Click a sample in the legend to show or hide it. ',
         'Sample names follow the TapeStation sample table.</p><div class="pic-plot-grid">',
         paste(blocks, collapse = ""), '</div>')
}
