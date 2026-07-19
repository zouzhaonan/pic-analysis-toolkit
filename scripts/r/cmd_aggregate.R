#!/usr/bin/env Rscript

# 役割:
#   遺伝子全体 (TSS→TES) のメタジーン集計プロファイルを計算する。
#   pic の CPM 正規化 bigwig (bw/<count_prefix>_umi.cpm.bw) を、各遺伝子を N
#   分割したビンに対して UCSC bigWigAverageOverBed で平均取得し、全遺伝子で
#   平均してサンプルごとのプロファイルにする (deepTools 非依存で高速・省メモリ)。
# 入力:
#   --out-dir <dir>   pic mapping/all の出力 (bw/ と summary/deftable_<run>_<genome>.tsv)
#   --run-name <name> 省略時は summary/mapping_sum__<run>.tsv から推定
#   --genome <g>      特定 genome のみ (省略時は全 genome)
#   --threads <int>   サンプルごとの bigWigAverageOverBed を mclapply で並列実行
#   --bins <int>      遺伝子本体の分割数 (default: 100)
# 出力:
#   <out>/summary/aggregate_profile_<run>_<genome>.csv (PNG は生成しない)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

parse_args <- function(a) {
  out <- list(out_dir = NULL, run_name = NULL, genome = NULL, group = NULL,
              threads = 8L, bins = 100L, flank = 1000L, flank_bins = 10L)
  i <- 1
  while (i <= length(a)) {
    k <- a[[i]]
    if (k == "--out-dir") { out$out_dir <- a[[i + 1]]; i <- i + 2 }
    else if (k == "--run-name") { out$run_name <- a[[i + 1]]; i <- i + 2 }
    else if (k == "--genome") { out$genome <- a[[i + 1]]; i <- i + 2 }
    else if (k == "--group") { out$group <- a[[i + 1]]; i <- i + 2 }
    else if (k == "--threads") { out$threads <- as.integer(a[[i + 1]]); i <- i + 2 }
    else if (k == "--bins") { out$bins <- as.integer(a[[i + 1]]); i <- i + 2 }
    else if (k == "--flank") { out$flank <- as.integer(a[[i + 1]]); i <- i + 2 }
    else if (k == "--flank-bins") { out$flank_bins <- as.integer(a[[i + 1]]); i <- i + 2 }
    else stop(sprintf("Unknown option: %s", k), call. = FALSE)
  }
  if (is.null(out$out_dir)) stop("aggregate requires --out-dir", call. = FALSE)
  out
}

# chrom.sizes を named integer で読み込む (フランク領域の端をクリップするため)。
read_chrom_sizes <- function(genome) {
  f <- file.path(pic_lib_dir(), "chrom_size", paste0(genome, ".chrom.sizes"))
  if (!file.exists(f)) return(NULL)
  cs <- fread(f, header = FALSE, showProgress = FALSE, select = c(1, 2), col.names = c("chr", "size"))
  stats::setNames(as.numeric(cs$size), cs$chr)
}

pic_lib_dir <- function() {
  l <- Sys.getenv("PIC_LIB")
  if (nzchar(l)) return(l)
  file.path(Sys.getenv("HOME"), "local/lib/pic")
}

resolve_gtf <- function(genome) {
  g <- file.path(pic_lib_dir(), "gtf", paste0(genome, ".gtf"))
  if (file.exists(g)) return(g)
  if (file.exists(paste0(g, ".gz"))) return(paste0(g, ".gz"))
  ""
}

# GTF から遺伝子本体 (chr,start,end,strand,gene) を取り出す。
read_gene_bodies <- function(gtf) {
  reader <- if (grepl("\\.gz$", gtf)) "gzcat" else "cat"
  dt <- fread(cmd = sprintf("%s '%s' | grep -v '^#'", reader, gtf),
              sep = "\t", header = FALSE, quote = "", showProgress = FALSE,
              select = c(1, 3, 4, 5, 7, 9),
              col.names = c("chr", "feature", "start", "end", "strand", "attr"))
  dt[, gene := sub('.*gene_id "?([^";]+)"?.*', "\\1", attr)]
  g <- dt[feature == "gene"]
  if (nrow(g) == 0) {
    g <- dt[, .(chr = chr[1], start = min(start), end = max(end), strand = strand[1]), by = gene]
  } else {
    g <- g[, .(chr = chr[1], start = min(start), end = max(end), strand = strand[1]), by = gene]
  }
  g[!is.na(chr) & end > start & strand %in% c("+", "-")]
}

# globalpos (1..total) -> region / plot 位置 (pos) の対応表。
# TSS=0, TES=100。上流フランクは負、下流フランクは 100 超 (各フランクは
# body ビンと同じ視覚幅 = 100*flank_bins/nbins)。
position_map <- function(nbins, fbins) {
  total <- 2L * fbins + nbins
  span <- 100 * fbins / nbins
  gp <- seq_len(total)
  region <- ifelse(gp <= fbins, "up", ifelse(gp <= fbins + nbins, "body", "down"))
  pos <- numeric(total)
  up <- gp[region == "up"]; pos[up] <- -span + (up - 0.5) / fbins * span
  bd <- gp[region == "body"]; bidx <- bd - fbins; pos[bd] <- (bidx - 0.5) / nbins * 100
  dn <- gp[region == "down"]; didx <- dn - (fbins + nbins); pos[dn] <- 100 + (didx - 0.5) / fbins * span
  data.table(globalpos = gp, region = region, pos = pos)
}

# 遺伝子本体 (N ビン, TSS→TES で scale) + 上流/下流フランク (flank bp, fbins) を
# タイル化した BED を作る。name = "<gene index>|<globalpos>"。
# globalpos は TSS→TES→下流の順 (strand 補正済み, 1..total)。
build_flanked_bed <- function(g, nbins, flank, fbins, chrom_sizes, bedfile) {
  g <- g[(end - start + 1) >= nbins]
  if (nrow(g) == 0) return(0L)
  g[, gi := .I]
  N <- nrow(g)
  total <- 2L * fbins + nbins
  fw <- flank / fbins
  start0 <- g$start - 1L; end0 <- g$end

  mk <- function(cnt, gs, ge, ginum, idx) {
    data.table(gi = g$gi[idx], chr = g$chr[idx], strand = g$strand[idx],
               bstart = as.integer(floor(gs)), bend = as.integer(floor(ge)), ginum = ginum)
  }
  # 左フランク (genomic, s0 の左)
  idxL <- rep.int(seq_len(N), rep.int(fbins, N)); iL <- rep.int(0:(fbins - 1), N)
  left <- mk(NULL, start0[idxL] - flank + iL * fw, start0[idxL] - flank + (iL + 1) * fw, iL + 1L, idxL)
  # body (genomic, s0..e0, nbins に scale)
  idxB <- rep.int(seq_len(N), rep.int(nbins, N)); jB <- rep.int(0:(nbins - 1), N)
  wB <- (end0[idxB] - start0[idxB]) / nbins
  body <- mk(NULL, start0[idxB] + jB * wB, start0[idxB] + (jB + 1) * wB, fbins + jB + 1L, idxB)
  # 右フランク (genomic, e0 の右)
  idxR <- rep.int(seq_len(N), rep.int(fbins, N)); kR <- rep.int(0:(fbins - 1), N)
  right <- mk(NULL, end0[idxR] + kR * fw, end0[idxR] + (kR + 1) * fw, fbins + nbins + kR + 1L, idxR)

  dt <- rbindlist(list(left, body, right))
  # strand 補正: + はそのまま、- は genomic 昇順を反転して TSS→TES に。
  dt[, globalpos := ifelse(strand == "+", ginum, total - ginum + 1L)]
  # chrom 端でクリップ
  if (!is.null(chrom_sizes)) {
    cl <- chrom_sizes[dt$chr]; cl[is.na(cl)] <- .Machine$integer.max
    dt[, bend := pmin(bend, as.integer(cl))]
  }
  dt[, bstart := pmax(0L, bstart)]
  dt <- dt[bend > bstart]
  fwrite(dt[, .(chr, bstart, bend, name = paste0(gi, "|", globalpos))],
         bedfile, sep = "\t", col.names = FALSE, quote = FALSE)
  N
}

# 1 サンプルの bigwig → ビンごとの全遺伝子平均プロファイル (length nbins)。
sample_profile <- function(bwaob, bw, bed, nbins, tmp_dir) {
  outtab <- tempfile(tmpdir = tmp_dir, fileext = ".tab")
  st <- suppressWarnings(system2(bwaob, c(shQuote(bw), shQuote(bed), shQuote(outtab)),
                                 stdout = FALSE, stderr = FALSE))
  if (st != 0 || !file.exists(outtab)) { unlink(outtab); return(NULL) }
  tab <- fread(outtab, header = FALSE, showProgress = FALSE,
               col.names = c("name", "size", "covered", "sum", "mean0", "mean"))
  unlink(outtab)
  if (nrow(tab) == 0) return(NULL)
  parts <- tstrsplit(tab$name, "|", fixed = TRUE)
  tab[, gi := parts[[1]]]
  tab[, binpos := as.integer(parts[[2]])]
  tab[, gtot := sum(mean0), by = gi]          # 無信号遺伝子は除外 (skipZeros 相当)
  prof <- tab[gtot > 0, .(value = mean(mean0)), by = binpos][order(binpos)]
  setkey(prof, binpos)
  prof[.(seq_len(nbins))]$value
}

aggregate_one_genome <- function(out_dir, genome, project, deftable, bwaob, nbins, flank, fbins, tmp_dir, threads = 1L, group_pat = NULL) {
  def <- fread(deftable, sep = "\t", header = TRUE, showProgress = FALSE)
  # 期待列: count_prefix, barcode, sample, group
  setnames(def, tolower(names(def)))
  if (!all(c("count_prefix", "sample", "group") %in% names(def))) {
    message(sprintf("[WARN] aggregate: deftable の列が不足 (%s)", deftable)); return(invisible(NULL))
  }
  if (!is.null(group_pat) && nzchar(group_pat)) {
    keep <- grepl(group_pat, def$group) | grepl(group_pat, def$sample)
    def <- def[keep]
    if (nrow(def) == 0) { message(sprintf("[WARN] aggregate: --group '%s' に一致なし", group_pat)); return(invisible(NULL)) }
  }
  gtf <- resolve_gtf(genome)
  if (!nzchar(gtf)) { message(sprintf("[WARN] aggregate: GTF なし genome=%s", genome)); return(invisible(NULL)) }

  message(sprintf("[INFO] aggregate: gene bodies を読み込みます (genome=%s)", genome))
  genes <- read_gene_bodies(gtf)
  if (nrow(genes) == 0) { message("[WARN] aggregate: 有効な遺伝子なし"); return(invisible(NULL)) }
  chrom_sizes <- read_chrom_sizes(genome)
  total <- 2L * fbins + nbins
  pmap <- position_map(nbins, fbins)

  bed <- tempfile(tmpdir = tmp_dir, fileext = ".bed")
  n_used <- build_flanked_bed(genes, nbins, flank, fbins, chrom_sizes, bed)
  message(sprintf("[INFO] aggregate: %d 遺伝子 × %d bins (body %d + flank %d×2, ±%dbp)", n_used, total, nbins, fbins, flank))

  bw_dir <- file.path(out_dir, "bw")
  # サンプルは互いに独立 (各自の bigwig を read-only の bed に対して集計) なので、
  # bigWigAverageOverBed 呼び出しを mclapply で並列化する (bwaob 自体は 1 スレッド)。
  ncore <- max(1L, min(threads, nrow(def)))
  results <- parallel::mclapply(seq_len(nrow(def)), function(i) {
    cp <- as.character(def$count_prefix[i]); sample <- as.character(def$sample[i]); group <- as.character(def$group[i])
    bw <- file.path(bw_dir, paste0(cp, "_umi.cpm.bw"))
    if (!file.exists(bw)) bw <- file.path(bw_dir, paste0(sample, "_umi.cpm.bw"))
    if (!file.exists(bw)) { message(sprintf("[WARN] aggregate: bigwig なし (skip): %s", sample)); return(NULL) }
    v <- sample_profile(bwaob, bw, bed, total, tmp_dir)
    if (is.null(v)) { message(sprintf("[WARN] aggregate: プロファイル取得失敗 (skip): %s", sample)); return(NULL) }
    v[!is.finite(v)] <- 0
    message(sprintf("[INFO] aggregate: %s 完了", sample))
    data.table(sample = sample, group = group, globalpos = seq_len(total), value = v)
  }, mc.cores = ncore)
  rows <- results[vapply(results, data.table::is.data.table, logical(1))]
  unlink(bed)
  if (length(rows) == 0) { message("[WARN] aggregate: 出力なし"); return(invisible(NULL)) }

  d <- rbindlist(rows)
  d <- merge(d, pmap, by = "globalpos", all.x = TRUE)
  d[, flank_bp := flank]
  setnames(d, "globalpos", "binpos")
  setorder(d, sample, binpos)
  # PNG は不要 (図はレポートが CSV から plotly で動的生成)。CSV は集計テーブルとして
  # summary/ に出力 (mapping_sum/deftable と同じ集約先)。
  summary_dir <- file.path(out_dir, "summary")
  dir.create(summary_dir, recursive = TRUE, showWarnings = FALSE)
  out_csv <- file.path(summary_dir, sprintf("aggregate_profile_%s.csv", project))
  fwrite(d[, .(sample, group, binpos, region, pos, value, flank_bp)], out_csv)
  message(sprintf("[INFO] aggregate: 出力 -> %s", out_csv))
}

# summary/ を優先し、無ければ out_dir 直下も探す (旧レイアウト互換)。
list_from_summary <- function(out_dir, pattern) {
  hits <- list.files(file.path(out_dir, "summary"), pattern = pattern, full.names = TRUE)
  if (length(hits) == 0) hits <- list.files(out_dir, pattern = pattern, full.names = TRUE)
  hits
}

main <- function() {
  args <- parse_args(commandArgs(trailingOnly = TRUE))
  out_dir <- normalizePath(args$out_dir, mustWork = TRUE)

  bwaob <- Sys.which("bigWigAverageOverBed")
  if (!nzchar(bwaob)) stop("bigWigAverageOverBed が見つかりません (pic の conda 環境を確認)", call. = FALSE)

  run_name <- args$run_name
  if (is.null(run_name)) {
    ms <- list_from_summary(out_dir, "^mapping_sum__.*\\.tsv$")
    if (length(ms) > 0) run_name <- sub("\\.tsv$", "", sub("^mapping_sum__", "", basename(ms[[1]])))
  }
  if (is.null(run_name) || !nzchar(run_name)) stop("run-name を特定できません。--run-name を指定してください。", call. = FALSE)

  deftables <- list_from_summary(out_dir, sprintf("^deftable_%s_.*\\.tsv$", run_name))
  if (length(deftables) == 0) stop(sprintf("deftable が見つかりません: %s/summary/deftable_%s_*.tsv", out_dir, run_name), call. = FALSE)

  tmp_dir <- file.path(out_dir, "tmp", "aggregate")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  for (df in deftables) {
    genome <- sub(sprintf("^deftable_%s_", run_name), "", sub("\\.tsv$", "", basename(df)))
    if (!is.null(args$genome) && genome != args$genome) next
    message(sprintf("[INFO] aggregate: genome=%s", genome))
    project <- sprintf("%s_%s", run_name, genome)
    aggregate_one_genome(out_dir, genome, project, df, bwaob, args$bins, args$flank, args$flank_bins, tmp_dir, args$threads, args$group)
  }
  unlink(tmp_dir, recursive = TRUE)
  message(sprintf("[INFO] aggregate: 完了 (out-dir=%s)", out_dir))
}

main()
