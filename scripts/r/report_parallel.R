# 役割: report の重い所 (節ごとの組み立て・GSEA/ORA の CSV 読みと図の spec・埋め込みファイルの gzip+base64)
#       を fork で並べる。出力の HTML は 1 本で回したときと 1 バイトも違わないようにする。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が report_build.R 経由で source する。
#
# 並べ方:
#   - 子はそれぞれ新しい空のレジストリ (reg) に登録し、終わったら中身 (plots / expr / pca / files) を返す。
#     親は返ってきた中身を、1 本で回したときと同じ順 (節の順) にレジストリへ入れ直す
#     (pic_report_reg_merge)。plots / expr / pca は同じ id を後から入れると値だけ替わり位置は残る、
#     files は先に入ったものが残る、という register_plot / register_file の振る舞いをそのまま写す。
#   - 図の spec は子の中で JSON にしておき (class "json" の文字列)、親の render では jsonlite がそのまま書く。
#     要素ごとに書いた JSON をつないだものは、まとめて書いたものと同じ文字列になる。
#   - worker の数は --threads に従う。減らすのは「仕事の数より多く作らない」と
#     「1 本あたりのメモリ × 本数が --mem-gb (無ければ PIC_MEM_GB) を超えない」のときだけ。

# worker 1 本あたりのメモリ (GB)。2026-11-07 に ws で、fork した子ごとの私有メモリ (smaps_rollup の
# Private_Clean + Private_Dirty の最大。親と共有のままのページは数えない) を測った値に基づく:
#   pic193_2 (mm10 36 サンプル、DEG 7,988 行、GSEA CSV 105 本): heatmap の子 0.74 GB (dist が 行数^2/2 × 8 バイト
#   = 0.26 GB と hclust・HTML の文字列)、Expression / MA-volcano の子 0.37-0.39 GB、GSEA・ORA の束 0.13-0.16 GB。
#   dev035 (hg38 10 サンプル、DEG 2,665 行): 0.42-0.43 GB。親は 0.7-0.95 GB。
# 重い子は 1 本 (heatmap) で残りは軽いので、平均より少し上の 0.5 GB を 1 本あたりとする。
# (DEG の行数が 2 万を超えると heatmap の子だけで 3 GB を超える。そのときは --mem-gb に余裕を持たせる)
PIC_REPORT_WORKER_GB <- 0.5

pic_report_threads <- function() {
  n <- suppressWarnings(as.integer(getOption("pic.report.threads", 1L)))
  if (length(n) != 1 || is.na(n) || n < 1L) 1L else n
}

pic_report_mem_gb <- function() {
  m <- getOption("pic.report.mem_gb", NULL)
  if (is.null(m) || !nzchar(as.character(m))) m <- Sys.getenv("PIC_MEM_GB", "")
  m <- suppressWarnings(as.numeric(m))
  if (length(m) != 1 || is.na(m) || m <= 0) NA_real_ else m
}

# 同時に走らせる worker の数。
pic_report_n_workers <- function(n_jobs, per_worker_gb = PIC_REPORT_WORKER_GB) {
  w <- min(pic_report_threads(), n_jobs)
  mem <- pic_report_mem_gb()
  if (!is.na(mem) && per_worker_gb > 0) w <- min(w, floor(mem / per_worker_gb))
  max(1L, as.integer(w))
}

# --- fork した子の後始末 -------------------------------------------------------
# 親が SIGTERM などで先に死ぬと、子は結果を親へ送れず、parallel の mcexit が親の合図を待って眠ったまま残る。
# 子は仕事の前と後にこれを呼び、親が替わって (いなくなって) いたら自分を SIGKILL で終える。
# (scripts/r/analysis_runtime.R の pic_exit_if_orphaned と同じもの。report はそちらを source しないので写しを置く)
pic_report_getppid <- function() {
  st <- "/proc/self/stat"
  if (file.exists(st)) {
    x <- tryCatch(readLines(st, warn = FALSE), error = function(e) "")
    rest <- sub("^.*\\) ", "", x[[1]])
    return(suppressWarnings(as.integer(strsplit(rest, " ", fixed = TRUE)[[1]][[2]])))
  }
  out <- tryCatch(system2("ps", c("-o", "ppid=", "-p", Sys.getpid()), stdout = TRUE, stderr = FALSE),
                  error = function(e) "")
  suppressWarnings(as.integer(trimws(out[[1]])))
}

pic_report_exit_if_orphaned <- function(master_pid) {
  pp <- pic_report_getppid()
  if (!is.na(pp) && !identical(as.integer(pp), as.integer(master_pid))) {
    tools::pskill(Sys.getpid(), tools::SIGKILL)
  }
  invisible(NULL)
}

# 仕事 (引数なしの関数) の list を回し、同じ順の結果の list を返す。
# order: 先に始める順 (重いものを先に)。NULL なら与えた順。
# worker が 1 なら fork せずに順に回す。子でのエラーは親で stop にする (1 本で回したときと同じく report は失敗する)。
pic_report_run_jobs <- function(jobs, order = NULL, per_worker_gb = PIC_REPORT_WORKER_GB) {
  n <- length(jobs)
  if (n == 0) return(list())
  w <- pic_report_n_workers(n, per_worker_gb)
  if (w <= 1L) return(lapply(jobs, function(f) f()))
  if (is.null(order)) order <- seq_len(n)
  master_pid <- Sys.getpid()
  res <- parallel::mclapply(order, function(i) {
    pic_report_exit_if_orphaned(master_pid)
    v <- jobs[[i]]()
    pic_report_exit_if_orphaned(master_pid)
    list(ok = v)   # 仕事が NULL を返しても、子が死んだとき (NULL) と見分けられるように包む
  }, mc.cores = w, mc.preschedule = FALSE, mc.set.seed = FALSE)
  out <- vector("list", n)
  for (k in seq_along(order)) {
    r <- res[[k]]
    if (inherits(r, "try-error")) {
      stop(sprintf("report: 並べた仕事 %d がエラーで終わりました: %s", order[[k]],
                   conditionMessage(attr(r, "condition"))), call. = FALSE)
    }
    if (is.null(r)) {
      stop(sprintf("report: 並べた仕事 %d の子が結果を返さずに終わりました (メモリ不足で殺された可能性。--threads か --mem-gb を下げてください)",
                   order[[k]]), call. = FALSE)
    }
    out[order[[k]]] <- list(r$ok)
  }
  out
}

# costs に従って n 個を k 個の束に分ける (重い順に、いま一番軽い束へ入れる)。束の中は元の順。
pic_report_chunks <- function(costs, k) {
  n <- length(costs)
  k <- max(1L, min(as.integer(k), n))
  load <- numeric(k); bins <- vector("list", k)
  for (i in order(costs, decreasing = TRUE)) {
    b <- which.min(load); load[[b]] <- load[[b]] + costs[[i]]; bins[[b]] <- c(bins[[b]], i)
  }
  lapply(bins[lengths(bins) > 0], sort)
}

# X の各要素に FUN をかけた list (X と同じ順)。束にして並べる。costs は束の釣り合いにだけ使う。
pic_report_pmap <- function(X, FUN, costs = NULL, per_worker_gb = PIC_REPORT_WORKER_GB) {
  n <- length(X)
  if (n == 0) return(list())
  w <- pic_report_n_workers(n, per_worker_gb)
  if (w <= 1L) return(lapply(X, FUN))
  if (is.null(costs)) costs <- rep(1, n)
  bins <- pic_report_chunks(costs, w * 2L)
  jobs <- lapply(bins, function(idx) function() lapply(X[idx], FUN))
  bin_cost <- vapply(bins, function(idx) sum(costs[idx]), numeric(1))
  res <- pic_report_run_jobs(jobs, order = order(bin_cost, decreasing = TRUE), per_worker_gb = per_worker_gb)
  out <- vector("list", n)
  for (b in seq_along(bins)) out[bins[[b]]] <- res[[b]]
  out
}

# --- レジストリの写しと入れ直し -------------------------------------------------

# レジストリに入れる値を、render が書くのと同じ JSON にしておく (render の pic_to_json が verbatim で書く)。
pic_json_verbatim <- function(x, digits) {
  if (inherits(x, "json")) return(x)
  structure(as.character(pic_to_json(x, digits = digits)), class = "json")
}

# render_report_page の digits と揃える (plots / expr = 6, pca = 4)。
pic_report_reg_snapshot <- function(r) {
  list(plots = lapply(r$plots, pic_json_verbatim, digits = 6),
       expr  = lapply(r$expr,  pic_json_verbatim, digits = 6),
       pca   = lapply(r$pca,   pic_json_verbatim, digits = 4),
       files = r$files)
}

pic_report_reg_merge <- function(reg, snap) {
  for (k in names(snap$plots)) reg$plots[[k]] <- snap$plots[[k]]
  for (k in names(snap$expr))  reg$expr[[k]]  <- snap$expr[[k]]
  for (k in names(snap$pca))   reg$pca[[k]]   <- snap$pca[[k]]
  if (!is.null(snap$files)) {
    if (is.null(reg$files)) reg$files <- list()
    for (k in names(snap$files)) if (is.null(reg$files[[k]])) reg$files[[k]] <- snap$files[[k]]
  }
  invisible(reg)
}

# 節を組み立てる仕事 (function(reg) が html などを返す) を、新しいレジストリで回すように包む。
# 子でも親でも同じ包みで回す (1 本のときも同じ道を通るので、--threads で出力が変わらない)。
pic_report_section_job <- function(f) {
  force(f)
  function() {
    r <- pic_report_registry()
    old <- options(pic.report.reg = r)
    on.exit(options(old), add = TRUE)
    v <- f(r)
    list(value = v, reg = pic_report_reg_snapshot(r))
  }
}

# --- 埋め込みファイルの gzip+base64 を先に並べて作る ------------------------------

# register_file と同じやり方で gz を作る。失敗は NULL。
pic_report_file_gz <- function(abs) {
  tryCatch({
    raw <- readBin(abs, "raw", n = file.info(abs)$size)
    jsonlite::base64_enc(memCompress(raw, "gzip"))
  }, error = function(e) NULL)
}

# rel (projdir からの相対パス) の gz を並べて作り、register_file が使う置き場に入れる。
# 登録済みのもの・除外されるもの・無いものは作らない。
pic_report_prefetch_gz <- function(rels) {
  reg <- getOption("pic.report.reg"); proj <- getOption("pic.report.projdir")
  if (is.null(reg) || is.null(proj) || length(rels) == 0) return(invisible(NULL))
  cache <- pic_report_gz_cache()
  todo <- character(0)
  for (rel in unique(rels)) {
    if (is.null(rel) || !nzchar(rel) || file_excluded(rel)) next
    id <- gsub("[^A-Za-z0-9]+", "_", rel)
    if (!is.null(reg$files) && !is.null(reg$files[[id]])) next
    abs <- file.path(proj, rel)
    if (exists(abs, envir = cache, inherits = FALSE)) next
    if (!file.exists(abs) || dir.exists(abs)) next
    todo <- c(todo, abs)
  }
  if (length(todo) == 0) return(invisible(NULL))
  sizes <- as.numeric(file.info(todo)$size); sizes[is.na(sizes)] <- 0
  gz <- pic_report_pmap(todo, function(a) list(pic_report_file_gz(a)), costs = sizes + 1e4)
  for (i in seq_along(todo)) assign(todo[[i]], gz[[i]][[1]], envir = cache)
  invisible(NULL)
}

# 先に作った gz の置き場 (絶対パス -> gz)。値が NULL のものは作るのに失敗したもの。
pic_report_gz_cache <- function() {
  e <- getOption("pic.report.gz_cache")
  if (is.null(e)) { e <- new.env(parent = emptyenv()); options(pic.report.gz_cache = e) }
  e
}
