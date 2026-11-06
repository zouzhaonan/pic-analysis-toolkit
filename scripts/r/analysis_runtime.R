# 役割:
#   R 解析の実行時設定（固定デフォルト値）と
#   パス解決をまとめる。

pic_runtime_defaults <- function() {
  pic_plot_spec()$defaults
}

pic_runtime_paths <- function(analysis_root) {
  lib_dir <- Sys.getenv("PIC_LIB", unset = "")
  if (identical(lib_dir, "")) {
    home_dir <- Sys.getenv("HOME", unset = "~")
    lib_dir <- file.path(home_dir, "local", "lib", "pic")
    Sys.setenv(PIC_LIB = lib_dir)
  }
  if (!dir.exists(lib_dir)) {
    ok <- dir.create(lib_dir, recursive = TRUE, showWarnings = FALSE)
    if (!isTRUE(ok) && !dir.exists(lib_dir)) {
      stop(
        sprintf("Failed to create PIC_LIB directory: %s", lib_dir),
        call. = FALSE
      )
    }
  }
  register_dir <- file.path(lib_dir, "register")

  list(
    lib_dir = lib_dir,
    help_dir = file.path(analysis_root, "help"),
    register_dir = register_dir,
    biomart_lookup_file = file.path(register_dir, "biomart_lookup.tsv"),
    genome_map_file = file.path(register_dir, "genome_enrichment_map.tsv"),
    biomart_cache_dir = file.path(lib_dir, "biomart"),
    biomart_ortholog_dir = file.path(lib_dir, "biomart", "ortholog"),
    enrichment_lib_dir = file.path(lib_dir, "enrichment")
  )
}

# ---------------------------------------------------------------------------
# fork した子 (parallel::mclapply / mcparallel) が、親が殺されたあとに残らないようにする。
# 親が SIGTERM などで先に死ぬと、子は結果を親へ送れず (sendMaster の失敗)、そのまま眠って残る
# (2026-11-07 に enrich で確かめた。子の終わりの mcexit が親の合図を待ち続ける)。子は job の前と後に
# これを呼び、親が替わって (いなくなって) いたら、すぐに終わる。master_pid は fork の前に親で Sys.getpid() を取っておいたもの。
# ---------------------------------------------------------------------------
pic_getppid <- function() {
  st <- "/proc/self/stat"
  if (file.exists(st)) {
    x <- tryCatch(readLines(st, warn = FALSE), error = function(e) "")
    # pid (comm) state ppid ... 。comm に空白や括弧が入りうるので、最後の ")" の後ろを読む
    rest <- sub("^.*\\) ", "", x[[1]])
    return(suppressWarnings(as.integer(strsplit(rest, " ", fixed = TRUE)[[1]][[2]])))
  }
  out <- tryCatch(system2("ps", c("-o", "ppid=", "-p", Sys.getpid()), stdout = TRUE, stderr = FALSE),
                  error = function(e) "")
  suppressWarnings(as.integer(trimws(out[[1]])))
}

pic_exit_if_orphaned <- function(master_pid) {
  pp <- pic_getppid()
  if (!is.na(pp) && !identical(as.integer(pp), as.integer(master_pid))) {
    # parallel:::mcexit は親からの合図 (SIGUSR1) を待ってから終わるので、親がいないと眠ったまま残る。
    # 自分を SIGKILL で終える (書きかけの CSV は一時ファイルなので残らない)。
    tools::pskill(Sys.getpid(), tools::SIGKILL)
  }
  invisible(NULL)
}
