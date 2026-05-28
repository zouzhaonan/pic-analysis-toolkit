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
    help_dir = file.path(lib_dir, "help"),
    register_dir = register_dir,
    biomart_lookup_file = file.path(register_dir, "biomart_lookup.tsv"),
    genome_map_file = file.path(register_dir, "genome_enrichment_map.tsv"),
    biomart_cache_dir = file.path(lib_dir, "biomart"),
    biomart_ortholog_dir = file.path(lib_dir, "biomart", "ortholog"),
    enrichment_lib_dir = file.path(lib_dir, "enrichment")
  )
}
