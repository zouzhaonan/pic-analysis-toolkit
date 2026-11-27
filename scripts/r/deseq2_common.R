# 部分 (--part) ごとに読み込むパッケージ。読み込みは起動の時間の大半を占める
# (ws で測ると DESeq2 が約 7 秒、clusterProfiler がさらに約 2.5 秒。cluster は DESeq2 も要らない)。
# どのコードも pkg::fn で呼んでいるので、読み込む集合を絞っても結果は変わらない
# (2026-11-07、pic193_2 / dev035 の出力が改修前と 1 バイトも同じことを確かめた)。
# 図 (ggplot2 など) と clusterProfiler は deseq2 では使わない (図はレポートが CSV から描く)。
deseq2_packages_for_part <- function(part = "all") {
  base <- c("readr", "dplyr", "tidyr", "purrr", "tibble")
  if (identical(part, "cluster")) base else c(base, "DESeq2")
}

load_required_packages <- function(part = "all") {
  required_packages <- deseq2_packages_for_part(part)
  missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      sprintf("Missing required R packages: %s", paste(missing_packages, collapse = ", ")),
      call. = FALSE
    )
  }

  invisible(lapply(required_packages, function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))
}

parse_deseq2_args <- function(args) {
  cfg <- pic_runtime_defaults()
  defaults <- list(
    deftable = NULL,
    count_dir = NULL,
    genome = NULL,
    out_dir = NULL,
    fdr = cfg$fdr,
    threads = 1L,
    part = "all",
    mem_gb = NA_character_,
    covariate = character(0)
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[[i]]

    if (!startsWith(arg, "--")) {
      stop(sprintf("Unexpected argument: %s", arg), call. = FALSE)
    }

    if (identical(arg, "--help")) {
      help_file <- file.path(analysis_root, "help", "help_deseq2.txt")
      cat(readChar(help_file, file.info(help_file)$size, useBytes = TRUE))
      cat("\n")

      registered <- "  (none)"
      if (exists("biomart_lookup_file", inherits = TRUE) &&
          file.exists(get("biomart_lookup_file", inherits = TRUE))) {
        lookup <- suppressMessages(readr::read_tsv(
          get("biomart_lookup_file", inherits = TRUE),
          show_col_types = FALSE,
          progress = FALSE,
          na = c("", "NA")
        ))
        for (nm in colnames(lookup)) {
          lookup[[nm]] <- gsub("\r", "", as.character(lookup[[nm]]), fixed = TRUE)
        }
        if (nrow(lookup) > 0 && "genome" %in% colnames(lookup) && "biomart_dataset" %in% colnames(lookup)) {
          registered_tbl <- lookup |>
            dplyr::filter(!is.na(.data$genome), .data$genome != "", !is.na(.data$biomart_dataset), .data$biomart_dataset != "") |>
            dplyr::distinct(.data$genome, .data$biomart_dataset) |>
            dplyr::arrange(.data$genome)
          if (nrow(registered_tbl) > 0) {
            registered <- paste(sprintf("  - %s (%s)", registered_tbl$genome, registered_tbl$biomart_dataset), collapse = "\n")
          }
        }
      }

      cat(
        paste(
          "Registered biomart genomes:",
          registered,
          sep = "\n"
        ),
        "\n"
      )
      quit(save = "no", status = 0)
    }

    if (grepl("=", arg, fixed = TRUE)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=(.*)$", "\\1", arg)
    } else {
      key <- sub("^--", "", arg)
      i <- i + 1
      if (i > length(args)) {
        stop(sprintf("Option requires a value: --%s", key), call. = FALSE)
      }
      value <- args[[i]]
    }

    key <- gsub("-", "_", key)
    if (!key %in% names(defaults)) {
      stop(sprintf("Unknown option: --%s", gsub("_", "-", key)), call. = FALSE)
    }

    # --covariate は繰り返してよい (カンマ区切りでもよい)。ほかは後の指定で上書き。
    if (identical(key, "covariate")) {
      defaults$covariate <- c(defaults$covariate, value)
    } else {
      defaults[[key]] <- value
    }
    i <- i + 1
  }

  if (is.null(defaults$deftable) || is.null(defaults$count_dir) || is.null(defaults$genome) || is.null(defaults$out_dir)) {
    stop("Missing required options: --deftable, --count-dir, --genome, --out-dir", call. = FALSE)
  }

  defaults$deftable <- normalizePath(defaults$deftable, winslash = "/", mustWork = TRUE)
  defaults$count_dir <- normalizePath(defaults$count_dir, winslash = "/", mustWork = TRUE)
  defaults$out_dir <- normalizePath(defaults$out_dir, winslash = "/", mustWork = FALSE)
  defaults$genome <- as.character(defaults$genome)
  defaults$fdr <- as.numeric(defaults$fdr)
  defaults$threads <- suppressWarnings(as.integer(defaults$threads))
  if (is.na(defaults$threads) || defaults$threads < 1) {
    stop("--threads must be an integer >= 1.", call. = FALSE)
  }
  defaults$ntop <- as.integer(cfg$ntop)

  # design に足す共変量 (deftable の列名)。既定は無し (design は ~group)。
  covs <- trimws(unlist(strsplit(as.character(defaults$covariate), ",", fixed = TRUE)))
  defaults$covariates <- unique(covs[nzchar(covs)])
  defaults$covariate <- NULL

  defaults$part <- as.character(defaults$part)
  if (!defaults$part %in% c("all", "fit", "rlog", "cluster")) {
    stop("--part must be one of: all, fit, rlog, cluster.", call. = FALSE)
  }
  # メモリの上限 (GB)。--mem-gb、無ければ環境変数 PIC_MEM_GB、無ければ制限なし (Inf)。
  mem <- defaults$mem_gb
  if (is.na(mem) || !nzchar(mem)) mem <- Sys.getenv("PIC_MEM_GB", unset = "")
  defaults$mem_gb <- if (nzchar(mem)) suppressWarnings(as.numeric(mem)) else Inf
  if (is.na(defaults$mem_gb) || defaults$mem_gb <= 0) {
    stop("--mem-gb (or PIC_MEM_GB) must be a number > 0.", call. = FALSE)
  }

  if (is.na(defaults$fdr) || defaults$fdr <= 0 || defaults$fdr >= 1) {
    stop("--fdr must be a number between 0 and 1.", call. = FALSE)
  }

  if (is.na(defaults$ntop) || defaults$ntop < 2) {
    stop("Internal default ntop must be an integer >= 2.", call. = FALSE)
  }

  defaults$project_name <- basename(defaults$deftable)
  defaults$project_name <- sub("\\.tsv$", "", defaults$project_name)
  defaults$project_name <- sub("^deftable_", "", defaults$project_name)

  defaults
}

read_deftable <- function(path) {
  if (!file.exists(path)) {
    stop(sprintf("deftable was not found: %s", path), call. = FALSE)
  }

  if (dir.exists(path)) {
    stop(sprintf("deftable must be a TSV file, but a directory was given: %s", path), call. = FALSE)
  }

  def <- tryCatch(
    suppressMessages(
      readr::read_tsv(
        path,
        show_col_types = FALSE,
        progress = FALSE,
        locale = readr::locale(encoding = "UTF-8")
      )
    ),
    error = function(err) {
      stop(
        paste(
          sprintf("Failed to read deftable: %s", path),
          "Make sure --deftable points to a TSV file created by pic.",
          sprintf("Original error: %s", conditionMessage(err))
        ),
        call. = FALSE
      )
    }
  )

  required_columns <- c("count_prefix", "barcode", "sample", "group")
  missing_columns <- setdiff(required_columns, colnames(def))
  if (length(missing_columns) > 0) {
    stop(
      paste(
        sprintf("deftable is missing columns: %s", paste(missing_columns, collapse = ", ")),
        "Make sure --deftable points to a genome-specific deftable TSV created by pic."
      ),
      call. = FALSE
    )
  }

  dplyr::filter(def, grepl("[^.]", .data$group))
}

build_contrasts_from_deftable <- function(def) {
  groups <- unique(def$group)
  if (length(groups) < 2) {
    stop("At least two groups are required to build contrasts.", call. = FALSE)
  }

  contrast_pairs <- combn(groups, 2, simplify = FALSE)
  # 名前は build_stats_tables() などが contrasts[[名前]] で引くので、重複してはいけない。
  # 英数字以外を「_」にすると「label-_TUC+」と「label+_TUC-」が同じ「label_TUC_」になり、
  # 全部の比較が最初の比較の名前で書かれていた。重複したら make.unique で番号を付ける。
  names(contrast_pairs) <- make.unique(vapply(
    contrast_pairs,
    function(pair) paste0(gsub("[^[:alnum:]]+", "_", pair[[1]]), "_vs_", gsub("[^[:alnum:]]+", "_", pair[[2]])),
    character(1)
  ), sep = "_")

  lapply(contrast_pairs, function(pair) c("group", pair[[1]], pair[[2]]))
}

build_num_umi_gene <- function(mat) {
  sample_columns <- colnames(mat)[5:ncol(mat)]
  summary <- tibble::tibble(
    sample = sample_columns,
    `#UMIs` = colSums(mat[, 5:ncol(mat), drop = FALSE]),
    `#Genes` = colSums(mat[, 5:ncol(mat), drop = FALSE] > 0)
  )

  dplyr::mutate(summary, `UMIs/Genes` = `#UMIs` / `#Genes`)
}

# ---------------------------------------------------------------------------
# 書きかけのファイルを残さない書き出し。
# 同じディレクトリの隠しの一時ファイル (.<名前>.pictmp-<part>-<pid>) に書いてから rename する
# (rename は同じファイルシステムの中では一度に置き換わる)。一時ファイルは「.」で始まるので、
# レポートのダウンロード一覧 (list.files は隠しファイルを拾わない) にも入らない。
# SIGKILL などで一時ファイルが残っても、同じ part の次の実行の始めに消す (pic_clean_stale_tmp)。
# ---------------------------------------------------------------------------
.pic_deseq2_state <- new.env(parent = emptyenv())
.pic_deseq2_state$tag <- "all"

pic_set_tmp_tag <- function(tag) {
  .pic_deseq2_state$tag <- tag
  invisible(tag)
}

pic_tmp_path <- function(path) {
  file.path(dirname(path), sprintf(".%s.pictmp-%s-%d", basename(path), .pic_deseq2_state$tag, Sys.getpid()))
}

pic_commit_tmp <- function(tmp, path) {
  if (!file.rename(tmp, path)) {
    unlink(tmp)
    stop(sprintf("Failed to write: %s", path), call. = FALSE)
  }
  invisible(path)
}

pic_write_csv_atomic <- function(x, path) {
  tmp <- pic_tmp_path(path)
  readr::write_csv(x, tmp)
  pic_commit_tmp(tmp, path)
}

# 前の実行が途中で殺されて残った一時ファイルを消す。part を指定すると、その part のものだけ
# (fit と rlog は同時に走るので、相手の書きかけを消さない)。"all" はどの part のものも消す。
pic_clean_stale_tmp <- function(dirs, part) {
  pat <- if (identical(part, "all")) "\\.pictmp-[a-z]+-[0-9]+$" else sprintf("\\.pictmp-%s-[0-9]+$", part)
  for (d in dirs[dir.exists(dirs)]) {
    f <- list.files(d, pattern = pat, recursive = TRUE, all.files = TRUE, full.names = TRUE)
    if (length(f)) unlink(f)
  }
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# worker (fork する子) の数。
# 減らすのは次の 2 つだけ (決め打ちの上限は置かない):
#   1. 独立した仕事 (比較など) の数より多く作らない
#   2. 親 1 本 + worker 1 本あたりのメモリ × 本数 が mem_gb (--mem-gb / PIC_MEM_GB) を超えない
# ---------------------------------------------------------------------------
pic_n_workers <- function(threads, n_jobs, per_worker_gb, mem_gb = Inf, main_gb = 0) {
  n <- min(as.integer(threads), as.integer(n_jobs))
  if (is.finite(mem_gb)) {
    n <- min(n, floor((mem_gb - main_gb) / per_worker_gb))
  }
  max(1L, as.integer(n))
}

# メモリの目安 (GB)。カウント行列の大きさ (cells = 遺伝子数 × サンプル数 / 1e6) で決める。
# 根拠は ws で測った値 (2026-11-07。pic-downstream-20261107/deseq2/README.md):
#   - fork した子は、R の GC が親から受け継いだオブジェクトに印を付けるたびにそのページを写すので、
#     子 1 本でほぼ親のヒープ 1 つ分を持つ (smaps_rollup の Private で 0.4〜1.1 GB)。
#   - --part fit の cgroup memory.peak (--threads 1/2/4/8):
#       pic193_2 (33342 遺伝子 × 36、cells 1.20、15 比較): 1.33 / 3.03 / 5.03 / 8.34 GB
#         → results の worker 1 本あたり (peak - 1 本のとき) / 本数 = 0.85 / 0.93 / 0.88 GB
#       dev035 (32103 × 10、cells 0.32、28 比較): 1.42 / 2.78 / 4.46 / 3.66 GB → 0.68 / 0.76 / 0.28 GB
#   - rlog の子の Private: pic193_2 1.11 GB、dev035 0.73 GB。
# 余裕を持たせて上に丸めた式にしている (pic193_2 で main 1.76 / results 1.04 / rlog 1.18 GB)。
deseq2_mem_estimate_gb <- function(n_genes, n_samples) {
  cells <- as.numeric(n_genes) * as.numeric(n_samples) / 1e6
  list(
    main = 1.4 + 0.3 * cells,
    results_worker = 0.8 + 0.2 * cells,
    rlog_worker = 0.7 + 0.4 * cells
  )
}

# ---------------------------------------------------------------------------
# fork した子が、親が殺されたあとに残らないようにする見張り。
# 子が長く走るとき (rlog は十数秒) は、仕事の後の pic_exit_if_orphaned だけでは仕事が終わるまで残る。
# 小さな sh を裏で立て、親がいなくなったら子を SIGKILL する。子が先に終われば sh もすぐ終わる。
# ---------------------------------------------------------------------------
pic_start_child_guard <- function(master_pid) {
  me <- Sys.getpid()
  cmd <- sprintf(
    "while kill -0 %d 2>/dev/null && kill -0 %d 2>/dev/null; do sleep 1; done; kill -0 %d 2>/dev/null || kill -9 %d 2>/dev/null",
    as.integer(master_pid), me, as.integer(master_pid), me
  )
  system2("sh", c("-c", shQuote(cmd)), wait = FALSE, stdout = FALSE, stderr = FALSE)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# --part の済んだ印と中間ファイル。
# <out_dir>/.pic_deseq2_parts_<project>/ (隠し。レポート・compare の一覧に入らない) に置く:
#   fit.done / rlog.done / cluster.done … 済んだ印。入力の鍵・書いたファイルの一覧を持つ
#   fit_stats.rds  … cluster 用 (stats の ens_gene・padj。stats_<project>.csv と同じ行の並び)
#   rlog.rds       … cluster 用 (rlog(blind = TRUE) の行列とサンプルの表)
# ---------------------------------------------------------------------------
deseq2_parts_dir <- function(args) {
  file.path(args$out_dir, sprintf(".pic_deseq2_parts_%s", args$project_name))
}

# 入力の鍵: deftable の中身・genome・count ファイルの大きさと更新時刻。
# fit と rlog の結果が同じ入力から作られたかを cluster で確かめるのに使う。
deseq2_input_key <- function(args, def) {
  cf <- file.path(args$count_dir, paste0(unique(def$count_prefix), ".txt.gz"))
  fi <- file.info(cf)
  paste(c(unname(tools::md5sum(args$deftable)), args$genome, args$count_dir,
          paste(basename(cf), fi$size, format(as.numeric(fi$mtime), nsmall = 3), sep = ":"),
          # 共変量は指定したときだけ鍵に入れる (指定しないときの鍵は前と同じ)
          if (length(args$covariates)) paste0("covariates=", paste(args$covariates, collapse = ","))),
        collapse = "|")
}

part_marker_path <- function(parts_dir, part) file.path(parts_dir, paste0(part, ".done"))

# 部分の始め: 自分の印と、自分に依る部分 (cluster) の印を消す。
begin_part <- function(parts_dir, part) {
  dir.create(parts_dir, recursive = TRUE, showWarnings = FALSE)
  unlink(part_marker_path(parts_dir, part))
  if (part %in% c("fit", "rlog")) unlink(part_marker_path(parts_dir, "cluster"))
  invisible(NULL)
}

write_part_marker <- function(parts_dir, part, key, files, fdr = NA) {
  files <- normalizePath(files[!is.na(files)], winslash = "/", mustWork = TRUE)
  path <- part_marker_path(parts_dir, part)
  tmp <- pic_tmp_path(path)
  writeLines(c(paste0("part\t", part), paste0("key\t", key), paste0("fdr\t", format(fdr, digits = 15)),
               paste0("finished\t", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
               paste0("file\t", files)), tmp)
  pic_commit_tmp(tmp, path)
}

# 前の部分が済んでいなければ、はっきりしたメッセージで終了コード 3 で止める。
require_part <- function(parts_dir, part, key, args, fdr = NULL) {
  path <- part_marker_path(parts_dir, part)
  fail <- function(why) {
    message(sprintf("[ERROR] pic deseq2 --part %s には、先に --part %s が済んでいる必要があります: %s", args$part, part, why))
    message(sprintf("[ERROR] 同じ --deftable / --count-dir / --genome / --out-dir / --covariate%s で pic deseq2 --part %s を実行してください。",
                    if (is.null(fdr)) "" else " / --fdr", part))
    quit(save = "no", status = 3)
  }
  if (!file.exists(path)) fail(sprintf("済んだ印がありません (%s)", path))
  x <- readLines(path, warn = FALSE)
  kv <- strsplit(x, "\t", fixed = TRUE)
  get1 <- function(k) { v <- vapply(kv[vapply(kv, `[`, "", 1) == k], `[`, "", 2); if (length(v)) v[[1]] else NA_character_ }
  if (!identical(get1("key"), key)) fail("入力 (deftable / count ファイル / genome / --covariate) が、その部分を回したときと違います")
  if (!is.null(fdr) && !isTRUE(all.equal(as.numeric(get1("fdr")), as.numeric(fdr)))) {
    fail(sprintf("--fdr が違います (%s では %s、今は %s)", part, get1("fdr"), format(fdr)))
  }
  files <- vapply(kv[vapply(kv, `[`, "", 1) == "file"], `[`, "", 2)
  missing <- files[!file.exists(files)]
  if (length(missing)) fail(sprintf("出力が欠けています (%s)", paste(missing, collapse = ", ")))
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# 共変量 (--covariate)。design を ~ <共変量…> + group にする。比較は今まで通り group の全ペア。
#   - 値は deftable の列 (サンプルシートの group の後ろの見出し付きの列を mapping が写したもの)。
#   - どの共変量も因子 (カテゴリ) として扱う。数字だけの列 (例: 1〜5 の提供者番号) も因子にする。
#   - 指定しないときは何も変わらない (design ~group、出力は前とバイト単位で同じ)。
# ---------------------------------------------------------------------------
PIC_RESERVED_DEFTABLE_COLUMNS <- c("count_prefix", "barcode", "sample", "group")

# 共変量の列を確かめ、サンプルごとの値の表 (sample + 共変量の列。文字列) を返す。
# 問題があれば分かりやすいメッセージで止める。
validate_covariates <- function(def, covariates) {
  if (length(covariates) == 0) return(NULL)
  bad <- intersect(covariates, PIC_RESERVED_DEFTABLE_COLUMNS)
  if (length(bad)) {
    stop(sprintf("--covariate %s: '%s' cannot be a covariate (it is a fixed deftable column).",
                 paste(bad, collapse = ","), paste(bad, collapse = "', '")), call. = FALSE)
  }
  extra <- setdiff(colnames(def), PIC_RESERVED_DEFTABLE_COLUMNS)
  missing <- setdiff(covariates, colnames(def))
  if (length(missing)) {
    stop(paste0(
      sprintf("--covariate: column(s) not found in the deftable: %s. ", paste(missing, collapse = ", ")),
      if (length(extra)) sprintf("Extra columns in the deftable: %s. ", paste(extra, collapse = ", "))
      else "The deftable has no extra columns. ",
      "Add the column (with a header) after 'group' in the sample sheet and re-run pic mapping --collect, ",
      "or add the column to summary/deftable_<project>.tsv."), call. = FALSE)
  }
  out <- data.frame(sample = as.character(def$sample), stringsAsFactors = FALSE)
  for (cv in covariates) {
    v <- as.character(def[[cv]])
    v <- trimws(v)
    empty <- is.na(v) | !nzchar(v)
    if (any(empty)) {
      stop(sprintf("--covariate %s: empty value for sample(s): %s", cv,
                   paste(unique(out$sample[empty]), collapse = ", ")), call. = FALSE)
    }
    # 同じサンプルの複数行 (複数 barcode の合算) で値が食い違っていないか
    per <- tapply(v, out$sample, function(x) length(unique(x)))
    if (any(per > 1)) {
      stop(sprintf("--covariate %s: sample(s) have different values in different rows: %s", cv,
                   paste(names(per)[per > 1], collapse = ", ")), call. = FALSE)
    }
    out[[cv]] <- v
  }
  out[!duplicated(out$sample), , drop = FALSE]
}

# design に渡す水準名。R の名前として正しい (make.names() で変わらない) ならそのまま、
# そうでなければ <prefix>1, <prefix>2, … (+ - などを含むとき。DESeq2 は make.names で重複を調べる)。
pic_design_levels <- function(levels_orig, prefix = "g") {
  if (identical(make.names(levels_orig, unique = TRUE), levels_orig)) return(levels_orig)
  sprintf("%s%d", prefix, seq_along(levels_orig))
}

# design の式に使う共変量の列名。R の名前として正しく、ほかと重ならなければそのまま、
# そうでなければ cov1, cov2, … (式・colData 用。表示には元の名前を使う)。
pic_covariate_design_names <- function(covariates) {
  ok <- make.names(covariates) == covariates & !covariates %in% c("group", "sample", "sizeFactor")
  ifelse(ok, covariates, sprintf("cov%d", seq_along(covariates)))
}

# 表示用の design の式 (元の列名。R の名前として正しくない列名は `` で囲む)。
pic_design_label <- function(covariates) {
  if (length(covariates) == 0) return("~group")
  nm <- ifelse(make.names(covariates) == covariates, covariates, paste0("`", covariates, "`"))
  paste0("~", paste(c(nm, "group"), collapse = " + "))
}

# design が full rank か、残差の自由度があるかを確かめる。だめなら原因の共変量を挙げて止める。
check_covariate_design <- function(col_df, cov_names, cov_labels) {
  mm_of <- function(terms) stats::model.matrix(stats::as.formula(paste("~", paste(terms, collapse = " + "))), col_df)
  full_rank <- function(mm) qr(mm)$rank == ncol(mm)
  mm <- mm_of(c(cov_names, "group"))
  if (!full_rank(mm)) {
    confounded <- cov_labels[vapply(cov_names, function(cn) !full_rank(mm_of(c(cn, "group"))), logical(1))]
    stop(paste0(
      sprintf("The design %s is not full rank: ", pic_design_label(cov_labels)),
      if (length(confounded))
        sprintf("covariate(s) %s are confounded with group (each value occurs in only one group, or the values separate the groups). ",
                paste(confounded, collapse = ", "))
      else "the covariates are confounded with each other or with group. ",
      "Group differences cannot be separated from the covariate; drop --covariate or fix the values."), call. = FALSE)
  }
  if (nrow(mm) - ncol(mm) < 1) {
    stop(sprintf(paste0("The design %s leaves no residual degrees of freedom (%d samples, %d coefficients); ",
                        "dispersion cannot be estimated. Drop --covariate."),
                 pic_design_label(cov_labels), nrow(mm), ncol(mm)), call. = FALSE)
  }
  invisible(TRUE)
}

# 共変量の効果を発現行列 (log の尺度。rlog など) から除く。group の差は design で残す。
#   limma::removeBatchEffect(x, batch = cov1, batch2 = cov2, covariates = <3 つ目以降>, design = model.matrix(~group))
# 3 つ目以降の共変量は、limma が batch に使うのと同じ和が 0 の対比 (contr.sum) の行列にして covariates に渡す。
pic_remove_covariate_effects <- function(x, label, covariates) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("R package 'limma' is required for --covariate (covariate-removed PCA / correlation).", call. = FALSE)
  }
  lab <- as.data.frame(label, stringsAsFactors = FALSE)
  lab <- lab[match(colnames(x), as.character(lab$sample)), , drop = FALSE]
  grp <- factor(as.character(lab$group), levels = unique(as.character(lab$group)))
  design <- stats::model.matrix(~grp)
  fac <- lapply(covariates, function(cv) factor(as.character(lab[[cv]]), levels = unique(as.character(lab[[cv]]))))
  more <- NULL
  if (length(fac) > 2) {
    more <- do.call(cbind, lapply(fac[-(1:2)], function(f) {
      stats::contrasts(f) <- stats::contr.sum(levels(f))
      stats::model.matrix(~f)[, -1, drop = FALSE]
    }))
  }
  limma::removeBatchEffect(x, batch = fac[[1]], batch2 = if (length(fac) >= 2) fac[[2]] else NULL,
                           covariates = more, design = design)
}

# Materials & Methods に書く、補正の式。
pic_remove_covariate_label <- function(covariates) {
  nm <- ifelse(make.names(covariates) == covariates, covariates, paste0("`", covariates, "`"))
  args <- c(sprintf("batch = %s", nm[1]),
            if (length(nm) >= 2) sprintf("batch2 = %s", nm[2]),
            if (length(nm) >= 3) sprintf("covariates = <%s>", paste(nm[-(1:2)], collapse = ", ")),
            "design = model.matrix(~group)")
  sprintf("limma::removeBatchEffect(x, %s)", paste(args, collapse = ", "))
}
