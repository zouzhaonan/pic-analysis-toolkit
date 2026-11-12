#!/usr/bin/env Rscript

# 役割:
#   ChIP-Atlas の Enrichment Analysis (EA) と PAGE を DDBJ の WABI に投げ、WABI が返す tsv / html / log を
#   そのまま保存する (まとめの表・図・向きの付け直しは作らない)。
#   - EA  : 比較ごとに「分子で高い DEG (A)」と「分母で高い DEG (B)」を互いの背景にして投げる
#           (WABI の typeA=gene / typeB=userlist)。どちらかが 0 個なら計算できないので飛ばし、理由を残す。
#   - PAGE: 比較の 2 群の生カウント表を投げる (typeA=count)。列は「分母の群 → 分子の群」の順に並べる。
#           本家の対比は「2 番目の群 − 最初に出てきた群」なので、分子 − 分母 の向きになる。
#   - ジョブの記録 (ID・投げた日時・引数) を <out-dir>/run/ に残し、同じ <out-dir> で再び実行すると続きから再開する。
# 入力:
#   stats CSV (deseq2 の stats_<project>.csv)、UMI_count CSV、deftable TSV、genome、出力先。
# 出力:
#   <out-dir>/<比較>/<比較>__<ea|page>__<種類>.{tsv,html,log} … WABI が返したもの、そのまま
#   <out-dir>/run/  … 再開用の記録 (jobs.tsv・jobs_history.tsv・run_params.tsv・chipatlas.log) と WABI に渡した入力 (input/)

script_path <- grep("^--file=", commandArgs(), value = TRUE)[[1]]
script_dir <- dirname(normalizePath(sub("^--file=", "", script_path)))
analysis_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/", mustWork = TRUE)

wabi_url <- "https://dtn1.ddbj.nig.ac.jp/wabi/chipatlas/"
wabi_sbatch_options <- "-p epyc -t 180"
supported_genomes <- c("hg38", "hg19", "mm10", "mm9", "rn6", "dm6", "dm3", "ce11", "ce10", "sacCer3")
known_classes <- c("TFs and others", "Histone", "ATAC-Seq", "DNase-seq", "RNA polymerase",
                   "Input control", "Bisulfite-Seq")

log_info <- function(...) cat(sprintf("[INFO] %s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
log_warn <- function(...) cat(sprintf("[WARN] %s %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), sprintf(...)))
fail <- function(...) stop(sprintf(...), call. = FALSE)
now_str <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S %z")

# ---------------------------------------------------------------------------
# 引数
# ---------------------------------------------------------------------------
parse_args <- function(args) {
  o <- list(stats = NULL, umi = NULL, deftable = NULL, genome = NULL, out_dir = NULL,
            comparisons = NULL, classes = "TFs and others,Histone,ATAC-Seq",
            cell_class = "All cell types", threshold = "100",
            distance_up = "5000", distance_down = "5000", fdr = "0.1",
            analyses = "ea,page", max_jobs = "0", timeout_min = "120", retries = "2",
            poll_sec = "60", dry_run = FALSE)
  value_keys <- c("stats", "umi", "deftable", "genome", "out_dir", "comparisons", "classes", "cell_class",
                  "threshold", "distance_up", "distance_down", "fdr", "analyses", "max_jobs",
                  "timeout_min", "retries", "poll_sec")
  i <- 1
  while (i <= length(args)) {
    a <- args[[i]]
    if (identical(a, "--help")) {
      help_file <- file.path(analysis_root, "help", "help_chipatlas.txt")
      cat(readChar(help_file, file.info(help_file)$size, useBytes = TRUE))
      quit(save = "no", status = 0)
    }
    key <- gsub("-", "_", sub("^--", "", a))
    if (identical(a, "--dry-run")) { o$dry_run <- TRUE; i <- i + 1; next }
    if (!startsWith(a, "--") || !(key %in% value_keys)) fail("chipatlas の知らないオプション: %s", a)
    if (i == length(args)) fail("%s に値がありません", a)
    o[[key]] <- args[[i + 1]]
    i <- i + 2
  }
  for (k in c("stats", "umi", "deftable", "genome", "out_dir")) {
    if (is.null(o[[k]])) fail("--stats, --umi, --deftable, --genome, --out-dir は必須です")
  }
  if (!(o$genome %in% supported_genomes)) {
    fail("ChIP-Atlas が対応していないゲノム: %s (対応: %s)", o$genome, paste(supported_genomes, collapse = ", "))
  }
  o$classes <- trimws(strsplit(o$classes, ",", fixed = TRUE)[[1]])
  bad <- setdiff(o$classes, known_classes)
  if (length(bad)) fail("--classes の分類名が違います: %s (使えるもの: %s)", paste(bad, collapse = ", "),
                        paste(known_classes, collapse = " / "))
  o$analyses <- tolower(trimws(strsplit(o$analyses, ",", fixed = TRUE)[[1]]))
  if (!length(o$analyses) || any(!o$analyses %in% c("ea", "page"))) fail("--analyses は ea, page, ea,page のどれか")
  if (!(o$threshold %in% c("50", "100", "200", "500"))) fail("--threshold は 50 / 100 / 200 / 500 のどれか")
  for (k in c("fdr", "max_jobs", "timeout_min", "retries", "poll_sec", "distance_up", "distance_down")) {
    v <- suppressWarnings(as.numeric(o[[k]]))
    if (is.na(v) || v < 0) fail("--%s が数ではありません: %s", gsub("_", "-", k), o[[k]])
    o[[k]] <- v
  }
  if (o$fdr <= 0 || o$fdr >= 1) fail("--fdr は 0 と 1 の間")
  for (k in c("stats", "umi", "deftable")) {
    if (!file.exists(o[[k]])) fail("ファイルがありません: %s", o[[k]])
    o[[k]] <- normalizePath(o[[k]])
  }
  dir.create(o$out_dir, recursive = TRUE, showWarnings = FALSE)
  o$out_dir <- normalizePath(o$out_dir)
  o
}

# ---------------------------------------------------------------------------
# 比較と入力の用意
# ---------------------------------------------------------------------------
normalize_cmp <- function(x) {
  x <- trimws(x)
  x <- sub("^(.+?)_vs_(.+)$", "\\1 / \\2", x)
  parts <- trimws(strsplit(x, "/", fixed = TRUE)[[1]])
  if (length(parts) != 2 || any(parts == "")) fail("比較の書き方が違います: %s (例: hp_r1/cc_r1)", x)
  paste(parts, collapse = " / ")
}

safe_name <- function(x) gsub("[^A-Za-z0-9]+", "_", x)

prepare_comparisons <- function(o) {
  stats <- utils::read.csv(o$stats, check.names = FALSE, stringsAsFactors = FALSE)
  need <- c("aspect", "ens_gene", "log2FoldChange", "padj")
  if (!all(need %in% colnames(stats))) fail("stats CSV に列がありません: %s", paste(setdiff(need, colnames(stats)), collapse = ","))
  stats$aspect <- vapply(stats$aspect, normalize_cmp, "")
  sym <- if ("ext_gene" %in% colnames(stats)) stats$ext_gene else rep(NA_character_, nrow(stats))
  miss <- is.na(sym) | sym == "" | sym == "NA"
  sym[miss] <- stats$ens_gene[miss]
  stats$gene_id <- sym

  avail <- unique(stats$aspect)
  if (is.null(o$comparisons)) {
    cmps <- avail
  } else {
    cmps <- vapply(strsplit(o$comparisons, ",", fixed = TRUE)[[1]], normalize_cmp, "", USE.NAMES = FALSE)
    bad <- setdiff(cmps, avail)
    if (length(bad)) fail("stats CSV に無い比較: %s\n  ある比較: %s", paste(bad, collapse = ", "), paste(avail, collapse = ", "))
  }

  def <- utils::read.delim(o$deftable, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("sample", "group") %in% colnames(def))) fail("deftable に sample / group の列が要ります")
  umi_head <- colnames(utils::read.csv(o$umi, nrows = 1, check.names = FALSE))
  map_group <- function(g) {
    hit <- unique(def$group[tolower(def$group) == tolower(g)])
    if (length(hit) != 1) fail("比較の群 %s が deftable の group に 1 つだけ対応しません (%s)", g, paste(hit, collapse = ","))
    hit
  }

  lapply(cmps, function(cmp) {
    parts <- strsplit(cmp, " / ", fixed = TRUE)[[1]]
    num <- map_group(parts[[1]]); den <- map_group(parts[[2]])
    s <- stats[stats$aspect == cmp, , drop = FALSE]
    deg <- !is.na(s$padj) & s$padj < o$fdr
    up <- s[deg & s$log2FoldChange > 0, , drop = FALSE]
    dn <- s[deg & s$log2FoldChange < 0, , drop = FALSE]
    # 1 つの symbol が両方の向きに出る (別の ens_gene が逆向き) ときは、どちらのリストからも外す
    both <- intersect(up$gene_id, dn$gene_id)
    up_ids <- setdiff(unique(up$gene_id), both)
    dn_ids <- setdiff(unique(dn$gene_id), both)
    samp_den <- unique(def$sample[def$group == den])
    samp_num <- unique(def$sample[def$group == num])
    lost <- setdiff(c(samp_den, samp_num), umi_head)
    if (length(lost)) fail("UMI count CSV に標本の列がありません: %s", paste(lost, collapse = ","))
    list(comparison = cmp, id = safe_name(paste0(parts[[1]], "_vs_", parts[[2]])),
         numerator = num, denominator = den,
         genes_a = up_ids, genes_b = dn_ids, n_both_dropped = length(both),
         n_tested = sum(!is.na(s$padj)),
         samples_den = samp_den, samples_num = samp_num)
  })
}

write_page_input <- function(cmp, umi, path) {
  m <- umi[, c(cmp$samples_den, cmp$samples_num), drop = FALSE]
  # 本家は標本名の末尾 _数字 を落としたものを群名にするので、列名を <群>_<番号> にそろえる
  colnames(m) <- c(paste0(safe_name(cmp$denominator), "_", seq_along(cmp$samples_den)),
                   paste0(safe_name(cmp$numerator), "_", seq_along(cmp$samples_num)))
  gene <- if ("ext_gene" %in% colnames(umi)) umi$ext_gene else umi$ens_gene
  miss <- is.na(gene) | gene == "" | gene == "NA"
  gene[miss] <- umi$ens_gene[miss]
  out <- data.frame(gene = gene, m, check.names = FALSE)
  utils::write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# ---------------------------------------------------------------------------
# ジョブの表
# ---------------------------------------------------------------------------
job_cols <- c("job", "comparison", "comparison_id", "numerator", "denominator", "analysis", "antigen_class",
              "n_genes_a", "n_genes_b", "status", "reason", "attempt", "request_id", "submitted_at",
              "finished_at", "params", "input_a", "input_b", "result_prefix")
result_formats <- c("tsv", "html", "log")

params_string <- function(o, cls, analysis) {
  sprintf("genome=%s;antigenClass=%s;cellClass=%s;threshold=%s;distanceUp=%d;distanceDown=%d;typeA=%s;typeB=%s;fdr=%s",
          o$genome, cls, o$cell_class, o$threshold, as.integer(o$distance_up), as.integer(o$distance_down),
          if (analysis == "ea") "gene" else "count", if (analysis == "ea") "userlist" else "empty", format(o$fdr))
}

build_plan <- function(o, cmps) {
  in_root <- file.path(o$out_dir, "run", "input")
  umi <- NULL
  rows <- list()
  for (cmp in cmps) {
    cdir <- file.path(in_root, cmp$id)
    dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
    fa <- file.path(cdir, sprintf("%s__up_in_%s.txt", cmp$id, safe_name(cmp$numerator)))
    fb <- file.path(cdir, sprintf("%s__up_in_%s.txt", cmp$id, safe_name(cmp$denominator)))
    writeLines(cmp$genes_a, fa); writeLines(cmp$genes_b, fb)
    fp <- file.path(cdir, sprintf("%s__page_counts.tsv", cmp$id))
    if ("page" %in% o$analyses) {
      if (is.null(umi)) umi <- utils::read.csv(o$umi, check.names = FALSE, stringsAsFactors = FALSE)
      write_page_input(cmp, umi, fp)
    }
    for (an in o$analyses) for (cls in o$classes) {
      job <- sprintf("%s__%s__%s", cmp$id, an, safe_name(cls))
      status <- "pending"; reason <- ""
      if (an == "ea" && (length(cmp$genes_a) == 0 || length(cmp$genes_b) == 0)) {
        status <- "skipped"
        reason <- sprintf("EA は 2 つのリストが互いの背景なので、片方が 0 個だと計算できない (%s で高い %d 個、%s で高い %d 個)",
                          cmp$numerator, length(cmp$genes_a), cmp$denominator, length(cmp$genes_b))
      }
      if (an == "page" && (length(cmp$samples_den) < 2 || length(cmp$samples_num) < 2)) {
        status <- "skipped"; reason <- "PAGE は各群に 2 標本以上が要る"
      }
      rows[[length(rows) + 1]] <- data.frame(
        job = job, comparison = cmp$comparison, comparison_id = cmp$id, numerator = cmp$numerator,
        denominator = cmp$denominator, analysis = toupper(an), antigen_class = cls,
        n_genes_a = if (an == "ea") length(cmp$genes_a) else NA,
        n_genes_b = if (an == "ea") length(cmp$genes_b) else NA,
        status = status, reason = reason, attempt = 0L, request_id = "", submitted_at = "", finished_at = "",
        params = params_string(o, cls, an),
        input_a = if (an == "ea") fa else fp, input_b = if (an == "ea") fb else "",
        result_prefix = file.path(o$out_dir, cmp$id, job),
        stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

state_path <- function(o) file.path(o$out_dir, "run", "jobs.tsv")
history_path <- function(o) file.path(o$out_dir, "run", "jobs_history.tsv")

save_state <- function(o, jobs) {
  p <- state_path(o)
  tmp <- paste0(p, ".tmp")
  utils::write.table(jobs[, job_cols], tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  file.rename(tmp, p)
}

append_history <- function(o, j, event, note = "") {
  p <- history_path(o)
  line <- data.frame(time = now_str(), job = j$job, event = event, attempt = j$attempt,
                     request_id = j$request_id, params = j$params, note = gsub("[\t\n]", " ", note))
  utils::write.table(line, p, sep = "\t", quote = FALSE, row.names = FALSE,
                     col.names = !file.exists(p), append = file.exists(p))
}

load_or_init_state <- function(o, plan) {
  p <- state_path(o)
  dir.create(dirname(p), recursive = TRUE, showWarnings = FALSE)
  if (!file.exists(p)) return(plan)
  old <- utils::read.delim(p, check.names = FALSE, stringsAsFactors = FALSE, colClasses = "character", na.strings = "")
  old[is.na(old)] <- ""
  log_info("前の記録 %s を読んで再開します (%d 件)", p, nrow(old))
  for (i in seq_len(nrow(plan))) {
    k <- match(plan$job[[i]], old$job)
    if (is.na(k)) next
    if (!identical(old$params[[k]], plan$params[[i]])) {
      fail("同じ <out-dir> に別の引数の記録があります (%s)。別の --out-dir を使ってください", plan$job[[i]])
    }
    if (plan$status[[i]] == "skipped") next
    st <- old$status[[k]]
    keep <- c("status", "reason", "attempt", "request_id", "submitted_at", "finished_at")
    if (st %in% c("done", "submitted")) {
      for (c in keep) plan[[c]][[i]] <- old[[c]][[k]]
      # 終わったのに結果のファイルが無いもの (消した・置き場所が変わった) は、WABI から取り直す
      have <- file.exists(paste0(plan$result_prefix[[i]], ".", result_formats))
      if (st == "done" && !all(have)) plan$status[[i]] <- "submitted"
    } else if (st %in% c("failed", "timeout")) {
      # 前に失敗したものは、再試行の回数を戻して投げ直す
      plan$attempt[[i]] <- 0L
      plan$reason[[i]] <- paste("前回:", old$reason[[k]])
    }
  }
  plan$attempt <- as.integer(plan$attempt)
  plan
}

# ---------------------------------------------------------------------------
# WABI
# ---------------------------------------------------------------------------
curl_text <- function(args) {
  out <- suppressWarnings(tryCatch(system2("curl", args, stdout = TRUE, stderr = TRUE),
                                   error = function(e) paste("curl error:", conditionMessage(e))))
  paste(out, collapse = "\n")
}

submit_job <- function(o, j) {
  cls <- j$antigen_class
  common <- c("-sS", "--max-time", "120", "-X", "POST",
              "-d", "format=text", "-d", "result=www",
              "-d", paste0("genome=", o$genome),
              "-d", paste0("antigenClass=", cls),
              "-d", paste0("cellClass=", o$cell_class),
              "-d", paste0("threshold=", o$threshold),
              "-d", "permTime=1",
              "-d", paste0("title=", j$job),
              "-d", paste0("distanceUp=", as.integer(o$distance_up)),
              "-d", paste0("distanceDown=", as.integer(o$distance_down)),
              "-d", paste0("sbatchOptions=", wabi_sbatch_options))
  if (j$analysis == "EA") {
    body <- c("-d", "typeA=gene", "--data-urlencode", paste0("bedAFile@", j$input_a),
              "-d", "typeB=userlist", "--data-urlencode", paste0("bedBFile@", j$input_b),
              "-d", paste0("descriptionA=up_in_", safe_name(j$numerator)),
              "-d", paste0("descriptionB=up_in_", safe_name(j$denominator)))
  } else {
    body <- c("-d", "typeA=count", "--data-urlencode", paste0("bedAFile@", j$input_a),
              "-d", "typeB=empty", "-d", "bedBFile=empty",
              "-d", "descriptionA=empty", "-d", "descriptionB=empty")
  }
  args <- shQuote(c(common, body, wabi_url))
  resp <- curl_text(args)
  id <- regmatches(resp, regexpr("wabi_chipatlas_[0-9A-Za-z._-]+", resp))
  if (!length(id)) return(list(ok = FALSE, msg = substr(gsub("\n", " ", resp), 1, 300)))
  list(ok = TRUE, id = id[[1]])
}

job_status <- function(id) {
  resp <- curl_text(shQuote(c("-sS", "--max-time", "90", paste0(wabi_url, id))))
  st <- regmatches(resp, regexpr("(?m)^status:\\s*\\S+", resp, perl = TRUE))
  if (!length(st)) return(list(status = "unknown", raw = resp))
  list(status = tolower(sub("^status:\\s*", "", st[[1]])), raw = resp)
}

# WABI が返す tsv / html / log をそのまま保存する (中身は加工しない)
download_result <- function(j) {
  dir.create(dirname(j$result_prefix), recursive = TRUE, showWarnings = FALSE)
  base <- paste0(wabi_url, j$request_id, "?info=result&format=")
  for (f in result_formats) {
    dest <- paste0(j$result_prefix, ".", f)
    tmp <- paste0(dest, ".part")
    rc <- system2("curl", shQuote(c("-fsS", "--max-time", "600", "-o", tmp, paste0(base, f))))
    if (rc != 0 || !file.exists(tmp) || file.size(tmp) == 0) {
      unlink(tmp)
      return(list(ok = FALSE, msg = sprintf("%s を取れなかった (curl %d)", f, rc)))
    }
    file.rename(tmp, dest)
  }
  list(ok = TRUE)
}

run_jobs <- function(o, jobs) {
  max_attempts <- 1L + as.integer(o$retries)
  timeout_sec <- o$timeout_min * 60
  todo <- function() which(jobs$status %in% c("pending", "submitted"))
  last_progress <- ""
  while (length(todo())) {
    # 1) 状態を見る
    for (i in which(jobs$status == "submitted")) {
      j <- jobs[i, ]
      st <- job_status(j$request_id)
      if (st$status == "finished") {
        dl <- download_result(j)
        if (isTRUE(dl$ok)) {
          jobs$status[i] <- "done"; jobs$reason[i] <- ""
          if (!nzchar(jobs$finished_at[i])) jobs$finished_at[i] <- now_str()
          append_history(o, jobs[i, ], "done")
          log_info("終わり: %s (%s)", j$job, j$request_id)
        } else {
          jobs$status[i] <- "retry"; jobs$reason[i] <- dl$msg
          append_history(o, jobs[i, ], "bad_result", dl$msg)
          log_warn("結果を取れない: %s (%s): %s", j$job, j$request_id, dl$msg)
        }
      } else if (grepl("error|fail|cancel|timeout|not.?found", st$status)) {
        jobs$status[i] <- "retry"; jobs$reason[i] <- paste("WABI の状態:", st$status)
        append_history(o, jobs[i, ], "wabi_failed", st$status)
        log_warn("WABI で失敗: %s (%s): %s", j$job, j$request_id, st$status)
      } else {
        sub_t <- as.numeric(as.POSIXct(j$submitted_at, format = "%Y-%m-%d %H:%M:%S"))
        if (!is.na(sub_t) && as.numeric(Sys.time()) - sub_t > timeout_sec) {
          jobs$status[i] <- "retry"
          jobs$reason[i] <- sprintf("%g 分で終わらなかった (最後の状態: %s)", o$timeout_min, st$status)
          append_history(o, jobs[i, ], "timeout", st$status)
          log_warn("時間切れ: %s (%s)", j$job, j$request_id)
        }
      }
    }
    # 2) 再試行の回数を使い切ったものを止める
    for (i in which(jobs$status == "retry")) {
      if (jobs$attempt[i] >= max_attempts) {
        jobs$status[i] <- if (grepl("終わらなかった", jobs$reason[i])) "timeout" else "failed"
        append_history(o, jobs[i, ], jobs$status[i], jobs$reason[i])
      } else {
        jobs$status[i] <- "pending"
      }
    }
    # 3) 投げる (--max-jobs 0 なら全部を一度に)
    running <- sum(jobs$status == "submitted")
    for (i in which(jobs$status == "pending")) {
      if (o$max_jobs > 0 && running >= o$max_jobs) break
      jobs$attempt[i] <- jobs$attempt[i] + 1L
      r <- submit_job(o, jobs[i, ])
      if (isTRUE(r$ok)) {
        jobs$status[i] <- "submitted"; jobs$request_id[i] <- r$id; jobs$submitted_at[i] <- now_str()
        append_history(o, jobs[i, ], "submitted")
        log_info("投げた: %s → %s (試行 %d)", jobs$job[i], r$id, jobs$attempt[i])
        running <- running + 1
      } else {
        jobs$reason[i] <- paste("投げられなかった:", r$msg)
        append_history(o, jobs[i, ], "submit_failed", r$msg)
        log_warn("投げられなかった: %s: %s", jobs$job[i], r$msg)
        jobs$status[i] <- if (jobs$attempt[i] >= max_attempts) "failed" else "pending"
      }
      Sys.sleep(1)
    }
    save_state(o, jobs)
    tab <- table(factor(jobs$status, levels = c("pending", "submitted", "done", "failed", "timeout", "skipped")))
    prog <- paste(sprintf("%s=%d", names(tab), as.integer(tab)), collapse = " ")
    if (!identical(prog, last_progress)) { log_info("状態: %s", prog); last_progress <- prog }
    if (length(todo())) Sys.sleep(o$poll_sec)
  }
  jobs
}

main <- function() {
  o <- parse_args(commandArgs(trailingOnly = TRUE))
  cmps <- prepare_comparisons(o)
  for (c in cmps) {
    log_info("比較 %s: %s で高い DEG %d 個、%s で高い DEG %d 個 (FDR < %g)、PAGE の標本 %s %d / %s %d",
             c$comparison, c$numerator, length(c$genes_a), c$denominator, length(c$genes_b), o$fdr,
             c$denominator, length(c$samples_den), c$numerator, length(c$samples_num))
  }
  plan <- build_plan(o, cmps)
  jobs <- load_or_init_state(o, plan)
  save_state(o, jobs)
  utils::write.table(data.frame(key = c("run_started", "stats", "umi", "deftable", "genome", "classes",
                                        "cell_class", "threshold", "distance_up", "distance_down", "fdr",
                                        "wabi_url", "sbatch_options"),
                                value = c(now_str(), o$stats, o$umi, o$deftable, o$genome,
                                          paste(o$classes, collapse = ","), o$cell_class, o$threshold,
                                          o$distance_up, o$distance_down, o$fdr, wabi_url, wabi_sbatch_options)),
                     file.path(o$out_dir, "run", "run_params.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (o$dry_run) {
    log_info("--dry-run: 入力と計画だけを書きました: %s", state_path(o))
    print(table(jobs$analysis, jobs$status))
    return(invisible())
  }
  jobs <- run_jobs(o, jobs)
  bad <- jobs$status %in% c("failed", "timeout")
  if (any(bad)) {
    log_warn("失敗・時間切れが %d 件あります。同じ引数でもう一度実行すると、その分だけ投げ直します", sum(bad))
    quit(save = "no", status = 1)
  }
  log_info("完了: %s", o$out_dir)
}

main()
