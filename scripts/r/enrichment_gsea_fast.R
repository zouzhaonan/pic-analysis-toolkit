# 役割: enrich の GSEA で、DOSE の leading edge の計算を同じ値の速い版に置き換える。
#   無効にするには PIC_GSEA_FAST=0。

# GSEA の leading edge を速く求める (DOSE 4.0.0 の gseaScores / leading_edge と同じ値)。
# clusterProfiler::GSEA は fgsea のあと、すべての遺伝子セットについて gseaScores で長さ N (全遺伝子)
# の running score を data.frame で作り、leading_edge で読む。pvalueCutoff = 1 では全セットなので、
# GO_BP では GSEA 1 本の半分近くがここになる。
# 同じ値を、セットのヒットの位置と「その手前の位置」だけで求める:
#   - Phit は R の cumsum (long double で足す) なので、0 を挟んでも挟まなくても同じ値になる。
#   - Pmiss は 1/(N - Nh) を m 回 cumsum した値。Nh ごとに 1 回だけ cumsum を作って引く。
#   - running score は ヒットで上がり、ミスで下がるので、最大はヒットか先頭、最小はヒットの直前・末尾・先頭に
#     しか来ない。which.max / which.min と同じ「最初の位置」をその候補から選ぶ。
# DOSE の版が違うときは置き換えない (pic_gsea_fast_install が FALSE を返す)。
pic_gsea_fast_env <- new.env(parent = emptyenv())

pic_gsea_pos_lookup <- function(gene_names) {
  e <- pic_gsea_fast_env
  if (!is.null(e$names) && identical(e$names, gene_names)) return(e$pos)
  pos <- new.env(hash = TRUE, parent = emptyenv(), size = length(gene_names))
  list2env(stats::setNames(as.list(seq_along(gene_names)), gene_names), envir = pos)
  e$names <- gene_names
  e$pos <- pos
  e$miss_cumsum <- list()
  pos
}

pic_gsea_miss_cumsum <- function(N, Nh) {
  key <- as.character(Nh)
  s <- pic_gsea_fast_env$miss_cumsum[[key]]
  if (is.null(s)) {
    # gseaScores と同じ: Pmiss[!hits] <- 1/(N - Nh); cumsum(Pmiss)。ミスの数だけ足した値の列。
    s <- cumsum(rep(1 / (N - Nh), N - Nh))
    pic_gsea_fast_env$miss_cumsum[[key]] <- s
  }
  s
}

pic_gsea_scores_fast <- function(geneList, geneSet, exponent = 1, fortify = FALSE) {
  if (isTRUE(fortify)) return(pic_gsea_fast_env$orig_gseaScores(geneList, geneSet, exponent, fortify))
  nm <- names(geneList)
  N <- length(geneList)
  pos <- pic_gsea_pos_lookup(nm)
  gs <- unique(as.character(geneSet))
  hp <- unlist(mget(gs, envir = pos, ifnotfound = list(NA_integer_)), use.names = FALSE)
  hp <- sort(hp[!is.na(hp)])
  Nh <- length(hp)
  w <- abs(geneList[hp])^exponent
  # sum(Phit) と cumsum(Phit / NR): 0 の要素は足しても変わらない (long double で足す)。
  NR <- sum(w)
  phit <- cumsum(w / NR)
  S <- c(0, pic_gsea_miss_cumsum(N, Nh))
  j <- seq_len(Nh)
  hv <- phit - S[hp - j + 1L]                       # ヒットの位置の running score
  # 最大・最小の候補: 先頭 (1)、ヒット、ヒットの直前のミス、末尾 (N)
  prev_phit <- c(0, phit)
  cand_pos <- c(1L, hp, hp - 1L, N)
  cand_val <- c(
    if (Nh > 0 && hp[[1]] == 1L) hv[[1]] else 0 - S[2L],
    hv,
    prev_phit[j] - S[pmax(hp - 1L - (j - 1L), 0L) + 1L],
    phit[Nh] - S[N - Nh + 1L]
  )
  # ミスでない候補 (ヒットの直前がヒット、先頭/末尾がヒット) の値は hv に入っているので、ヒットの値に差し替える
  is_hit <- cand_pos %in% hp
  cand_val[is_hit] <- hv[match(cand_pos[is_hit], hp)]
  ok <- cand_pos >= 1L
  cand_pos <- cand_pos[ok]; cand_val <- cand_val[ok]
  o <- order(cand_pos)
  cand_pos <- cand_pos[o]; cand_val <- cand_val[o]
  dup <- duplicated(cand_pos)
  cand_pos <- cand_pos[!dup]; cand_val <- cand_val[!dup]
  max_ES <- max(cand_val); min_ES <- min(cand_val)
  ES <- if (abs(max_ES) > abs(min_ES)) max_ES else min_ES
  list(ES = ES, pic_fast = list(
    hv = hv, genes = nm[hp], N = N, Nh = Nh,
    pos_max = cand_pos[which.max(cand_val)], pos_min = cand_pos[which.min(cand_val)]
  ))
}

pic_leading_edge_fast <- function(observed_info) {
  if (length(observed_info) == 0 || is.null(observed_info[[1]]$pic_fast)) {
    return(pic_gsea_fast_env$orig_leading_edge(observed_info))
  }
  core_enrichment <- lapply(observed_info, function(x) {
    f <- x$pic_fast
    if (x$ES >= 0) {
      f$genes[1:which.max(f$hv)]
    } else {
      i <- which.min(f$hv)
      f$genes[-c(1:(i - 1))]                      # DOSE と同じ (i = 1 のときは先頭が落ちる)
    }
  })
  rank <- sapply(observed_info, function(x) {
    f <- x$pic_fast
    if (x$ES >= 0) f$pos_max else f$N - f$pos_min + 1
  })
  tags <- sapply(observed_info, function(x) {
    f <- x$pic_fast
    if (x$ES >= 0) which.max(f$hv) / f$Nh else (f$Nh - which.min(f$hv) + 1) / f$Nh
  })
  ll <- sapply(observed_info, function(x) {
    f <- x$pic_fast
    if (x$ES >= 0) f$pos_max / f$N else (f$N - f$pos_min + 1) / f$N
  })
  N <- observed_info[[1]]$pic_fast$N
  setSize <- sapply(observed_info, function(x) x$pic_fast$Nh)
  signal <- tags * (1 - ll) * (N/(N - setSize))
  tags <- paste0(round(tags * 100), "%")
  ll <- paste0(round(ll * 100), "%")
  signal <- paste0(round(signal * 100), "%")
  leading_edge <- paste0("tags=", tags, ", list=", ll, ", signal=", signal)
  list(rank = rank, tags = tags, list = ll, signal = signal,
       leading_edge = leading_edge, core_enrichment = core_enrichment)
}

pic_gsea_fast_install <- function() {
  if (isTRUE(pic_gsea_fast_env$installed)) return(TRUE)
  if (!identical(as.character(utils::packageVersion("DOSE")), "4.0.0")) return(FALSE)
  if (identical(Sys.getenv("PIC_GSEA_FAST", "1"), "0")) return(FALSE)
  ns <- asNamespace("DOSE")
  pic_gsea_fast_env$orig_gseaScores <- get("gseaScores", envir = ns)
  pic_gsea_fast_env$orig_leading_edge <- get("leading_edge", envir = ns)
  utils::assignInNamespace("gseaScores", pic_gsea_scores_fast, ns = "DOSE")
  utils::assignInNamespace("leading_edge", pic_leading_edge_fast, ns = "DOSE")
  pic_gsea_fast_env$installed <- TRUE
  TRUE
}
