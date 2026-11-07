# 役割: レポートに埋め込む JSON を速く書く (report_page.R の render_report_page が使う)。
# 注記: report_build.R を責務別に分割したファイル。cmd_build_report.R が report_build.R 経由で source する。

# レポートに埋め込む JSON を速く作る。jsonlite::toJSON と 1 バイトも違わない出力にする。
# 遅いのは「長さ 1 の値を as.list で並べたリスト」(plotly の x / y / text、発現の counts など) を
# jsonlite が要素ごとに S4 の振り分けを通して書くこと。同じ並びをベクトルにまとめて 1 回で書く。
#   - 名前の無いリストで、要素がすべて属性の無い長さ 1 の同じ型の値 → I(ベクトル) (auto_unbox でも配列のまま)
#   - 名前の無いリストで、要素がすべて「名前の無い、長さ k の、長さ 1 の値のリスト」で、列ごとに型が揃う
#     → 行ごとの配列の配列。data.frame を dataframe = "values" で書いた JSON を、そのまま埋め込む (class "json")
# それ以外は中へ降りて同じことをする。
pic_json_scalar_type <- function(e) {
  if (is.atomic(e) && length(e) == 1L && is.null(attributes(e))) typeof(e) else NA_character_
}

pic_json_compact <- function(x, digits) {
  if (!is.list(x) || is.object(x) || length(x) == 0L) return(x)
  if (is.null(names(x))) {
    types <- vapply(x, pic_json_scalar_type, character(1))
    if (!anyNA(types) && all(types == types[[1]]) && types[[1]] %in% c("double", "integer", "character", "logical")) {
      return(I(unlist(x, use.names = FALSE)))
    }
    rows <- pic_json_rows(x, digits)
    if (!is.null(rows)) return(rows)
  }
  nm <- names(x)
  x <- lapply(x, pic_json_compact, digits = digits)
  names(x) <- nm
  x
}

pic_json_rows <- function(x, digits) {
  k <- length(x[[1]])
  if (!is.list(x[[1]]) || k == 0L) return(NULL)
  ok <- vapply(x, function(e) is.list(e) && !is.object(e) && is.null(names(e)) && length(e) == k, logical(1))
  if (!all(ok)) return(NULL)
  cols <- vector("list", k)
  for (j in seq_len(k)) {
    col <- lapply(x, `[[`, j)
    types <- vapply(col, pic_json_scalar_type, character(1))
    if (anyNA(types) || any(types != types[[1]]) || !(types[[1]] %in% c("double", "integer", "character", "logical"))) return(NULL)
    cols[[j]] <- unlist(col, use.names = FALSE)
  }
  df <- as.data.frame(stats::setNames(cols, paste0("V", seq_len(k))), stringsAsFactors = FALSE)
  structure(as.character(jsonlite::toJSON(df, dataframe = "values", na = "null", digits = digits)), class = "json")
}

pic_to_json <- function(x, digits = 4) {
  jsonlite::toJSON(pic_json_compact(x, digits), auto_unbox = TRUE, null = "null", na = "null", digits = digits,
                   json_verbatim = TRUE)
}

# --- 埋め込みデータの圧縮 (render_report_page が使う) ----------------------------
# 大きな埋め込み (plotly の js・PIC_PLOTS/EXPR/PCA の JSON・大きな表の HTML・ダウンロード用ファイル) は
# zlib (RFC1950。R の memCompress(type = "gzip") が出す形) で縮めて base64 にし、
# 実行されない <script type="application/x-pic-z" id="..."> に入れる。
# ブラウザは report_data.js が DecompressionStream("deflate") で展開する。
# JS の文字列として読ませないので、開いたときの JS の構文解析も軽くなる。

# 文字列 (UTF-8) または raw を zlib で縮めて base64 にする。
pic_z_b64 <- function(x) {
  raw <- if (is.raw(x)) x else charToRaw(enc2utf8(paste(as.character(x), collapse = "")))
  jsonlite::base64_enc(memCompress(raw, "gzip"))
}

# 実行されない script 要素に base64 を入れる。attrs は名前付きの文字列 (data-n など)。
pic_z_block <- function(id, b64, attrs = NULL) {
  a <- if (length(attrs) > 0) paste0(" ", names(attrs), '="', vapply(attrs, html_escape, character(1)), '"', collapse = "") else ""
  paste0('<script type="application/x-pic-z" id="', id, '"', a, '>', b64, '</script>')
}

# body の HTML から、大きな <table> (入れ子を含まない一番内側の表で min_bytes 以上) を
# <template class="pic-zfrag" data-z="<id>"></template> に置き換え、表の HTML を縮めた塊にする。
# report_data.js が展開して同じ場所に戻してから report.js の初期化を走らせるので、表示は変わらない。
# 返り値: list(html = 置き換えた body, blocks = script 要素の文字列のベクトル)
pic_z_fragments <- function(html, min_bytes = 65536L, prefix = "pic-z-frag-") {
  html <- enc2utf8(html)
  op <- gregexpr("<table[ >]", html, perl = TRUE, useBytes = TRUE)[[1]]
  cl <- gregexpr("</table>", html, fixed = TRUE, useBytes = TRUE)[[1]]
  if (op[[1]] < 0 || cl[[1]] < 0) return(list(html = html, blocks = character(0)))
  r <- charToRaw(html)
  cut_a <- integer(0); cut_b <- integer(0)
  for (o in op) {
    c1 <- cl[cl > o]
    if (length(c1) == 0) next
    e <- c1[[1]] + 8L - 1L                       # "</table>" の最後のバイト
    if (any(op > o & op < c1[[1]])) next        # 入れ子の外側は置き換えない
    if (e - o + 1L < min_bytes) next
    if (length(cut_b) > 0 && o <= cut_b[[length(cut_b)]]) next
    cut_a <- c(cut_a, o); cut_b <- c(cut_b, e)
  }
  if (length(cut_a) == 0) return(list(html = html, blocks = character(0)))
  parts <- character(0); blocks <- character(0); pos <- 1L
  for (k in seq_along(cut_a)) {
    frag <- r[cut_a[[k]]:cut_b[[k]]]
    s <- rawToChar(frag)
    if (grepl("<script", s, fixed = TRUE, useBytes = TRUE)) next   # script を含む表は innerHTML で動きが変わるので残す
    id <- paste0(prefix, k)
    if (cut_a[[k]] > pos) parts <- c(parts, rawToChar(r[pos:(cut_a[[k]] - 1L)]))
    parts <- c(parts, sprintf('<template class="pic-zfrag" data-z="%s"></template>', id))
    blocks <- c(blocks, pic_z_block(id, pic_z_b64(frag)))
    pos <- cut_b[[k]] + 1L
  }
  if (pos <= length(r)) parts <- c(parts, rawToChar(r[pos:length(r)]))
  out <- paste(parts, collapse = "")
  Encoding(out) <- "UTF-8"
  list(html = out, blocks = blocks)
}
