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
