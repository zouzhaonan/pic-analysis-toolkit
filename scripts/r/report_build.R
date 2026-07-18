# 役割:
#   解析出力ディレクトリ (mapping / deseq2 / enrich) を読み取り、
#   自己完結型の HTML レポートを 1 ファイル生成する補助関数群。
# 入力:
#   <out-dir> 配下の mapping_sum__*.tsv, deseq2/**/stats_*.csv 等。
# 出力:
#   <out-dir>/report_<project>.html (plotly + 画像を base64 で内包)。
# 注記:
#   PCA / heatmap / MA / volcano は plotly でインタラクティブ (hover でラベル)。
#   mapping QC は CSS で 100% 積み上げ棒 + データバー表。enrich plot は再生成して内包。
# 実体は責務別に分割済み。以下を順に source して全関数を読み込む。
# (source は globalenv で評価されるため script_dir が参照可能)
for (.pic_report_part in c(
  "report_helpers.R", "report_qc_dim.R", "report_deg_expr.R",
  "report_enrichment.R", "report_page.R")) {
  source(file.path(script_dir, .pic_report_part))
}
rm(.pic_report_part)

