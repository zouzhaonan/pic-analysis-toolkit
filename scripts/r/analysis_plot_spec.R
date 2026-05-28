# 役割:
#   R 解析スクリプトで使う定数を一元管理する。

pic_plot_spec <- function() {
  list(
    defaults = list(
      fdr = 0.1,
      ntop = 500L,
      threads = 8L,
      enrich_methods_csv = "GO_BP,KEGG,REACTOME,WIKIPATHWAYS,HDO,HPO,MPO"
    ),
    plot = list(
      palette_d3 = "category20",
      png_width = 10,
      png_height = 8,
      contrast = list(
        label_top_n = 10L,
        point_size = 1.2,
        point_alpha = 0.7,
        direction_colors = c("#d7301f", "#2166ac", "#bdbdbd"),
        guide_line_color = "#636363",
        deg_count_height_min = 8,
        deg_count_height_per_contrast = 0.8,
        heatmap = list(
          width = 12,
          height = 10,
          rowname_max = 80L,
          fontsize_row = 7,
          fontsize_col = 10,
          color_low = "#2166ac",
          color_mid = "#ffffff",
          color_high = "#b2182b",
          color_steps = 100L
        )
      ),
      violin = list(
        facet_nrow = 3L,
        facet_ncol = 4L,
        panels_per_page = 12L,
        pdf_width = 15,
        pdf_height = 9,
        box_width = 0.6,
        box_alpha = 0.7,
        box_linewidth = 0.4,
        point_size = 1.8,
        point_alpha = 0.7,
        legend_title_size = 10,
        legend_text_size = 9,
        empty_panel_prefix = "__empty_"
      ),
      dimension = list(
        pca_point_size = 4,
        pca_label_size = 4,
        corr_color_low = "#2166ac",
        corr_color_mid = "#ffffff",
        corr_color_high = "#d7301f",
        corr_width_min = 8,
        corr_width_per_sample = 1.2,
        corr_height_min = 7,
        corr_height_per_sample = 1.0
      ),
      enrichment = list(
        top_n_per_direction = 20L,
        point_alpha = 0.95,
        term_label_max_chars = 80L,
        gradient_low = "#FF0000",
        gradient_high = "#0000FF",
        save_width = 10,
        save_height = 10,
        save_height_min = 8,
        save_height_per_label = 0.15,
        save_height_max = 28
      )
    )
  )
}
