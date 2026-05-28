#!/usr/local/bin/Rscript

# 役割:
#   pic の集計結果をもとに、read 数と UMI 数の関係を図にする。
# 入力:
#   args[1]: pic が生成した集計 TSV。
#   args[2]: プロットの出力先ディレクトリ。
# 出力:
#   指定したディレクトリ内の `sim_umis.png`。

suppressWarnings(library(ggplot2))
suppressWarnings(library(scales))
suppressWarnings(library(rlang))
suppressWarnings(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
mapping_sum <- file.path(args[1])
plot_dir <- file.path(args[2])

# mapping_sum <- "/Users/zou/tmp/NextSeq_250514/mapping_liver/simulation/mapping_sum.tsv"
# plot_dir <- "/Users/zou/tmp/plot_test"


if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
# depth と UMI 数の関係を見るのに必要な列だけを残す。
df <- read.table(mapping_sum, header = TRUE, check.names = FALSE)[, c(1, 2, 3, 10)] %>%
  mutate(
    across(c("total", "umis"), as.numeric),
    group = factor(group, levels = unique(group)) # groupの順序を固定
  )

max_umis <- max(df$umis)/1000
for (i in c(1,10,100,1000)) {
  # 観測された UMI 数に応じて、見やすい上限値を決める。
  if (i > max_umis) {
    max_y <- i * 10
    break
  }
}

group_colors <- c(
  "#4F71BE", "#DE8344", "#A5A5A5", "#F5C242",
  "#B02418", "#7EAB55", "#68349A", "#68349A",
  "#583327", "#B08FB0"
)

plot_all <- ggplot() +
  # 実測（実線）
  geom_line(
    data = df,
    aes(x = total/1000, y = umis/1000, group = sample, color = group),
    linewidth = 0.5,
    linetype = "solid"
  ) +

  scale_color_manual(values = group_colors) +

  scale_y_log10(
    breaks = scales::log_breaks(base = 10),
    labels = scales::label_comma(),
    limits = c(min(df$umis)/1000, max_y)
  ) +

  # x軸：log10(total)の整数
  scale_x_log10(breaks = c(10,100,1000,10000),
                minor_breaks = c(20,50,200,500,2000,5000,20000),
                labels = scales::label_comma(),
                limits = c(10, 20000)
                ) +
  labs(
    x = expression("Total Reads (K)"),
    y = "UMIs (K)"
  ) +
  theme_bw(base_size = 13)

ggsave(paste(plot_dir, "sim_umis.png", sep = "/"), plot = plot_all, width = 10, height = 10, dpi = 300)
