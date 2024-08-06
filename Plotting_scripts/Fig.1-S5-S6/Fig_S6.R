library(ggplot2)
library(scales)
library(gridExtra)
library(cowplot)
rm(list = ls())
setwd('E:/experiment/STR2/01.STR_summary/major_allele_diff_from_ref')
df <- read.csv('majoralleledifffromref_resort.txt', sep='\t')

count_dataframes <- list()
for (i in 1:6) {
  bp <- df$majordifffromref[df$periodlen == i]
  count_data <- table(bp)
  count_dataframe <- as.data.frame(count_data)
  count_dataframe$bp <- as.factor(count_dataframe$bp)
  assign(paste0("bp", i), count_dataframe)
}

#count_dataframe$Var1 <- as.factor(count_dataframe$Var1)
plot_count_data <- function(count_dataframe,title) {
  plot <- ggplot(count_dataframe, aes(x = bp, y = Freq)) +
    geom_col(fill = "#34495e") +
    scale_y_log10(
      expand = c(0, 0),
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    # scale_x_continuous(breaks = seq(-15, 15, 5)) +
    ggtitle(title) +
    theme_classic()+
    theme(axis.text.x = element_text(size = 8, angle = 90,vjust=0.5))+  # 设置X轴标题字体大小
    labs(x = "", y = "") 

  return(plot)
}
# 存储bp和标题的列表
bp_list <- list(bp1, bp2, bp3, bp4, bp5, bp6)
title_list <- c("1bp", "2bp", "3bp", "4bp", "5bp", "6bp")

# 创建一个空的图表列表
plot_list <- list()

# 循环遍历并生成图表
for (i in 1:length(bp_list)) {
  plot <- plot_count_data(bp_list[[i]], title_list[i])
  plot_list[[i]] <- plot
}

# 绘制2x3的柱形图
plot_grid <- do.call(grid.arrange, c(plot_list, nrow = 2, ncol = 3))#绘制组合直方图

ggsave('16bpmajoralleledifffromref_resort.png',plot_grid,dpi=600)
ggsave('16bpmajoralleledifffromref_resort.pdf',plot_grid)

