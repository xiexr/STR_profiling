library(ggplot2)
library(scales)
df <- read.csv('majorallelefreqof16bp.txt', sep='\t')

# 提取每个length对应的majorallelefreq生成向量
bp1 <- df$majorallelefreq[df$length == 1]
bp2 <- df$majorallelefreq[df$length == 2]
bp3 <- df$majorallelefreq[df$length == 3]
bp4 <- df$majorallelefreq[df$length == 4]
bp5 <- df$majorallelefreq[df$length == 5]
bp6 <- df$majorallelefreq[df$length == 6]

library(ggplot2)
library(gridExtra)

# 创建直方图函数
plot_histogram <- function(data, title) {
  ggplot(data.frame(data), aes(x = data)) +
    geom_histogram(binwidth = .01, fill = "#34495e") +
    labs(x = "Major allele frequency", y = "# of STR") +
    ggtitle(title) +
    theme_classic()
}

# 绘制组合直方图
p1 <- plot_histogram(bp1, "1bp")
p2 <- plot_histogram(bp2, "2bp")
p3 <- plot_histogram(bp3, "3bp")
p4 <- plot_histogram(bp4, "4bp")
p5 <- plot_histogram(bp5, "4bp")
p6 <- plot_histogram(bp6, "6bp")

plot_grid = grid.arrange(p1, p2, p3, p4,p5, p6, nrow = 2, ncol = 3)

# 保存图片为png格式
ggsave("eachbp16.png", dpi=600,plot_grid, width = 10, height = 8, units = "in")
ggsave("eachbp16.pdf", plot_grid, width = 10, height = 8, units = "in")

