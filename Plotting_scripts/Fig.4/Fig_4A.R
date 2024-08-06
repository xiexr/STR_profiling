library(ggplot2)
rm(list=ls())
library(ggupset)
library(tidyverse)
library(ComplexUpset)
# 你的数据


setwd("E:/Rst/each_pop")
example = read.csv("output_Rst_0.25_upset_plot.csv",header=TRUE,check.names = FALSE)

# 将数据转换为长格式
data_long <- example %>%
  pivot_longer(cols = -strname, names_to = "category", values_to = "value") %>%
  filter(value == 1) %>%
  group_by(strname) %>%
  summarise(categories = list(category))

p <- ggplot(data_long,aes(x=categories)) +
  geom_bar(fill='black') +
  scale_x_upset(order_by='degree',scale_name = "Genres")+
  theme_combmatrix(combmatrix.panel.point.color.fill = "blue",
                   combmatrix.panel.line.color.fill = "blue",
                   combmatrix.panel.line.size = 1,
                   combmatrix.label.make_space = FALSE)+
  scale_y_continuous(expand = expansion(mult = c(0, 0.2)))+
  theme_classic()+
  theme(plot.margin = margin(5, 5, 5, 40))

ggsave("my_upset_plot.png", plot = p, width = 12, height = 8, dpi = 600)