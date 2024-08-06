library(ggplot2)
library(patchwork)
library(reshape2)
rm(list = ls())
setwd('E:heatmap')
# Create example data
heatmap_data <- read.csv('output(old).csv',sep = '\t')

heat_data=melt(heatmap_data,id.vars = 'X')
heat_data$variable <- factor(heat_data$variable, levels = c("X6bp", "X5bp", "X4bp", "X3bp", "X2bp", "X1bp"))
heatmap_plot <- ggplot(data =heat_data, aes(x = X, y = variable, fill = value)) +
  geom_tile() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  scale_fill_gradient(low = "white", high = "red")

line_data <- read.csv('totalSTRoutput(old).csv',sep = '\t')
line_plot <- ggplot(data = line_data, aes(x = X, y = totalSTR)) +
  geom_line(size =2) +
  theme_bw()+
  theme(panel.border = element_blank(), axis.line = element_line())+
  theme(panel.grid=element_blank())

combined_plot <- heatmap_plot / line_plot +
  plot_layout(ncol = 1, heights = c(4, 1))


print(combined_plot)



ggsave("combined_plot.pdf", plot = combined_plot, width = 10, height = 8)
