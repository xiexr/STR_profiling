library(data.table)
library(CMplot)
library(viridis)
setwd('E:/experiment/STR_midu_heatmap')
map1 = fread("res.txt", select = c(2, 1, 3), col.names = c("STR", "Chromosome", "Position"),header = F)
table(map1$Chromosome)

summary(map1$Position)
hist(map1$Position, breaks = 50, main = "Distribution of Positions", xlab = "Position")

gradient_colors <- viridis(100, option = "plasma") 

print(gradient_colors)
CMplot(map1, plot.type = "d", bin.size = 1e6, col = gradient_colors,file = "tiff",dpi = 600, file.output = TRUE, verbose = TRUE)

