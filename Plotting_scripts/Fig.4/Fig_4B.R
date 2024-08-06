library(tidyverse)
library(networkD3)
library(dplyr)
library(RColorBrewer)
library(webshot)


# Print all hexadecimal color values
print(colorArray)
setwd('E:/experiment/sangkitu')

data <- read.table("GWAS_pop_sangkitu.txt", sep = '\t', header = TRUE)
nodes <- data.frame(name = c(as.character(data$pop1_pop2), as.character(data$GWAS_trait)) %>% unique()) # ÖÆ×÷ nodes
data$IDsource <- match(data$pop1_pop2, nodes$name) - 1
data$IDtarget <- match(data$GWAS_trait, nodes$name) - 1
nodes$group <- ifelse(nodes$name %in% as.character(data$pop1_pop2), nodes$name, "right")

data$groupdata <- paste(data$pop1_pop2, data$GWAS_trait)

#
ColourScal <- 'd3.scaleOrdinal().range(["#66CCFF", "#aec7e8", "#FFCC00", "#ffbb78", "#CCCCFF", "#98df8a", "#A9A9A9"])'

a <- sankeyNetwork(
  Links = data, Nodes = nodes,
  Source = "IDsource", Target = "IDtarget",LinkGroup = 'groupdata', 
  Value = "value", NodeID = "name",NodeGroup = "group",colourScale = ColourScal,
  nodeWidth = 40, fontSize = 0, nodePadding = 20
)
print(a)
saveNetwork(a,"sankey.html")
webshot("sankey.html" , "sankey.png")
webshot("sankey.html" , "sankey.pdf")
webshot::webshot(a, file = "sankey_plot.pdf", cliprect = "viewport")
webshot::webshot(a, file = "sankey_plot.png", cliprect = "viewport")


###################################################

