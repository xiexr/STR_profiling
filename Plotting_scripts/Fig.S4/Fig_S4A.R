library(ggplot2)
library(dplyr)
library(extrafont)


setwd('E:/experiment/allelenum')
df <- read.csv("E:/res_of_each_allele.txt", sep="\t")
font_import(pattern = "arial", prompt = FALSE)

data <- data.frame(
  period = factor(df$period),
  allelenum = as.numeric(as.character(df$allelenum))  
)


grouped_data <- data %>%
  group_by(period) %>%
  summarise(allelenum = list(allelenum))
ggplot(data, aes(x = period, y = allelenum)) +
  geom_violin() +
  xlab("") +
  ylab("") +
  ggtitle("")+
  theme_classic()+
  theme(
    
    axis.text.x = element_text(size = 30, color = "black"),  
    axis.text.y = element_text(size = 30, color = "black", vjust = 0.5),  
    legend.text = element_text(size = 25, family = "Arial"),    
    legend.title = element_text(size = 25, family = "Arial"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks = element_line(linewidth = 1.2, color = "black"),      
    axis.line = element_line(linewidth = 1.2, color = "black")        
  )


ggsave('violin.pdf',device = "pdf",width = 8, height = 10)
