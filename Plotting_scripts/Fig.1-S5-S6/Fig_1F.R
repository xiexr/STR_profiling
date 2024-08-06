library(ggplot2)
library(scales)
rm(list = ls())
setwd('D:/major_allele_diff_from_ref/')
df<- read.csv('majoralleledifffromref.txt',sep='\t')
count_data <- table(df$majordifffromref)
count_dataframe<-data.frame(count_data)
count_dataframe$Var1 <- as.factor(count_dataframe$Var1)
class(count_dataframe$Var1)
class(count_dataframe)
ggplot(count_dataframe,aes(x = Var1, y = Freq)) +
  geom_col(fill = "#A9A9A9") +
  scale_y_log10(
    expand = c(0, 0),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  #scale_x_continuous(breaks = seq(-15, 15, 5)) +
  labs(
    x = "Major allele repeat differences from reference",
    y = "# of pSTR"
  )+
  theme_classic()
ggsave('majoralleledifffromref.png',dpi=600)
ggsave('majoralleledifffromref.pdf')

write.csv(count_dataframe, file = "repeat_dif_from_ref.summary.csv", quote = F, row.names = F)
