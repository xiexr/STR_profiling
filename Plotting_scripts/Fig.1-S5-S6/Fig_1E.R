library(ggplot2)
library(scales)
rm(list = ls())
setwd('E:/majorallelefreq')
df<- read.csv('majorallelefreqof16bp.txt',sep='\t')
p.numallele <- 
  ggplot(df,aes(x = majorallelefreq)) +
  geom_histogram(binwidth = .01, fill = "#3765A9") +
  scale_y_log10(
    expand = c(0, 0),
    breaks = trans_breaks("log10", function(x) 10^x),
    labels = trans_format("log10", math_format(10^.x))
  ) +
  labs(x = "major_allele_freq", y = "# of pSTR")+
  theme_classic()
p.numallele
ggsave('major_allele_freq.png',dpi=600)
ggsave('major_allele_freq.pdf')