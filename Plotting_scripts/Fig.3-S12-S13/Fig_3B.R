library(ggplot2)
library(scales)
library(patchwork)
library(cowplot)
library(ggpubr)
library(ggcorrplot)
library(ggrepel)
library(ggrastr)
library(EnvStats)
library(tidyverse)
library(magrittr)
library(dplyr)
rm(list = ls())
dfmock <- read.csv('MockPvalue_FPKM_1.txt')
dfreal <- read.csv('RealPvalue_FPKM_1.txt')
# Observed data
dfreal$p
y = sort(-log10(dfreal$p))

x = sort(-log10(runif(length(y))))

# Control data
yctrl = sort(-log10(dfmock$p))
xctrl = sort(-log10(runif(length(yctrl))))
d4p1 = data.frame(y2=dfreal$X,
  y1 = y,
                  x1 = x)
d4p2 = data.frame(y2=dfmock$X,yctrl=yctrl,xctrl=xctrl)

df2 <- full_join(d4p1,d4p2,by='y2')


p = ggplot(df2) +
  geom_point_rast(aes(x=x1, y=y1),color="red", size = .3, raster.dpi = 600, dev = "ragg") +
 geom_point_rast(aes(x=xctrl,y=yctrl),color = "gray", size = .3, raster.dpi = 600, dev = "ragg") +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = "Expected P value (-log10)",
       y = "Observed P value (-log10)")+
theme_bw()+theme(panel.grid=element_blank(),panel.border=element_blank())+
  theme(axis.line.x=element_line(linetype=1,color="black",linewidth=.5))+
  theme(axis.line.y=element_line(linetype=1,color="black",linewidth=.5))+
  theme(axis.ticks.x=element_line(color="black",linewidth=.5,lineend = 1))+
  theme(axis.ticks.y=element_line(color="black",linewidth=.5,lineend = 10))


ggsave(plot=p, file="Fig2.qqplot.png",dpi=300,
       height = 2.5, width = 3.5)

ggsave("Fig.qqplot.pdf",device = cairo_pdf,width =3.15, height =2.36)

