rm(list = ls())
setwd("D:repeat_unit_time_change/")
library(ggplot2)
Repeatunit = read.table("Repeat_unit_major_allele_repeat_res.txt",header=T,sep="\t")
library(RColorBrewer)
pMBN <- ggplot(Repeatunit, aes(x = reorder(repeat_unit, -majorallele_repeattime), y=majorallele_repeattime, fill=repeat_unit)) +
  geom_boxplot(fill = "#099963")+
  labs(title = "",x="Major allele unit",y="Repeat time")
print(pMBN)
mytheme2 <- theme(plot.title = element_text(family='Arial',face = "bold",
                                            size = "30", color = "brown",hjust = 0.5,vjust = 0.5),
                  axis.title = element_text(family='Arial',face = "bold",
                                            size = "25",color = "blue"),
                  axis.text.x = element_text(family='Arial',face = "bold",
                                             size = 25, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(family='Arial',face = "bold",size = 25),
                  axis.text = element_text(colour = "black"),
                  axis.ticks =element_line(color="black",size=2),
                  panel.background = element_rect(fill = "white",color = "black"),
                  panel.border = element_rect(color = "black", size = 2, fill = NA),
                  legend.position = "none",
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank())


pMBNfinal = pMBN + mytheme2

ggsave(file="out_box_plot.png",dpi=300,plot=pMBNfinal, width = 12, height = 8)
