rm(list = ls())
setwd("E:/Rst")
library(ggplot2)
Repeatunit = read.table("Rst_cal_last_res_for_box.txt",header=T,sep="\t")
library(RColorBrewer)
custom_colors <- c("#76D1FF","#D1D1FF", "#FFD11A","#FFC286",  "#A2E296","#B6CDEA"  )
custom_order <- c('AllIndica_AllJaponica', 'AllJaponica_Aus', 'AllIndica_Aus','AllJaponica_Aromatic','Aus_Aromatic', 'AllIndica_Aromatic')
Repeatunit$pop1_pop2 <- factor(Repeatunit$pop1_pop2, levels = custom_order)
pMBN <- ggplot(Repeatunit, aes(x = pop1_pop2, y=Rst, fill=pop1_pop2)) +
  geom_boxplot(fill = custom_colors,width = 0.5, color = "black", outlier.shape = NA,lwd = 1)+
  labs(title = "",x="pop1_pop2",y="Rst")
print(pMBN)
mytheme2 <- theme(plot.title = element_text(family='Arial',face = "bold",
                                            size = "30", color = "brown",hjust = 0.5,vjust = 0.5),
                  axis.title = element_text(family='Arial',face = "bold",
                                            size = "25",color = "blue"),
                  axis.text.x = element_text(family='Arial',face = "bold",
                                             size = 20, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(family='Arial',face = "bold",size = 15),
                  panel.background = element_rect(fill = "white",color = "black"),
                  legend.position = "none",
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank())


pMBNfinal = pMBN + mytheme2

ggsave(file="out_box_plot.png",dpi=300,plot=pMBNfinal, width = 12, height = 10)
######################################  小提琴加箱线图

