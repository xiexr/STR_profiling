######################################  
rm(list = ls())
setwd("E:/Heter")
library(ggplot2)
library(extrafont)
phenotype = read.table("STRhetall_first_25.txt",header=T,sep="\t")
library(RColorBrewer)
colourCount = length(unique(phenotype$Species))
getPalette = colorRampPalette(brewer.pal(9, "Set3"))
vMBN <- ggplot(phenotype, aes(x=factor(Species,level=c("All","All_Indica","All_Japonica","Aus","Aromatic","Intermediate")), y=Het, fill=Species))+
  scale_y_continuous(labels=c(0,0.25,0.5,0.75,1), limits=c(0,1))+
  geom_boxplot(width = 0.6,color="black", outlier.shape = NA,lwd = 1)+
  scale_fill_manual(values = getPalette(colourCount))+ 
  labs(title = "Heterozygosity of STR in each population",x=NULL,y="Heterozygosity")+
  guides(fill = 'none')+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title = element_text(size = 10, color = "black", face = "bold"),
        axis.text.x = element_text(size = 5, angle = 45, hjust = 1, vjust = 1, face = "bold"))
mytheme2 <- theme(plot.title = element_text(family='Arial',face = "bold",
                                            size = "30", color = "black",hjust = 0.5,vjust = 0.5),
                  axis.title = element_text(family='Arial',face = "bold",
                                            size = "25",color = "black"),
                  axis.text.x = element_text(family='Arial',face = "bold",
                                             size = 20, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(family='Arial',face = "bold",size = 20),
                  panel.background = element_rect(fill = "white",color = "black"),
                  legend.position = "none",
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  panel.grid.minor.x = element_blank())
vMBNfinal = vMBN +  mytheme2
ggsave(file="STRhet.png",dpi=300,plot=vMBNfinal, width = 12, height = 10)
#######################################
