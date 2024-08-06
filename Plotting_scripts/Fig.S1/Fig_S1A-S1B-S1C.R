library(ggside)
library(ggplot2)
library(tidyverse)
library(ggrastr)
library(ggpmisc)


setwd('G:/experiment/callrateandcoverage')
rm(list=ls())
data<- data.table::fread('HipSTRtotalcallrateandcoverage.txt',sep = '\t')

ggplot(data,aes(x = coverage, y = callrate)) +
     geom_point_rast(size = 0.5, raster.dpi = 600, color = "#0096d6")+
     geom_xsideviolin(aes(y = coverage), orientation = "y",fill='#00b388')+
     scale_xsidey_discrete(labels = NULL) +
     geom_ysideviolin(aes(y = callrate), orientation = "x",fill='#00b388')+
    scale_ysidex_discrete(guide = guide_axis(angle = 45),
                          labels = NULL) +
    stat_correlation(label.y = "bottom",
                     label.x = "right") +
    labs(y = "Sample call rate (%)", x= "Autosome coverage")+
     ylim(0, 100)+
  geom_xsidepoint(data = data.frame(x=mean(data$coverage),y = 30), 
                  aes(x=x,y = y), color='#cc0000',size=5)+
  geom_ysidepoint(data = data.frame(x=32,y = mean(data$callrate)), 
                  aes(x=x,y = y), color='#cc0000',size=5)+
  theme_classic()

 
 ggsave('HipSTRcallrate.pdf')
####################################################
 ######################################
 #######################
 library(ggside)
 library(ggplot2)
 library(tidyverse)
 library(ggrastr)
 library(ggpmisc)

 setwd('G:/experiment/callrateandcoverage')
 rm(list=ls())
 data<- data.table::fread('GATKtotalcallrateandcoverage.txt',sep = '\t')
 ggplot(data,aes(x = coverage, y = callrate)) +
   geom_point_rast(size = 0.5, raster.dpi = 600, color = "#0096d6")+
   geom_xsideviolin(aes(y = coverage), orientation = "y",fill='#00b388')+
   scale_xsidey_discrete(labels = NULL) +
   geom_ysideviolin(aes(y = callrate), orientation = "x",fill='#00b388')+
   scale_ysidex_discrete(guide = guide_axis(angle = 45),
                         labels = NULL) +
   stat_correlation(label.y = "bottom",
                    label.x = "right") +
   labs(y = "Sample call rate (%)", x= "Autosome coverage")+
   ylim(0,100)+
 theme_classic()+
   geom_xsidepoint(data = data.frame(x=mean(data$coverage),y = 30), 
                   aes(x=x,y = y), color='#cc0000',size=5)+
   geom_ysidepoint(data = data.frame(x=32,y = mean(data$callrate)), 
                   aes(x=x,y = y), color='#cc0000',size=5)
 
 
 ggsave('GATKSTRcallrate.pdf')
###################################################
 ################################
 library(ggside)
 library(ggplot2)
 library(tidyverse)
 library(ggrastr)
 library(ggpmisc)

 
 setwd('G:/experiment/callrateandcoverage')
 rm(list=ls())
 data<- data.table::fread('LobSTRtotalcallrateandcoverage.txt',sep = '\t')
 mean_coverage <- mean(data$coverage)
 mean_callrate <- mean(data$callrate)

 ggplot(data,aes(x = coverage, y = callrate)) +
   geom_point_rast(size = 0.5, raster.dpi = 600, color = "#0096d6")+
   geom_xsideviolin(aes(y = coverage), orientation = "y",fill='#00b388')+
   scale_xsidey_discrete(labels = NULL) +
   geom_ysideviolin(aes(y = callrate), orientation = "x",fill='#00b388')+
   scale_ysidex_discrete(guide = guide_axis(angle = 45),
                         labels = NULL) +
   stat_correlation(label.y = "bottom",
                    label.x = "right") +
   labs(y = "Sample call rate (%)", x= "Autosome coverage")+
   theme_classic()+
   geom_xsidepoint(data = data.frame(x=mean(data$coverage),y = 30), 
                   aes(x=x,y = y), color='#cc0000',size=5)+
   geom_ysidepoint(data = data.frame(x=30,y = mean(data$callrate)), 
                   aes(x=x,y = y), color='#cc0000',size=5)
 
 
 ggsave('LobSTRcallrate.pdf')
 