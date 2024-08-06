setwd("E:/gatkhiplobveen")
rm(list = ls())
library (VennDiagram)

library(extrafont)
hipSTR <- read.csv('hipstr.txt', header = FALSE, sep = ',')
lobSTR <- read.csv('lobstr.txt', header = FALSE, sep = ',')
gatk <- read.csv('gatk.txt', header = FALSE, sep = ',')
venn_list <- list(lobstr = lobSTR$V1, gatk = gatk$V1,hipstr=hipSTR$V1)
venn.diagram(venn_list,filename = 'STRVeenpiasdc.png',category.names=c("LobSTR","GATK","HipSTR"), cat.pos=c(0,0,180),main.pos=c(0.5,1.0,1.5),cat.dist=c(0.035,0.035,0.035),cat.default.pos='outer',height = 6000, width = 6000, resolution =
               600,imagetype = 'png',main='',main.fontfamily = "Times New Roman",main.cex = 5.5,
             fill = c(colors()[468], colors()[616],colors()[25]), alpha = 0.85,
             cat.col = c('black', 'black','black'), cat.cex = 0.000001, cat.fontfamily = 'Times New Roman',
             lwd=3,col = 'white', cex = 0, fontfamily = "Times New Roman",display.mode = "none")
inter <- get.venn.partitions(venn_list)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
write.table(inter[-c(3, 4)], 'venn_inter.txt', row.names = FALSE, sep = '/t', quote = FALSE)
a<-union(JapArovenn$lobstr,JapArovenn$hipstr)
b<- union(a,JapArovenn$gatk)