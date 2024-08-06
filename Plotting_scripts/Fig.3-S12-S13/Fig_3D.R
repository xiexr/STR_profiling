rm(list = ls())
testsuffix<-function(fn,suff){
  parts<-strsplit(fn,".",fixed = TRUE)
  nparts<-length(parts[[1]])
  return(parts[[1]][nparts]==suff)
}
library(ggplot2)
library(hash)
library(EnvStats)
library(ggforce)
setwd("D:/expression_example")
file = list.files()
for (i in file){if(testsuffix(i,'txt')){
  print(i)
  parts1<-strsplit(i,".",fixed = TRUE)
  parts2<-strsplit(parts1[[1]][1],"=",fixed = TRUE)
  if (parts2[[1]][length(parts2[[1]])]=='None'){
    titlename<- parts2[[1]][2]}
  else{titlename<- parts2[[1]][3]}
  phenotype = read.table(i,header=T,sep="\t")
  a<-phenotype$Repeattimes[!duplicated(phenotype$Repeattimes)]
  ###
  h=hash()
  b <- c()
  j=1
  for(k in a){
    parts3<-strsplit(k," ",fixed = TRUE)
    key=parts3[[1]][1]
    b <- append(b,key)
    .set(h,keys=key,values=k)
    j=j+1}
  b<- as.numeric(b)
  b <- sort(b)
  d<- c()
  for (m in b){
    m<- as.character(m)
    ankey<- h[[m]]
    d<- append(d,ankey)
  }
  parts4 <- strsplit(parts2[[1]][1],"-",fixed = TRUE)
  titlename <- paste(titlename,' (Unit=',parts4[[1]][2],')')
  phenotype$Repeattimes<-as.character(phenotype$Repeattimes)
  phenotype$Repeattimes <- factor(phenotype$Repeattimes,levels=d,order=TRUE)
  vMBN <- ggplot(phenotype, aes(x=Repeattimes, y=Expression, fill=Repeattimes))+
    geom_violin(fill='#e6e6e6') + 
    geom_boxplot(width=0.1,color="black",fill="#b3dcff",outlier.colour = NA)+
    #geom_jitter(width = 0.2, color='#58595b', size =3)+
    stat_summary(fun=median, geom="line", aes(group=1),color='#f79400',width=100,size=2)  + 
    geom_sina(size=2.3,color='#56b881')+
    #(fun=median, geom="point")+
    #stat_n_text() + 
    ylim(0,2)+
    labs(title = titlename,x='Repeat times',y="Gene expression")+
    guides(fill = FALSE)
  mytheme2 <- theme(plot.title = element_text(family='Arial',face = "bold",
                                              size = "30", color = "black",hjust = 0.5,vjust = 0.5),
                    axis.title = element_text(family='Arial',face = "bold",
                                              size = "25",color = "black"),
                    axis.text.x = element_text(family='Arial',face = "bold",
                                               size = 20, angle = 0, hjust = 0.5, vjust = 1),
                    axis.text.y = element_text(family='Arial',face = "bold",size = 20),
                    panel.background = element_rect(fill = "white",color = "black"),
                    legend.position = "none",
                    panel.grid.major.y = element_blank(),
                    panel.grid.minor.y = element_blank(),
                    panel.grid.minor.x = element_blank())
  filename<-paste(parts1[[1]][1],".png",sep="")
  vMBNfinal = vMBN +  mytheme2
  ggsave(file=filename,dpi=300,plot=vMBNfinal, width = 12, height = 6)
  
}
}
