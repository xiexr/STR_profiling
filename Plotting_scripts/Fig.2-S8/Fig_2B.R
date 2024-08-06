rm(list = ls())
setwd('E:/fold-enrichment')
dat <- read.csv('STR_densitytocento.dat.txt',sep='\t')
library(dplyr)
library(zoo)
cPalette <- c('black',"#5B9BD5", "#ED7D31", "#A5A5A5", "#FFC000", "#4472C4", "#70AD47")
names(cPalette) = c('All','p1', 'p2', 'p3', 'p4', 'p5', "p6")



alt_types = c('All','p1', 'p2', 'p3', 'p4', 'p5', "p6")
plot.metaByContext <- function(dat, meta.svtypes, colors){
  dat = as.data.frame(dat)
  #Get point estimates and CIs
  plot.vals <- calc.meanByContext(dat, meta.svtypes)
  #print(plot.vals)
  
  #Prep plot area
  ylims <- c(1.1*min(as.numeric(unlist(lapply(plot.vals, range, na.rm=T)))), 
             1.1*max(as.numeric(unlist(lapply(plot.vals, range, na.rm=T)))))
  par(bty="n", mar=c(0.25, 2.8, 1.8, 0.25))
  plot(x=c(0, length(plot.vals)-0.4), y=ylims, type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(h=1, col="gray50")
  #Add points per svtype
  sapply(1:length(plot.vals), function(i){
    segments(x0=i-c(1, 0.75, 0.5), 
             x1=i-c(1, 0.75, 0.5), 
             y0=as.numeric(plot.vals[[i]][1, ]), 
             y1=as.numeric(plot.vals[[i]][3, ]), 
             col=colors[i], lwd=2.5)
    print(plot.vals[[i]][1, ])
    points(x=i-c(1, 0.75, 0.5), 
           y=as.numeric(plot.vals[[i]][2, ]), 
           pch=15, col=colors[i], cex=3,bg=colors[i])
    text(x=i-c(1, 0.75, 0.5), 
         y=as.numeric(plot.vals[[i]][2, ]), 
         labels=c("T", "I", "C"), cex=1.3, font=2, col="white")
    #Add category label
    axis(3, at=i-c(1, 0.5), tck=0, labels=NA, line=0.1,lwd=2)
    axis(3, at=i-0.75, tick=F, line=-0.9, cex.axis=2, labels=meta.svtypes[i], col.axis=colors[i])
    #Add p-values
    par(xpd=T)
    sapply(1:3, function(k){
      if(plot.vals[[i]][4, k]<0.05/(3*length(plot.vals))){
        text(x=(i-c(1, 0.75, 0.5))[k], 
             y=plot.vals[[i]][2, k]-0.035, 
             pos=1, labels="*",cex=2)
      }
      if(plot.vals[[i]][5, k]<0.05/(3*length(plot.vals))){
        text(x=(i-c(1, 0.75, 0.5))[k], 
             y=plot.vals[[i]][2, k]+0.01, 
             pos=3, labels="*",cex=2)
      }
    })
    par(xpd=F)
  })
  #Clean up
  axis(2, at=axTicks(2), labels=NA, tck=-0.015)
  sapply(axTicks(2), function(y){
    axis(2, at=y, tick=F, las=2, cex.axis=1.5, line=-0.5, 
         labels=y)
  })
  mtext(2, line=1.75, text="SV Fold-Enrichment")
}
################################################################################
calc.meanByContext <- function(dat, meta.svtypes, ter.buf=0.05, cen.buf=0.05) {

  get.ci <- function(vals) {
    vals <- vals[which(!is.na(vals))] 
    k <- 1.96 * (sd(vals, na.rm=T) / sqrt(length(vals)))  
    p.less <- t.test(vals, mu=1, alternative="less")$p.value  
    p.greater <- t.test(vals, mu=1, alternative="greater")$p.value  
    return(c(c(mean(vals)-k, mean(vals), mean(vals)+k), p.less, p.greater)) 
  }
  
  res <- lapply(meta.svtypes, function(svtype) {
    
    
    vals <- as.numeric(dat[, which(colnames(dat) == svtype)])  
    vals <- vals / mean(vals, na.rm=T)  
    
    
    ter.idx <- which(dat$cdist.norm <= -1 + ter.buf | dat$cdist.norm >= 1 - ter.buf)  
    cen.idx <- which(dat$cdist.norm >= -cen.buf & dat$cdist.norm <= cen.buf)
    int.idx <- which(!(1:nrow(dat) %in% c(ter.idx, cen.idx)))  
    
    
    ter.stats <- get.ci(vals[ter.idx])
    int.stats <- get.ci(vals[int.idx])  
    cen.stats <- get.ci(vals[cen.idx])  
    
    return(data.frame("ter" = ter.stats, "int" = int.stats, "cen" = cen.stats))  # 创建包含计算出的统计数据的数据框
    })
  
  names(res) <- meta.svtypes  
  return(res)  
}


pdf("enrichmentwithpval.pdf", width = 8, height = 6)
plot.metaByContext(dat, alt_types, 
                   colors=cPalette)
dev.off()

