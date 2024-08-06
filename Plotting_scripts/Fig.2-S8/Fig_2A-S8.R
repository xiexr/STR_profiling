rm(list = ls())
setwd('E:/experiment/cento')
dat <- read.csv('STR_densitytocento.dat.txt',sep='\t')
library(dplyr)
library(zoo)
cPalette <- c("#5B9BD5", "#ED7D31", "#A5A5A5", "#FFC000", "#4472C4", "#70AD47")
names(cPalette) = c('p1', 'p2', 'p3', 'p4', 'p5', "p6")

alt_types = c('All','p1', 'p2', 'p3', 'p4', 'p5', "p6")
#Each period
plot.metaDist <- function(meta.dat, color, fill=T, norm=F, xlabel=NULL){
  #Clean meta dat
  meta.dat <- meta.dat[which(!is.na(meta.dat$mean)), ]
  if(norm==T){
    meta.dat$mean <- meta.dat$mean.norm
  }
  plot.vals <- rollapply(data=meta.dat$mean.norm, width=6, mean, partial=T)
  #Set parameters
  ymax <- 1.075*max(max(plot.vals, na.rm=T), 
                    2*mean(plot.vals, na.rm=T))
  lpos <- 1
  #Prep plot area
  par(bty="n", bg="white")
  plot(x=1.015*range(meta.dat$norm.pos), y=c(-0.05*ymax, ymax), type="n", 
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i")
  #Add background shading rectangle
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4], 
       border="gray99", col="gray99")
  box(which="plot", col="white", lwd=3)
  #Add contig positions & labels
  if(is.null(xlabel)){
    axis(lpos, at=mean(par("usr")[1:2]), 
         tick=F, labels="Meta-chromosome", line=-0.6)
  }else{
    axis(lpos, at=mean(par("usr")[1:2]), 
         tick=F, labels=xlabel, line=-0.6, cex.axis=0.8)
  }
  #Add y axis & gridlines
  y.at <- axTicks(2)
  x.at <- axTicks(1)
  print(x.at)
  y.at <- c(y.at, max(y.at)+(y.at[2]-y.at[1]))
  axis(2, at=y.at, labels=NA, tck=-0.05, col="gray40")

  axis(2, at=y.at, tick=F, las=2, cex.axis=1, line=-0.1, labels=round(y.at, 2))

  if(fill==T){
    polygon(x=c(meta.dat$norm.pos, rev(meta.dat$norm.pos)), 
            y=c(plot.vals, rep(0, times=length(plot.vals))), 
            border=NA, col=adjustcolor(color, alpha=0.7))
  }
  points(meta.dat$norm.pos, plot.vals, type="l", lwd=1.25, col=color)
  #Cleanup
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=0, 
       bty="n", border=NA, col="white")
  abline(h=0)
}

metaAverage <- function(dat, SVTYPE="All", n.bins=500){
  dat = as.data.frame(dat)
  p.bins <- seq(-1, 0, by=1/n.bins)
  q.bins <- seq(0, 1, by=1/n.bins)
  col.idx <- which(colnames(dat)==SVTYPE)
  p.means <- sapply(1:(length(p.bins)-1), function(i){
    mean(dat[which(dat$cdist.norm>=p.bins[i] & dat$cdist.norm<p.bins[i+1]), col.idx], na.rm=T)
  })

  q.means <- sapply(1:(length(q.bins)-1), function(i){
    mean(dat[which(dat$cdist.norm>q.bins[i] & dat$cdist.norm<=q.bins[i+1]), col.idx], na.rm=T)
  })
  means <- c(p.means, q.means)
  means.norm <- means/mean(means, na.rm=T)
  out.df <- data.frame("norm.pos"=c(p.bins[-length(p.bins)], q.bins[-1]), 
                       "mean"=means, "mean.norm"=means.norm)
  return(out.df)
}


meta.means <- lapply(alt_types, function(svtype){
  m <- metaAverage(dat=dat, SVTYPE=svtype, n.bins=500)
})
names(meta.means) <- alt_types



#{r p-4, fig.width=7, fig.height=7}
png("eachbpplot_image.png", width = 800, height = 600,res=300)
par(mfrow=c(3, 2), 
    mar=c(1.5, 3, 2, 1))
lapply(alt_types[-1], function(svtype){
  svt.col <- cPalette[which(names(cPalette) == svtype)]
  plot.metaDist(meta.dat=meta.means[[which(names(meta.means)==svtype)]], 
                color=svt.col, norm=T,xlabel=svtype)
})
dev.off()

pdf("eachbpplot_image.pdf", width = 8, height = 6)
par(mfrow=c(3, 2), 
    mar=c(1.5, 3, 2, 1))
lapply(alt_types[-1], function(svtype){
  svt.col <- cPalette[which(names(cPalette) == svtype)]
  plot.metaDist(meta.dat=meta.means[[which(names(meta.means)==svtype)]], 
                color=svt.col, norm=T,xlabel=svtype)
})
dev.off()







###totalSTRpic



png("Allplot_image.png", width = 800, height = 600,res=300)
par(mar=c(2, 3, 2, 1))
plot.metaDist(meta.dat=meta.means[[which(names(meta.means)=="All")]], 
              color='#3874A7', norm=T, 
              xlabel="Meta-chromosome (Mean of Autosomes)")
dev.off()

pdf("Allplot_image.pdf", width = 8, height = 6)
par(mar=c(2, 3, 2, 1))
plot.metaDist(meta.dat=meta.means[[which(names(meta.means)=="All")]], 
              color='#3874A7', norm=T, 
              xlabel="Meta-chromosome (Mean of Autosomes)")
dev.off()


#alt_types



