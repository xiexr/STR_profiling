setwd('E:/experiment/SNP_STR_Pair_LD')
dat=read.csv('output_for_SNP_STR_LD_paint.txt',sep = '\t',header=T)
library(zoo)
plotvals <- rollapply(data=dat$averLD, width=10, mean, partial=T)
pdf("OutputPrefix_SNP_STR.pdf")
plot(dat[,1]/1000,plotvals,type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,250),ylab=expression(r^{2}),bty="n",lwd=4,axes=FALSE,ylim=c(0, 0.2),xaxs="i")

axis(side = 1, at = pretty(dat$dist/1000), labels = TRUE,pos=0)  
axis(side = 2, at = seq(0, 0.3, 0.05), labels = TRUE,pos=-4)
dev.off()

png("OutputPrefix_SNP_STR.png")
plot(dat[,1]/1000,plotvals,type="l",col="blue",main="LD decay",xlab="Distance(Kb)",xlim=c(0,250),ylab=expression(r^{2}),bty="n",lwd=4,axes=FALSE,ylim=c(0, 0.2),xaxs="i")


axis(side = 1, at = pretty(dat$dist/1000), labels = TRUE,pos=0)  
axis(side = 2, at = seq(0, 0.35, 0.05), labels = TRUE,pos=-4)
dev.off()

