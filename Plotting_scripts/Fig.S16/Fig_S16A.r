##################################################################
setwd('E:/experiment/STR_pairs_LD')
data <- read.table("OutputPrefix.bin")


pdf("OutputPrefix_STR_STR.pdf")
plot(data[,1]/1000, data[,2], type="l", col="#606060", main="LD decay", 
     xlab="Distance(Kb)", xlim=c(0,250), ylab=expression(r^{2}), bty="n", 
     lwd=4, axes=FALSE, ylim=c(0, max(data[,2])), xaxs="i")


axis(side = 1, at = pretty(data$V1/1000), labels = TRUE, pos=0)  
axis(side = 2, at = seq(0, 1, 0.1), labels = TRUE, pos=0)  

abline(h = 0.1, col = "black", lty = 2)
dev.off()


png("OutputPrefix_STR_STR.png")
plot(data[,1]/1000, data[,2], type="l", col="#606060", main="LD decay", 
     xlab="Distance(Kb)", xlim=c(0,250), ylab=expression(r^{2}), bty="n", 
     lwd=4, axes=FALSE, ylim=c(0, max(data[,2])), xaxs="i")


axis(side = 1, at = pretty(data$V1/1000), labels = TRUE, pos=0)  
axis(side = 2, at = seq(0, 1, 0.1), labels = TRUE, pos=0)  

abline(h = 0.1, col = "black", lty = 2)
dev.off()
