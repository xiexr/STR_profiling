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
    return(data.frame("ter" = ter.stats, "int" = int.stats, "cen" = cen.stats))  
  })
  
  names(res) <- meta.svtypes  
  return(res) 
}

setwd('E:/experiment/STR')
dat <- read.csv('STR_densitytocento.dat.txt',sep='\t')
dat = as.data.frame(dat)
plot.vals <- calc.meanByContext(dat, meta.svtypes)



all_results <- list()

plot.vals <- calc.meanByContext(dat, alt_types)

for (i in 1:length(plot.vals)) {
  svtype <- alt_types[i]
  svtype_result <- as.data.frame(plot.vals[[i]])
  svtype_result$svtype <- svtype
  all_results[[i]] <- svtype_result
}

final_result <- do.call(rbind, all_results)

write.csv(final_result, file = "output_results.csv", row.names = FALSE)
