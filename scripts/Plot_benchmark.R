args <- commandArgs(trailingOnly = TRUE) # 1=path

path=args[1]

setwd(paste(path))

blIND = as.data.frame(t(read.table(file = 'haploblocksIND.bench')))
blSNP = as.data.frame(t(read.table(file = 'haploblocksSNP.bench')))

biIND = as.data.frame(t(read.table(file = 'hapbinIND.bench')))
biSNP = as.data.frame(t(read.table(file = 'hapbinSNP.bench')))

options(scipen = 8)

MEMmin = min(biIND$V3, blIND$V3, biSNP$V3, blSNP$V3)/1e3
MEMmax = max(biIND$V3, blIND$V3, biSNP$V3, blSNP$V3)/1e3

WALLmin = min(biIND$V2, blIND$V2, biSNP$V2, blSNP$V2)
WALLmax = max(biIND$V2, blIND$V2, biSNP$V2, blSNP$V2)

### INDV ###
pdf('BenchINDV.pdf', pointsize = 10, width = 8, height = 4, )
#op<-par(no.readonly=TRUE)
#par(oma=c(2,2,1,1))
par(mfrow=c(1,2), mai=c(0.55,0.7,0.15,0.1))
plot(blIND$V1, blIND$V2, type = 'n', pch=0, 
     xlim = c(min(biIND$V1)/2, max(biIND$V1)/2),
     ylim = c(WALLmin, WALLmax), 
     xlab = '', 
     ylab = 'Runtime (s)',
     log = 'xy',
     xaxt = 'n')
for (i in seq(0.1,1, by=0.1)) {
        abline(h=i, col='grey')  
}
for (i in seq(2,10, by=1)) {
        abline(h=i, col='grey')  
}
for (i in seq(20,100, by=10)) {
        abline(h=i, col='grey')  
}
for (i in seq(200,1000, by=100)) {
        abline(h=i, col='grey')  
}
for (i in seq(2000,10000, by=1000)) {
        abline(h=i, col='grey')  
}
lines(biIND$V1/2, biIND$V2, type = 'o', cex=1.5)
lines(blIND$V1, blIND$V2, type = 'o', pch=0, cex=1.5)
legend(25, max(biIND$V2)*1, 
       legend = c("hapbin", "haploblocks"), 
       lty = c(1, 1),
       pch = c(1, 0), bg = 'white')
axis(side = 1, at = c(25, 100, 500, 1000, 2500), labels = FALSE)
text(cex=1, x=c(25, 100, 500, 1000, 2500), 
     y=rep(WALLmin-WALLmin/1.5, 5), 
     c("25", "100", "500", "1000", "2500"), 
     las = 3, xpd=TRUE)
title(xlab="Number of individuals", line=2)
        

plot(blIND$V1, blIND$V3/1e3, type = 'n', pch=4, 
     lty=1, 
     ylim = c(MEMmin, MEMmax),
     xlab = '',
     ylab = 'Peak memory usage (mb)',
     log = 'xy',
     xaxt = 'n')
for (i in seq(10,100, by=10)) {
        abline(h=i, col='grey')  
}
for (i in seq(200,1000, by=100)) {
        abline(h=i, col='grey')  
}

lines(biIND$V1/2, biIND$V3/1e3, type = 'o', lty=1, cex=1.5)
lines(blIND$V1, blIND$V3/1e3, type = 'o', pch=0, lty=1, cex=1.5)
legend(25, MEMmax*1, 
       legend = c("hapbin", "haploblocks"), 
       lty = c(1, 1),
       pch = c(1, 0), bg = 'white')
axis(side = 1, at = c(25, 100, 500, 1000, 2500), labels = FALSE)
text(cex=1, x=c(25, 100, 500, 1000, 2500), 
     y=rep((MEMmin-log(MEMmin))/1.085, 5), 
     c("25", "100", "500", "1000", "2500"), 
     las = 3, xpd=TRUE)
title(xlab="Number of individuals", line=2)
dev.off()

### SNPS ###
pdf('BenchSNPS.pdf', pointsize = 10, width = 8, height = 4, )
#op<-par(no.readonly=TRUE)
#par(oma=c(2,2,1,1))
par(mfrow=c(1,2), mai=c(0.55,0.7,0.15,0.1))
plot(blSNP$V1, blSNP$V2, type = 'n', 
     pch=4, 
     xlim = c(min(biSNP$V1), max(biSNP$V1)),
     ylim = c(WALLmin, WALLmax), 
     xlab = '', 
     ylab = 'Runtime (s)', 
     log = 'xy',
     xaxt = 'n')
for (i in seq(0.1,1, by=0.1)) {
        abline(h=i, col='grey')  
}
for (i in seq(2,10, by=1)) {
        abline(h=i, col='grey')  
}
for (i in seq(20,100, by=10)) {
        abline(h=i, col='grey')  
}
for (i in seq(200,1000, by=100)) {
        abline(h=i, col='grey')  
}
for (i in seq(2000,10000, by=1000)) {
        abline(h=i, col='grey')  
}
lines(biSNP$V1, biSNP$V2, type = 'o', cex=1.5)
lines(blSNP$V1, blSNP$V2, type = 'o', pch=0, cex=1.5)
axis(side = 1, at = c(100, 1000, 10000, 100000), labels = FALSE)
text(cex=1, x=c(100, 1000, 10000, 100000), 
     y=rep(WALLmin-WALLmin/1.5, 5), 
     c("100", "1000", "10000", "100000"), 
     las = 3, xpd=TRUE)
legend(100, WALLmax*1, 
       legend = c("hapbin", "haploblocks"), 
       lty = c(1, 1),
       pch = c(1, 0), bg = 'white')
title(xlab="Number of SNPs", line=2)

plot(blSNP$V1, blSNP$V3/1e3, type = 'n', pch=0, 
     lty=1, 
     ylim = c(MEMmin, MEMmax),
     xlab = '',
     ylab = 'Peak memory usage (mb)',
     log = 'xy',
     xaxt = 'n')
for (i in seq(10,100, by=10)) {
        abline(h=i, col='grey')  
}
for (i in seq(200,1000, by=100)) {
        abline(h=i, col='grey')  
}
lines(biSNP$V1, biSNP$V3/1e3, type = 'o', lty=1, cex=1.5)
lines(blSNP$V1, blSNP$V3/1e3, type = 'o', pch=0, lty=1, cex=1.5)
axis(side = 1, at = c(100, 1000, 10000, 100000), labels = FALSE)
text(cex=1, x=c(100, 1000, 10000, 100000), 
     y=rep((MEMmin-log(MEMmin))/1.084, 5), 
     c("100", "1000", "10000", "100000"), 
     las = 3, xpd=TRUE)
#mtext("Chromosome 22", side = 3, line = -3, outer = TRUE)
legend(100, MEMmax*1, 
       legend = c("hapbin", "haploblocks"), 
       lty = c(1, 1),
       pch = c(1, 0), bg = 'white')
title(xlab="Number of SNPs", line=2)

dev.off()

