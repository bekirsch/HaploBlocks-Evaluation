# loading libraries
library(latex2exp)

args <- commandArgs(trailingOnly = TRUE) # 1=path
path=args[1]
files=as.numeric(args[2])
setwd(paste(path))

# Naming frequencies for pattern matching
frequencies <- c("0.0", "0.02", "0.05", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0")

# Create matrix with filenames as entries
s002 <- matrix(data = NA, nrow = files, ncol = 13)
for (j in 1:length(frequencies)) {
  s002[,j] <- c(list.files(pattern = paste("sC0.02_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = '')), rep(NA, length.out=files-length(list.files(pattern = paste("sC0.02_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = '')))))
}

# Create empty matrix for numeric entries
s002S = s002
class(s002S) = "numeric"

# Get maximum s-hat of blocks overlapping 5mb (the centre of the chromosome)
for (i in 1:13) {
  for (j in 1:files) {
    if (file.exists(paste(s002[j,i]))){
      filtered=read.csv(file = paste(s002[j,i]))
      neighbours=which(filtered$bp.end >= 5e6 & filtered$bp.start <= 5e6)
      if (length(neighbours) != 0) {
        s002S[j,i] <- max(filtered$sHat[neighbours])
      } else {
        s002S[j,i] <- 0
      }
    } else {
      s002S[j,i] = 0
    }
    
  }
}

# Set position for y labels
yrange=max(s002S)+0.011
perc=0.001/yrange
ylabPos=-2*perc*yrange*10

# Plotting
pdf('evaluation.pdf', pointsize = 10, title = 'evaluation', width = 10, height = 8)
op<-par(no.readonly=TRUE)
par(oma=c(2,2,2,2))
par(mai=c(0.5,0.5,0,0.01))

boxplot.matrix(s002S,
               las=2,
               col = 'grey',
               ylim=c(-0.001,max(s002S)+0.01),
               at=c(0, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
               boxwex=c(1, 1.5, 2,4,6,6,6,6,6,6,6,6,6),
               xaxt='n',
               boxlwd=0.01,
               medlwd=1.5,
               whisklwd=0.01,
               staplelwd=0.01,
               outlwd=0.01,
               outpch=1,
               cex=0.8)
abline(h=0.02, lty=2, col='red', lwd=0.5)
abline(v=1, lty=1, col='grey', lwd=0.5)

axis(side = 1,
     at = c(0, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
     labels = FALSE)
text(x = c(0, 2.5, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100), 
     y = ylabPos,
     labels = c('0%', '2%', '5%', '10%', '20%', '30%', '40%', '50%','60%', '70%', '80%', '90%', '100%'),
     xpd = NA,
     srt = 45,
     adj = 1,
     cex = 1)

mtext(text = TeX('Inferred selection coefficient $\\hat{s}'), side = 2, outer = TRUE, line = 0.5)
mtext('Frequency in population', side=1, line=3)

dev.off()
