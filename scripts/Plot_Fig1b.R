# loading libraries
library(latex2exp)

args <- commandArgs(trailingOnly = TRUE)  # 1=path
path=args[1]
files=as.numeric(args[2])                 # 2=Number of files
setwd(paste(path))

cbcol  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#0072B2", "#D55E00", "#CC79A7")

###################################################################

gravel_EU=list.files(path = '../output_gravel_CEU/', pattern = paste("_filtered.sHat.csv", sep = ''), recursive = TRUE)
gravel_EU_count=str_replace(gravel_EU, '.vcf_filtered.sHat.csv', '.vcf.gz.count')

files = length(gravel_EU)


gravel_EUS = gravel_EU
gravel_EU_countS = gravel_EU_count
class(gravel_EUS) = "numeric"
class(gravel_EU_countS) = "numeric"

for (i in 1:length(gravel_EU)) {
  gravel_EU_countS[i]=read.table(paste('../output_gravel_CEU/', gravel_EU_count[i], sep = ''))$V1/2000
}

for (i in 1:length(gravel_EU)) {
  
  if (file.exists(paste('../output_gravel_CEU/', gravel_EU[i], sep = ''))){
    filtered=read.csv(file = paste('../output_gravel_CEU/', gravel_EU[i], sep = ''))
    neighbours=which(filtered$bp.start <= 5e6 & filtered$bp.end >= 5e6)#which(filtered$bp.end >=4000000 & filtered$bp.end <=6000000 & filtered$bp.start >=4000000 & filtered$bp.start <=6000000)
    if (length(neighbours) != 0) {
      gravel_EUS[i] <- max(filtered$sHat[neighbours])
    } else {
      gravel_EUS[i] <- 0
    }
  } else {
    gravel_EUS[i] = 0
  }
}

maxg=max(gravel_EU_countS)
nbins=10
bins <- as.data.frame(matrix(data = NA, nrow = 200, ncol = 11))
names(bins)=c('0', '(0-0.1]', '(0.1-0.2]', '(0.2-0.3]', '(0.3-0.4]', '(0.4-0.5]', '(0.5-0.6]', '(0.6-0.7]', '(0.7-0.8]', '(0.8-0.9]', '(0.9-1]')

bins$`0`[1:length(which(gravel_EU_countS <= 0))]=gravel_EUS[which(gravel_EU_countS <= 0)]
bins$`(0-0.1]`[1:length(which(gravel_EU_countS > 0 & gravel_EU_countS <= 0.1))]=gravel_EUS[which(gravel_EU_countS > 0 & gravel_EU_countS <= 0.1)]
bins$`(0.1-0.2]`[1:length(which(gravel_EU_countS > 0.1 & gravel_EU_countS <= 0.2))]=gravel_EUS[which(gravel_EU_countS > 0.1 & gravel_EU_countS <= 0.2)]
bins$`(0.2-0.3]`[1:length(which(gravel_EU_countS > 0.2 & gravel_EU_countS <= 0.3))]=gravel_EUS[which(gravel_EU_countS > 0.2 & gravel_EU_countS <= 0.3)]
bins$`(0.3-0.4]`[1:length(which(gravel_EU_countS > 0.3 & gravel_EU_countS <= 0.4))]=gravel_EUS[which(gravel_EU_countS > 0.3 & gravel_EU_countS <= 0.4)]
bins$`(0.4-0.5]`[1:length(which(gravel_EU_countS > 0.4 & gravel_EU_countS <= 0.5))]=gravel_EUS[which(gravel_EU_countS > 0.4 & gravel_EU_countS <= 0.5)]
bins$`(0.5-0.6]`[1:length(which(gravel_EU_countS > 0.5 & gravel_EU_countS <= 0.6))]=gravel_EUS[which(gravel_EU_countS > 0.5 & gravel_EU_countS <= 0.6)]
bins$`(0.6-0.7]`[1:length(which(gravel_EU_countS > 0.6 & gravel_EU_countS <= 0.7))]=gravel_EUS[which(gravel_EU_countS > 0.6 & gravel_EU_countS <= 0.7)]
bins$`(0.7-0.8]`[1:length(which(gravel_EU_countS > 0.7 & gravel_EU_countS <= 0.8))]=gravel_EUS[which(gravel_EU_countS > 0.7 & gravel_EU_countS <= 0.8)]
bins$`(0.8-0.9]`[1:length(which(gravel_EU_countS > 0.8 & gravel_EU_countS <= 0.9))]=gravel_EUS[which(gravel_EU_countS > 0.8 & gravel_EU_countS <= 0.9)]
bins$`(0.9-1]`[1:length(which(gravel_EU_countS > 0.9 & gravel_EU_countS <= 1))]=gravel_EUS[which(gravel_EU_countS > 0.9 & gravel_EU_countS <= 1)]

cbcol= c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888
         ")



options(scipen=10)

inch <- function(mm) {mm*0.0393701}

##########################################################################################################################################################################################################
cex=1
pdf('Fig_1b.pdf', 
    #paper = 'a4', 
    pointsize = 7, 
    height = inch(140), 
    width = inch(180))

boxplot(bins, xlim=c(0,10), yaxt='n', ylab='', xaxt='n', xlab='', boxwex=c(0.3,rep(0.6, length.out=10)), at=(c(0.3, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)-0.3))
#mui=boxplot(bins$`0`, add = TRUE, at = 0, boxwex=1.01)
abline(h=0.03, lty=2, col=cbcol[6], lwd=0.8)

axis(side = 1, 
     at=c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),  
     labels = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"))

axis(side = 2, at=seq(0,0.04, by=0.01))

mtext(TeX('Frequency of selected variant in population'), side = 1, outer = TRUE, line = -3.5, adj=0.5)
mtext('b', side = 3, at=-1.5, font = 2)
#abline(v=c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10))

#legend(5.5, 0.025, xjust = 0.5, yjust = 0.5, legend = '!Model Misspecification Not Final!    ', bg = 'red')

dev.off()


