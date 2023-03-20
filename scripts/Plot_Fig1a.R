# loading libraries
library(latex2exp)

args = commandArgs(trailingOnly = TRUE) # 1=path
path = args[1]
files = as.numeric(args[2])
setwd(paste(path))

cbcol  = c("#000000", "#E69F00", "#56B4E9", "#009E73",
"#0072B2", "#D55E00", "#CC79A7")

# Naming frequencies for pattern matching
frequencies = c("0.0", "0.02", "0.05", "0.1", "0.2", "0.3", "0.4",
"0.5", "0.6", "0.7", "0.8", "0.9", "1.0")

####### s = 0.0075
s00075 = matrix(data = NA, nrow = files, ncol = 13)
for (j in 1:length(frequencies)) {
  s00075[,j] = c(list.files(pattern = paste("sC0.0075_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")), rep(NA, length.out=files-length(list.files(pattern = paste("sC0.0075_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")))))
}

s00075S = s00075
class(s00075S) = "numeric"


# Get maximum s-hat of blocks overlapping 5mb (the centre of the chromosome)
for (i in 1:13) {
  for (j in 1:files) {
    if (file.exists(paste(s00075[j, i]))) {
      filtered = read.csv(file = paste(s00075[j, i]))
      neighbours = which(filtered$bp.end >= 5e6 & filtered$bp.start <= 5e6)
      if (length(neighbours) != 0) {
        s00075S[j, i] = max(filtered$sHat[neighbours])
      } else {
        s00075S[j, i] = 0
      }
    } else {
      s00075S[j, i] = 0
    }
    
  }
}

####### s = 0.01
s001 = matrix(data = NA, nrow = files, ncol = 13)
for (j in 1:length(frequencies)) {
  s001[, j] = c(list.files(pattern = paste("sC0.01_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")), rep(NA, length.out=files-length(list.files(pattern = paste("sC0.01_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")))))
}

s001S = s001
class(s001S) = "numeric"


# Get maximum s-hat of blocks overlapping 5mb (the centre of the chromosome)
for (i in 1:13) {
  for (j in 1:files) {
    if (file.exists(paste(s001[j, i]))) {
      filtered = read.csv(file = paste(s001[j, i]))
      neighbours = which(filtered$bp.end >= 5e6 & filtered$bp.start <= 5e6)
      if (length(neighbours) != 0) {
        s001S[j, i] = max(filtered$sHat[neighbours])
      } else {
        s001S[j, i] = 0
      }
    } else {
      s001S[j, i] = 0
    }
    
  }
}

####### s = 0.02
s002 = matrix(data = NA, nrow = files, ncol = 13)
for (j in 1:length(frequencies)) {
  s002[, j] = c(list.files(pattern = paste("sC0.02_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")), rep(NA, length.out=files-length(list.files(pattern = paste("sC0.02_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")))))
}

s002S = s002
class(s002S) = "numeric"


# Get maximum s-hat of blocks overlapping 5mb (the centre of the chromosome)
for (i in 1:13) {
  for (j in 1:files) {
    if (file.exists(paste(s002[j, i]))) {
      filtered = read.csv(file = paste(s002[j, i]))
      neighbours = which(filtered$bp.end >= 5e6 & filtered$bp.start <= 5e6)
      if (length(neighbours) != 0) {
        s002S[j, i] = max(filtered$sHat[neighbours])
      } else {
        s002S[j, i] = 0
      }
    } else {
      s002S[j, i] = 0
    }
    
  }
}

####### s = 0.05
s005 = matrix(data = NA, nrow = files, ncol = 13)
for (j in 1:length(frequencies)) {
  s005[, j] = c(list.files(pattern = paste("sC0.05_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")), rep(NA, length.out=files-length(list.files(pattern = paste("sC0.05_mF", frequencies[j], ".trees.uniform.vcf.sHat.csv", sep = "")))))
}

s005S = s005
class(s005S) = "numeric"


# Get maximum s-hat of blocks overlapping 5mb (the centre of the chromosome)
for (i in 1:13) {
  for (j in 1:files) {
    if (file.exists(paste(s005[j, i]))) {
      filtered = read.csv(file = paste(s005[j, i]))
      neighbours = which(filtered$bp.end >= 5e6 & filtered$bp.start <= 5e6)
      if (length(neighbours) != 0) {
        s005S[j, i] = max(filtered$sHat[neighbours])
      } else {
        s005S[j, i] = 0
      }
    } else {
      s005S[j, i] = 0
    }
    
  }
}

MAX = max(s00075S, s001S, s002S, s005S)
op = par(no.readonly = TRUE)
par(oma = c(4.7, 2, 1, 1))
par(mfrow = c(4, 1), mai = c(0, 0.4, 0, 0))

first = boxplot.matrix(s00075S,
                      las = 2,
                      xlab = "",
                      col = "grey",
                      ylim = c(0, MAX),
                      at = c(0, 2, 5, 9, 16, 24, 32, 40, 48, 56, 64, 72, 80),
                      boxwex = c(1, 1.5, 2, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6),
                      xaxt = "n",
                      boxlwd = 0.5,
                      medlwd = 1.5,
                      whisklwd = 0.5,
                      staplelwd = 0.5,
                      outlwd = 0.5,
                      outpch = 20,
                      cex = 0.8)
abline(h = 0.0075, lty = 2, col = "red", lwd = 0.5)
abline(v = 1, lty = 1, col = "grey", lwd = 0.5)

second = boxplot.matrix(s001S,
                      las = 2,
                      col = "grey",
                      ylim = c(0, MAX),
                      at = c(0, 2, 5, 9, 16, 24, 32, 40, 48, 56, 64, 72, 80),
                      boxwex = c(1, 1.5, 2, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6),
                      xaxt = "n",
                      boxlwd = 0.5,
                      medlwd = 1.5,
                      whisklwd = 0.5,
                      staplelwd = 0.5,
                      outlwd = 0.5,
                      outpch = 20,
                      cex = 0.8)
abline(h = 0.01, lty = 2, col = "red", lwd = 0.5)
abline(v = 1, lty = 1, col = "grey", lwd = 0.5)

third = boxplot.matrix(s002S,
                       las = 2,
                       col = "grey",
                       ylim = c(0, MAX),
                       at = c(0, 2, 5, 9, 16, 24, 32, 40, 48, 56, 64, 72, 80),
                       boxwex = c(1, 1.5, 2, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6),
                       xaxt = "n",
                       boxlwd = 0.5,
                       medlwd = 1.5,
                       whisklwd = 0.5,
                       staplelwd = 0.5,
                       outlwd = 0.5,
                       outpch = 20,
                       cex = 0.8)
abline(h = 0.02, lty = 2, col = "red", lwd = 0.5)
abline(v = 1, lty = 1, col = "grey", lwd = 0.5)

fourth = boxplot.matrix(s005S,
                        las = 2,
                        col = "grey",
                        ylim = c(0, MAX),
                        at = c(0, 2, 5, 9, 16, 24, 32, 40, 48, 56, 64, 72, 80),
                        boxwex = c(1, 1.5, 2, 4, 6, 6, 6, 6, 6, 6, 6, 6, 6),
                        xaxt = "n",
                        boxlwd = 0.5,
                        medlwd = 1.5,
                        whisklwd = 0.5,
                        staplelwd = 0.5,
                        outlwd = 0.5,
                        outpch = 20,
                        cex = 0.8)
abline(h = 0.05, lty = 2, col = "red", lwd = 0.5)
abline(v = 1, lty = 1, col = "grey", lwd = 0.5)

options(scipen = 10)

inch = function(mm) {
  mm * 0.0393701
  }

######################
cex = 1
pdf("Fig1a.pdf",
    pointsize = 7,
    height = inch(140),
    width = inch(180))

xval = c(0, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)

plot(xval, (1 - abs(0.0075 - first$stats[3, ]) / 0.0075), type = "l", col = cbcol[1], xlim = c(0, 100), ylim = c(0,1), ylab = "", xaxt = "n", xlab = "")
lines(xval, (1 - abs(0.01 - second$stats[3, ]) / 0.01), col = cbcol[2])
lines(xval, (1 - abs(0.02 - third$stats[3, ]) / 0.02), col = cbcol[3])
lines(xval, (1 - abs(0.05 - fourth$stats[3, ]) / 0.05), col = cbcol[4])

mtext("a", side = 3, at = -15, font = 2)

legend("topleft", title = "Selection coefficient", legend = c("0.0075", "0.01", "0.02", "0.05"), col = c(cbcol[1:4]), lty = c(1, 1, 1, 1), xjust = 0, bg = "grey92", box.col = "black", cex = cex)


mtext(TeX("$1 - \\frac{s-\\tilde{\\hat{s}}}{s}"), side = 2, outer = TRUE, line = -2.5, padj=0.5, at=0.8)

mtext(TeX("Inferred selection coefficient $\\hat{s}"), side = 2, outer = TRUE, line = -2.5, padj=0.5, at=0.35)

axis(side = 1,
  at = c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10),
  labels = c("0%", "10%", "20%", "30%", "40%", "50%", "60%", "70%", "80%", "90%", "100%"))

axis(side = 2, at = seq(0, 0.04, by = 0.01))

mtext(TeX("Frequency of selected variant in population"), side = 1, outer = TRUE, line = -3.5, adj = 0.5)
mtext("b", side = 3, at = -1.5, font = 2)

dev.off()