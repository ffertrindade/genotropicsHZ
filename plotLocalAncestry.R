#!/usr/bin/env Rscript
# plotar ancestralidade local ao longo de varios cromossomos
# usa como input arquivo gerado pelo script getAvrLAI.py
# script para rodar por linha de comando
# genotropics workshop 2023

### command arguments, libraries and wd
args = commandArgs(trailingOnly=TRUE)

if (length(args)<2) {
  stop("Usage: plotLocalAncestry.R id path", call.=FALSE)
}

id <- args[1]
path <- args[2]

if(require(doMC)){
  registerDoMC(5)
} else {
  registerDoSEQ()
}

setwd(path)

### data input
file <- paste(id, ".posterior.parent1.csv", sep = '')
individual <- read.table(file, header=TRUE, sep = "\t")
ylab = "parent1 posterior"

color <- sapply(individual, is.factor)
individual[color] <- lapply(individual[color], as.character)

#head(print(individual))
### plotting geoffroyi posterior
file <- paste(id, ".posterior.parent1.png", sep = "")
png(filename=file, width=900, height=200)
plot( x = individual$pos,
      y = individual$avr_parent1,
      ylab = ylab,
      xlab = "SNPs",
      main = id,
      col = individual$color,
      ylim = c(0,1),
      xaxt="n",
      pch = '.')
axis(2, at=0:2)
dev.off()
