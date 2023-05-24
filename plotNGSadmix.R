# plotar resultados NGSadmix (adaptado direto do site)
# genotropics workshop 2023

setwd("G:/OneDrive - PUCRS - BR/cursos/genotropics2023")

# Reading arguments
qopt <- "leopardus_63ind_chrA1_minInd21_k2.qopt" # .qopt file name
k <- 2 # Number of K from input qopt, ex. 2
metadata <- "leopardus_63ind_chrA1.popinfo" # Metadata of the population 

# Get ID and pop info for each individual
name <- gsub(".qopt","",qopt)
popinfo <- read.csv(metadata, sep = ";", header = TRUE)
ind <- popinfo$ID
pop <- popinfo$POP
reg <- popinfo$REG
lat <- popinfo$LAT
lon <- popinfo$LON
latlon <- as.numeric(popinfo$LAT)*as.numeric(popinfo$LON)

# Read inferred admixture proportions file
q <- read.table(qopt)

# Plot the results ordered by q val
png(filename=paste(name, ".ordq.png", sep = ""), width=1000, height=400)
ylab=paste("Admixture proportions (K=",k,")", sep = "")
main=paste(name,"(ordered by q value)", sep = " ")
par(mar=c(5,4,1,1))
barplot(t(q)[,order(q$V2)],
        col=c('darkred','blue4','orange','darkgreen','black','gray','pink'),
        space=0,
        #border=NA,
        names=ind[order(q$V2)],
        las=2,
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
dev.off()

# Plot the results ordered by population name
png(filename=paste(name, ".ordp.png", sep = ""), width=1000, height=400)
main=paste(name,"(ordered by populations name)", sep = " ")
par(mar=c(9,4,1,1))
barplot(t(q)[,order(pop)],
        col=c('darkred','blue4','darkgreen','black','orange','gray','pink'),
        space=0,
        border=NA,
        las=2,
        names=ind[order(pop)],
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
text(tapply(0.5:length(pop),pop[order(pop)],mean),-0.3,unique(pop[order(pop)]),xpd=T, srt=90, cex=0.75)
abline(v=cumsum(sapply(unique(pop[order(pop)]),function(x){sum(pop[order(pop)]==x)})),col=1,lwd=1)
dev.off()

# Plot the results ordered by region name
png(filename=paste(name, ".ordr.png", sep = ""), width=1000, height=400)
main=paste(name,"(ordered by regions name)", sep = " ")
par(mar=c(9,4,1,1))
barplot(t(q)[,order(reg)],
        col=c('darkred','blue4','black','orange','darkgreen','gray','pink'),
        space=0,
        border=NA,
        las=2,
        names=ind[order(reg)],
        ylab=ylab,
        main=main,
        cex.axis=1, 
        cex.names=0.7)
text(tapply(0.5:length(reg),reg[order(reg)],mean),-0.3,unique(reg[order(reg)]),xpd=T, srt=90, cex=0.75)
abline(v=cumsum(sapply(unique(reg[order(reg)]),function(x){sum(reg[order(reg)]==x)})),col=1,lwd=1)
dev.off()
