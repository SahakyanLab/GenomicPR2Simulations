load("solnosym.Rdata")

solnosym$C <- (solnosym$GC/100)/(1+solnosym$GCratio)
solnosym$G <- (solnosym$GC/100)-solnosym$C
solnosym$T <- (1-(solnosym$GC/100)) / (1+solnosym$ATratio)
solnosym$A <- (1-(solnosym$GC/100)) - solnosym$T

pdf(width=7, height=7, file="nosym_GC_hist.pdf")
hist(solnosym$GC, col="cornflowerblue", breaks=50, freq=FALSE,
     xlab="G+C content, %", ylab="density",
     main=paste("non-symmetric, t[0]=25/25/25/25, mean = ",round(mean(solnosym$GC),2),", sd = ",round(sd(solnosym$GC),2),sep=""))
dev.off()

#require(geneplotter)
#colfun = colorRampPalette(c("white","blue","skyblue",
#                            "chartreuse3","green","yellow",
#                            "orange","red","darkred"))
#plot(y=solnosym$GC, x=solnosym$kag, cex=0.3)
#plot(y=solnosym$GCratio, x=solnosym$ATratio, cex=0.3)

#kag = ktc
#kat = kta
#kac = ktg
#kct = kga
#kca = kgt
#kcg = kgc

require(geneplotter)
colfun = colorRampPalette(c("white","blue","skyblue",
                            "chartreuse3","green","yellow",
                            "orange","red","darkred"))

y = solnosym$G-solnosym$C
x = solnosym$A-solnosym$T

pdf(width=6, height=6, file="nosym_G-C_2Dhist.pdf")
smoothScatter(y=y, x=x,
              nrpoints=100, nbin=400,
              bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="A-T", ylab="G-C",
              xlim=c(-1,1), ylim=c(-1,1),
              colramp=colfun, main=NULL)
dev.off()


#y = solnosym$GCratio
#x = solnosym$ATratio
#smoothScatter(y=y, x=x,
#              nrpoints=100, nbin=400,
#              bandwidth=c(diff(range(x))/200, diff(range(y))/700),
#              xlab="AT ratio", ylab="GC ratio",
#              xlim=c(0,60), ylim=c(0,60),
#              colramp=colfun, main=NULL)
#
#smoothScatter(y=y, x=x,
#              nrpoints=100, nbin=400,
#              bandwidth=c(diff(range(x))/200, diff(range(y))/1900),
#              xlab="AT ratio", ylab="GC ratio",
#              xlim=c(0,10), ylim=c(0,10),
#              colramp=colfun, main=NULL)


#kag = ktc
#kat = kta
#kac = ktg
#kct = kga
#kca = kgt
#kcg = kgc

ind <- 1:length(solnosym[,1])
pdf(height=13, width=5, file="nosym_allind.pdf")
par(mfrow=c(6,1))

hist(solnosym$kag[ind]/solnosym$ktc[ind], main="kAG/kTC", xlab="Value",
     breaks=15000000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kat[ind]/solnosym$kta[ind], main="kAT/kTA", xlab="Value",
     breaks=5000000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kac[ind]/solnosym$ktg[ind], main="kAC/kTG", xlab="Value",
     breaks=3500000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kct[ind]/solnosym$kga[ind], main="kCT/kGA", xlab="Value",
     breaks=1000000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kca[ind]/solnosym$kgt[ind], main="kCA/kGT", xlab="Value",
     breaks=1500000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kcg[ind]/solnosym$kgc[ind], main="kCG/kGC", xlab="Value",
     breaks=5000000, xlim=c(0,5), freq=FALSE, col="skyblue")

dev.off()


y = solnosym$G-solnosym$C
x = solnosym$A-solnosym$T
ind <- which(abs(y)<0.01 & abs(x)<0.01)

pdf(height=13, width=5, file="nosym_pm001_subtr.pdf")
par(mfrow=c(6,1))

hist(solnosym$kag[ind]/solnosym$ktc[ind], main="kAG/kTC", xlab="Value",
     breaks=3000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kat[ind]/solnosym$kta[ind], main="kAT/kTA", xlab="Value",
     breaks=4000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kac[ind]/solnosym$ktg[ind], main="kAC/kTG", xlab="Value",
     breaks=2000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kct[ind]/solnosym$kga[ind], main="kCT/kGA", xlab="Value",
     breaks=12000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kca[ind]/solnosym$kgt[ind], main="kCA/kGT", xlab="Value",
     breaks=11000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kcg[ind]/solnosym$kgc[ind], main="kCG/kGC", xlab="Value",
     breaks=1500, xlim=c(0,5), freq=FALSE, col="skyblue")

dev.off()



## RATIO #######################################################################
ind <- which(solnosym$GCratio<1.05 & solnosym$GCratio>0.95 &
      solnosym$ATratio<1.05 & solnosym$ATratio>0.95)
ind <- (1:length(solnosym[,1]))[-ind]

pdf(height=13, width=5, file="nosym_else.pdf")
par(mfrow=c(6,1))

hist(solnosym$kag[ind]/solnosym$ktc[ind], main="kAG/kTC", xlab="Value",
     breaks=15000000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kat[ind]/solnosym$kta[ind], main="kAT/kTA", xlab="Value",
     breaks=5000000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kac[ind]/solnosym$ktg[ind], main="kAC/kTG", xlab="Value",
     breaks=3500000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kct[ind]/solnosym$kga[ind], main="kCT/kGA", xlab="Value",
     breaks=1000000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kca[ind]/solnosym$kgt[ind], main="kCA/kGT", xlab="Value",
     breaks=1500000, xlim=c(0,5), freq=FALSE, col="skyblue")

hist(solnosym$kcg[ind]/solnosym$kgc[ind], main="kCG/kGC", xlab="Value",
     breaks=5000000, xlim=c(0,5), freq=FALSE, col="skyblue")

dev.off()
## RATIO #######################################################################






smoothScatter(y=solnosym$GCratio,
              x=solnosym$ATratio,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/700),
              xlab="AT ratio", ylab="GC ratio",
              xlim=c(0,60), ylim=c(0,60),
              colramp=colfun, main=NULL)



################################################################################
# Screenshotted for the better quality of the colours!
par(mfrow=c(3,2))

y=solnosym$GC
x=solnosym$kag
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_AG&TC, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solnosym$GC
x=solnosym$kac
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_AC&TG, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solnosym$GC
x=solnosym$kct
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_CT&GA, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solnosym$GC
x=solnosym$kca
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_CA&GT, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solnosym$GC
x=solnosym$kat
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_AT&TA, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solnosym$GC
x=solnosym$kcg
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_CG&GC, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)
################################################################################



















#pdf(height=17, width=4, file=plotfile)
jpeg(height=1700/1.5, width=500/1.5, file=plotfile, quality=99)
par(mfrow=c(6,1))

#$$$$$$$$$$$$
for(k in ks){
  kmerDB = paste(kmerDB.path,"GenomeKmer_",
                 k,"/chr_ALL_kmer_",k,"_tbl.txt",sep="")
  kmer.data <- read.table(kmerDB, header=TRUE)

  # Determining the flanks around the center of the segments to be extracted from mut.seq11:
  flank <- (k-1)/2
  mut.seqk <- sapply(mut.seq11, USE.NAMES=FALSE, simplify=TRUE,
                     FUN=function(i){substr(i, start=(6-flank), stop=(6+flank))} )
  counts <- kmer.data[match(mut.seqk, kmer.data[,"Var1"]), "Freq"]

  ##############################################################################
  smoothScatter(y=counts,
                x=1000*mut.rate,
                nrpoints=100, nbin=700,
                xlab="mutation rate, mut/byr", ylab="k-mer count",
                colramp=colfun, main=paste("k = ",k,sep=""))
  ##############################################################################
}
#$$$$$$$$$$$$

dev.off()











CAA = 134408249/1366670324
CAC = 67825715/1366670324
CAG = 96146324/1366670324
CAT = 108265406/1366670324
CCA = 97345439/1366670324
CCC = 69625606/1366670324
CCG = 13066878/1366670324
CCT = 96176526/1366670324
CGA = 81863001/1366670324
CGC = 56875422/1366670324
CGG = 69709260/1366670324
CGT = 68024076/1366670324
CTA = 93207971/1366670324
CTC = 81851281/1366670324
CTG = 97573424/1366670324
CTT = 134705746/1366670324











