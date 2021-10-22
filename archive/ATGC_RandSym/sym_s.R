load("solsym.Rdata")


pdf(width=7, height=7, file="GC_hist.pdf")
hist(solsym$GC, col="cornflowerblue", breaks=50, freq=FALSE,
     xlab="G+C content, %", ylab="density",
     main=paste("strand-symmetric, t[0]=25/25/25/25, mean = ",round(mean(solsym$GC),2),", sd = ",round(sd(solsym$GC),2),sep=""))
dev.off()


plot(y=solsym$GC, x=solsym$kag, cex=0.3)




kag = ktc
kat = kta
kac = ktg
kct = kga
kca = kgt
kcg = kgc




require(geneplotter)
colfun   = colorRampPalette(c("white","blue","skyblue",
                              "chartreuse3","green","yellow","orange","red","darkred"))


################################################################################
# Screenshotted for the better quality of the colours!
par(mfrow=c(3,2))

y=solsym$GC
x=solsym$kag
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_AG&TC, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solsym$GC
x=solsym$kac
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_AC&TG, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solsym$GC
x=solsym$kct
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_CT&GA, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solsym$GC
x=solsym$kca
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_CA&GT, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solsym$GC
x=solsym$kat
smoothScatter(y=y,
              x=x,
              nrpoints=100, nbin=400,bandwidth=c(diff(range(x))/200, diff(range(y))/200),
              xlab="k_AT&TA, byr-1", ylab="G+C content, %",
              colramp=colfun, main=NULL)

y=solsym$GC
x=solsym$kcg
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








