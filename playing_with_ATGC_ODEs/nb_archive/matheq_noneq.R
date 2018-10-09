###############################################################################
get.complementar <- function(seq="GCATTCCACA"){
  seq <- unlist(strsplit(seq,""))
  seq2 <- seq
  seq2[which(seq=="A")] <- "T"
  seq2[which(seq=="G")] <- "C"
  seq2[which(seq=="T")] <- "A"
  seq2[which(seq=="C")] <- "G"
  return(paste(seq2,collapse=""))
}
###############################################################################
###############################################################################
get.reverse <- function(seq="GCATTCCACA"){
  seq <- unlist(strsplit(seq,""))
  return(paste(rev(seq),collapse=""))
}
###############################################################################




# This script generates the equations for the hypercube system.
# source("/Users/sirius/GIT/TrantoR/GEN_SinglePointMutate.R")
library("gtools")
permuts <- permutations(n=4, r=2, v=c("A","G","T","C"), repeats.allowed=TRUE)
dibases <- sapply(1:16, function(i){paste(permuts[i,],collapse="")})

k <- matrix(NA, ncol=16, nrow=16) # row=FROM, column=TO
dimnames(k)[[1]] <- dimnames(k)[[2]] <- dibases

I <- J <- "dummyChar"

count <- 1
for(i in 1:16){
  for(j in 1:16){
    if( (dibases[i] != dibases[j]) &
        ((permuts[i,1]==permuts[j,1]) | (permuts[i,2]==permuts[j,2])) ){
      db.from <- dibases[i]
      cp.from <- get.reverse(get.complementar(db.from))
      db.to <- dibases[j]
      cp.to <- get.reverse(get.complementar(db.to))
      i.cp <- which(dibases==cp.from)
      j.cp <- which(dibases==cp.to)
      if(sum(I==i.cp & J==j.cp)==0 & sum(I==i & J==j)==0){
        k[i,j] <- paste("k", dibases[i],"2",dibases[j], sep="")
        I <- c(I, i); J <- c(J, j)
        k[i.cp, j.cp] <- paste("k", dibases[i],"2",dibases[j], sep="")
        I <- c(I, i.cp); J <- c(J, j.cp)
        count <- count+1
      }

    }

  }
}
# symbolic rate matrix is generated


EQS <- NULL
for(p in 1:length(permuts[,1])){
  di <- permuts[p, ]; di.joint <- paste(di, collapse="")
  b1 <- di[1]
  b2 <- di[2]
  eq <- "0 == " #paste("C",di.joint,"'[t] == ", sep="")

  # Adding all the cases where the second base is set constant, hence varying
  # the first base only.
  to <- b1
  from <- c("A","C","T","G")[which(c("A","C","T","G")!=to)]
  term <- 1
  for(fr in from){
    if(term==1){
      eq <- paste(eq,      k[paste(fr,b2,sep=""),paste(to,b2,sep="")],"*C",fr,b2,
                     " - ",k[paste(to,b2,sep=""),paste(fr,b2,sep="")],"*C",to,b2,sep="")
    } else {
      eq <- paste(eq," + ",k[paste(fr,b2,sep=""),paste(to,b2,sep="")],"*C",fr,b2,
                     " - ",k[paste(to,b2,sep=""),paste(fr,b2,sep="")],"*C",to,b2,sep="")
    }
    term <- 0
  }

  to <- b2
  from <- c("A","C","T","G")[which(c("A","C","T","G")!=to)]
  for(fr in from){
    eq <- paste(eq," + ",k[paste(b1,fr,sep=""),paste(b1,to,sep="")],"*C",b1,fr,
                   " - ",k[paste(b1,to,sep=""),paste(b1,fr,sep="")],"*C",b1,to,sep="")
  }

  EQS <- c(EQS, eq)

}

# paste(sapply(1:16, function(i){paste("C",paste(permuts[i,],collapse=""),sep="")}),collapse=" + ")


SinglePointMutate(sequence="AG", seqtype="DNA", mutation.position=2)


















EQS <- NULL
for(p in 1:length(permuts[,1])){
  di <- permuts[p, ]; di.joint <- paste(di, collapse="")
  b1 <- di[1]
  b2 <- di[2]
  eq <- paste("C",di.joint,"'[t] == ", sep="")


  EQS <- c(EQS, eq)

}

# paste(sapply(1:16, function(i){paste("C",paste(permuts[i,],collapse=""),sep="")}),collapse=" + ")








