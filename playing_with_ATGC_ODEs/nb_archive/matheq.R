# This script generates the equations for the hypercube system.
# source("/Users/sirius/GIT/TrantoR/GEN_SinglePointMutate.R")
library("gtools")


kCA <- kGT <- "i"
kAC <- kTG <- "j"
kCG <- kGC <- "k"
kAT <- kTA <- "l"
kCT <- kGA <- "m"
kAG <- kTC <- "n"





permuts <- permutations(n=4, r=2, v=c("A","G","T","C"), repeats.allowed=TRUE)


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
      eq <- paste(eq,      eval(parse(text=paste("k",fr,to,sep=""))),"*C",fr,b2,
                     " - ",eval(parse(text=paste("k",to,fr,sep=""))),"*C",to,b2,sep="")
    } else {
      eq <- paste(eq," + ",eval(parse(text=paste("k",fr,to,sep=""))),"*C",fr,b2,
                     " - ",eval(parse(text=paste("k",to,fr,sep=""))),"*C",to,b2,sep="")
    }
    term <- 0
  }

  to <- b2
  from <- c("A","C","T","G")[which(c("A","C","T","G")!=to)]
  for(fr in from){
    eq <- paste(eq," + ",eval(parse(text=paste("k",fr,to,sep=""))),"*C",b1,fr,
                   " - ",eval(parse(text=paste("k",to,fr,sep=""))),"*C",b1,to,sep="")
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








