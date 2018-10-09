source("fun.R")
library("deSolve")
library("foreach")
library("doMC")

GenomeSize = 3000000000 # nt  # Irrelevant!
Acont      = 25 # %
Gcont      = 25 # %
Ccont      = 25 # %
span       = 4.0   # byr
step       = 0.001 # byr
muttype    = "Strand-Symmetric" # "Non Symmetric"
minr = 0
maxr = 1.5
NCPU = 12


registerDoMC( cores = NCPU )
Tcont <- 100 - Acont - Gcont - Ccont
state <- c(Ca=GenomeSize*Acont/100,
           Cg=GenomeSize*Gcont/100,
           Ct=GenomeSize*Tcont/100,
           Cc=GenomeSize*Ccont/100)


### foreach ##############################################
solsym <- foreach(i=1:500000, .combine="rbind", .inorder=FALSE)%dopar%{

if(muttype == "Non Symmetric"){
  rates <- runif(n=12, min=minr, max=maxr)
  kag = rates[1]
  kat = rates[2]
  kac = rates[3]
  kga = rates[4]
  kgt = rates[5]
  kgc = rates[6]
  kta = rates[7]
  ktg = rates[8]
  ktc = rates[9]
  kca = rates[10]
  kcg = rates[11]
  kct = rates[12]
}


if(muttype == "Strand-Symmetric"){
  rates <- runif(n=6, min=minr, max=maxr)
  kag = ktc = rates[1]
  kat = kta = rates[2]
  kac = ktg = rates[3]
  kct = kga = rates[4]
  kca = kgt = rates[5]
  kcg = kgc = rates[6]
}

parameters <- c(kag=kag, kat=kat, kac=kac,
                kga=kga, kgt=kgt, kgc=kgc,
                kta=kta, ktg=ktg, ktc=ktc,
                kca=kca, kcg=kcg, kct=kct) #mut/byr
                
atgc <- SolveATGC(parameters = parameters, state = state, step = step, span = span,
                  EQtolerance = 1, CHtolerance = 0.001)

return(data.frame(GC=atgc$gc.content[atgc$length.out],
       GCratio=atgc$gc.ratio[atgc$length.out],
       ATratio=atgc$at.ratio[atgc$length.out],
       kag=kag, kat=kat, kac=kac,
       kga=kga, kgt=kgt, kgc=kgc,
       kta=kta, ktg=ktg, ktc=ktc,
       kca=kca, kcg=kcg, kct=kct)
      )

}
### end of foreach ######################################
save(solsym, file="solsym.Rdata")




muttype    = "Non Symmetric"
### foreach ##############################################
solnosym <- foreach(i=1:500000, .combine="rbind", .inorder=FALSE)%dopar%{

if(muttype == "Non Symmetric"){
  rates <- runif(n=12, min=minr, max=maxr)
  kag = rates[1]
  kat = rates[2]
  kac = rates[3]
  kga = rates[4]
  kgt = rates[5]
  kgc = rates[6]
  kta = rates[7]
  ktg = rates[8]
  ktc = rates[9]
  kca = rates[10]
  kcg = rates[11]
  kct = rates[12]
}


if(muttype == "Strand-Symmetric"){
  rates <- runif(n=6, min=minr, max=maxr)
  kag = ktc = rates[1]
  kat = kta = rates[2]
  kac = ktg = rates[3]
  kct = kga = rates[4]
  kca = kgt = rates[5]
  kcg = kgc = rates[6] 
}

parameters <- c(kag=kag, kat=kat, kac=kac,
                kga=kga, kgt=kgt, kgc=kgc,
                kta=kta, ktg=ktg, ktc=ktc,
                kca=kca, kcg=kcg, kct=kct) #mut/byr

atgc <- SolveATGC(parameters = parameters, state = state, step = step, span = span,
                  EQtolerance = 1, CHtolerance = 0.001)

return(data.frame(GC=atgc$gc.content[atgc$length.out],
       GCratio=atgc$gc.ratio[atgc$length.out],
       ATratio=atgc$at.ratio[atgc$length.out],
       kag=kag, kat=kat, kac=kac,
       kga=kga, kgt=kgt, kgc=kgc,
       kta=kta, ktg=ktg, ktc=ktc,
       kca=kca, kcg=kcg, kct=kct)
      )

}
### end of foreach ######################################
save(solnosym, file="solnosym.Rdata")

