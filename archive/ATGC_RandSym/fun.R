################################################################################
SolveATGC <- function(parameters=parameters, state=state, step=step, span=span,
                      EQtolerance=EQtolerance, CHtolerance=CHtolerance){

  times <- seq(0, span, by=step) # byr

  ##########################################
  mut.mdl <- function(t, state, parameters){

   with(as.list(c(state, parameters)),{
   # rate of change
     dCa <- kca*Cc+kta*Ct+kga*Cg-(kac+kat+kag)*Ca
     dCg <- kag*Ca+kcg*Cc+ktg*Ct-(kga+kgt+kgc)*Cg
     dCt <- kat*Ca+kgt*Cg+kct*Cc-(kta+ktc+ktg)*Ct
     dCc <- kac*Ca+ktc*Ct+kgc*Cg-(kca+kct+kcg)*Cc

     # return the rate of change
     list(c(dCa, dCg, dCt, dCc))
   }) # end with(as.list ...

  }
  ##########################################

  library(deSolve)
  out <- ode(y=state, times=times, func=mut.mdl, parms=parameters)
  length.out <- length(out[,1])

  length.genome <-  sum(state)


  dif.nucl <- c( sum( out[1,c("Ca","Cg","Ct","Cc")] ),
                sapply(2:length.out, function(i){
                  sum(abs(out[i,c("Ca","Cg","Ct","Cc")] - out[(i-1),c("Ca","Cg","Ct","Cc")]))
                  }, USE.NAMES=FALSE, simplify=TRUE))

  at.ratio <- sapply(1:length.out, function(i){
                  out[i,"Ca"]/out[i,"Ct"]
                  }, USE.NAMES=FALSE, simplify=TRUE)

  gc.ratio <- sapply(1:length.out, function(i){
                  out[i,"Cg"]/out[i,"Cc"]
                  }, USE.NAMES=FALSE, simplify=TRUE)

  gc.content <- sapply(1:length.out, function(i){
                  (out[i,"Cg"] + out[i,"Cc"])*100/length.genome 
                  }, USE.NAMES=FALSE, simplify=TRUE)

  Eq.time <- times[suppressWarnings(min(which((dif.nucl<=EQtolerance)==TRUE)))]
  Ch.time <- times[suppressWarnings(
      min(which((abs(at.ratio-1)<=CHtolerance)&(abs(gc.ratio-1)<=CHtolerance)))
                 )]
  Fin.GC  <- as.vector(gc.content[length.out])

  RESULTS <- NULL
  RESULTS$inp              <- NULL # input data
  RESULTS$inp$parameters   <- parameters
  RESULTS$inp$state        <- state     
  RESULTS$inp$step         <- step    
  RESULTS$inp$span         <- span      
  RESULTS$inp$EQtolerance  <- EQtolerance
  RESULTS$inp$CHtolerance  <- CHtolerance 
  RESULTS$out              <- out                 # numerical solution to the diff eqs
  RESULTS$dif.nucl         <- dif.nucl
  RESULTS$at.ratio         <- as.vector(at.ratio) # at.ratio dynamics vector
  RESULTS$gc.ratio         <- as.vector(gc.ratio)
  RESULTS$gc.content       <- as.vector(gc.content)
  RESULTS$Eq.time          <- Eq.time  # time to reach the equilibration tolerance limit
  RESULTS$Ch.time          <- Ch.time  # time to reach Chargaff's tolerance limit
  RESULTS$Fin.GC           <- Fin.GC   # the G+C content at the end of the simulation
  RESULTS$length.genome    <- length.genome
  RESULTS$length.out       <- length.out

  return(RESULTS)
  
}
################################################################################

################################################################################
PlotATGC <- function(atgc=atgc, time.unit="byr", xlim=c(0,4), ylim=c(0,60)){

  plot(x=atgc$out[,"time"],
       y=100*atgc$out[, "Cg"]/atgc$length.genome,
       ylab="base content, %", xlab=paste("time, ",time.unit,sep=""),
       main=paste("genome length = ", atgc$length.genome, sep=" "),
       pch = 1, cex=0.01, type="l", lwd=3, col="orange",
       xlim=xlim, ylim=ylim)
      
  lines(x=atgc$out[,"time"],
        y=100*atgc$out[, "Ca"]/atgc$length.genome,
        pch = 1, cex=0.01, type="l", lwd=3, col="forestgreen")

  lines(x=atgc$out[,"time"],
        y=100*atgc$out[, "Cc"]/atgc$length.genome,
        pch = 1, cex=0.01, type="l", lwd=3, col="blue")

  lines(x=atgc$out[,"time"],
        y=100*atgc$out[, "Ct"]/atgc$length.genome,
        pch = 1, cex=0.01, type="l", lwd=3, col="red")
      
  text(x=min(xlim)+diff(xlim)/15, y=max(ylim), paste("t = 0 ",time.unit,"\n",
                                      "nG = ",atgc$inp$state["Cg"],"\n",
                                      "nC = ",atgc$inp$state["Cc"],"\n",
                                      "nA = ",atgc$inp$state["Ca"],"\n",
                                      "nT = ",atgc$inp$state["Ct"],"\n",
                                      "G+C = ",round(100*(atgc$inp$state["Cg"]+atgc$inp$state["Cc"])/
                                                    atgc$length.genome,2),"%","\n",
                                      "G/C = ",format(round(atgc$inp$state["Cg"]/atgc$inp$state["Cc"],3),nsmall=3),"\n",   
                                      "A/T = ",format(round(atgc$inp$state["Ca"]/atgc$inp$state["Ct"],3),nsmall=3), sep=""),             
                                      pos = 1, cex=0.6)

  text(x=max(xlim)-diff(xlim)/15, y=max(ylim), paste("t = ",round(atgc$inp$span,2)," ",time.unit,"\n",
                                      "nG = ",round(as.vector(atgc$out[atgc$length.out,"Cg"]),0),"\n",
                                      "nC = ",round(as.vector(atgc$out[atgc$length.out,"Cc"]),0),"\n",
                                      "nA = ",round(as.vector(atgc$out[atgc$length.out,"Ca"]),0),"\n",
                                      "nT = ",round(as.vector(atgc$out[atgc$length.out,"Ct"]),0),"\n",
                                      "G+C = ",round(atgc$gc.content[atgc$length.out],2),"%","\n",
                                      "G/C = ",format(round(atgc$out[atgc$length.out,"Cg"]/atgc$out[atgc$length.out,"Cc"],3),nsmall=3),"\n",   
                                      "A/T = ",format(round(atgc$out[atgc$length.out,"Ca"]/atgc$out[atgc$length.out,"Ct"],3),nsmall=3), sep=""),             
                                      pos = 1, cex=0.6)
  if(!is.na(atgc$Ch.time)){
    abline(v=atgc$Ch.time, lty="dashed")
    text(y=min(ylim),  x=atgc$Ch.time, pos=4, "Chargaff eq. reached", las=1, cex=0.6, srt=90)  
  }

  if(!is.na(atgc$Eq.time)){
    abline(v=atgc$Eq.time)
    text(y=min(ylim),  x=atgc$Eq.time, pos=4, "Genome eq. reached", las=1, cex=0.6, srt=90)  
  }
  
}
################################################################################
