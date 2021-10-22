################################################################################
SolveATGC <- function(parameters=parameters, state=state, step=step, span=span,
                      sim.evol = FALSE, 
                      EQtolerance=EQtolerance, CHtolerance=CHtolerance){

  # Function to numerically solve the kinetic rate equations for 
  # the simulation of the genome dynamics evolution

  # Dependencies
  #     SolveATGC.R, InputChecking.R

  # Flag         Format       Description
  # parameters   <numeric>    All mutation rate constants
  # state        <numeric>    All four base contents
  # step         <numeric>    Time step in evaluating the kinetic
  #                           mutation rate equations
  # span         <numeric>    Maximum time period for the time evolution
  #                           of genomic base content
  # sim.evol     <character>  Capture the 1 billion year time evolution of the change
  #                           of base content and meta-data
  # EQtolerance  <numeric>    Tolerance value for reaching genome equilibration
  # CHtolerance  <numeric>    Tolerance value for reaching Chargaff compliance

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
  suppressPackageStartupMessages(library(deSolve))

  out <- ode(y=state, times=times, func=mut.mdl, parms=parameters)

  length.out <- dim(out)[1]
  length.genome <- sum(state)
  dif.nucl <- c(sum(out[1,c("Ca","Cg","Ct","Cc")]),
                rowSums(abs(diff(out[,c("Ca","Cg","Ct","Cc")])))/4)
  at.ratio <- out[,"Ca"]/out[,"Ct"]
  gc.ratio <- out[,"Cg"]/out[,"Cc"]
  at.content <- (out[,"Ca"] + out[,"Ct"])*100/length.genome
  gc.content <- (out[,"Cg"] + out[,"Cc"])*100/length.genome
  at.skew <- ((out[,"Ca"] - out[,"Ct"])/(out[,"Ca"] + out[,"Ct"]))/length.genome
  gc.skew <- ((out[,"Cg"] - out[,"Cc"])/(out[,"Cg"] + out[,"Cc"]))/length.genome
  Fin.GC  <- as.vector(gc.content[length.out])

  # 1bn year time intervals
  when.one.bn <- which(out[,"time"] == 1)
  when.two.bn <- which(out[,"time"] == 2)
  when.three.bn <- which(out[,"time"] == 3)
  when.four.bn <- which(out[,"time"] == 4)

  RESULTS <- NULL
  RESULTS$inp                 <- NULL # input data
  RESULTS$inp$parameters      <- parameters
  RESULTS$inp$state           <- state
  RESULTS$inp$step            <- step
  RESULTS$inp$span            <- span
  RESULTS$out                 <- out                 # numerical solution to the diff eqs
  RESULTS$dif.nucl            <- dif.nucl
  RESULTS$fluctuation.mean    <- mean(dif.nucl[2:length(dif.nucl)]/4) # average fluctuation difference
  RESULTS$fluctuation.sd      <- sd(dif.nucl[2:length(dif.nucl)]/4) # st.dev. fluctuation difference
  RESULTS$at.ratio            <- as.vector(at.ratio) # at.ratio dynamics vector
  RESULTS$gc.ratio            <- as.vector(gc.ratio)
  RESULTS$gc.ratio            <- as.vector(gc.ratio)
  RESULTS$at.content          <- as.vector(at.content)
  RESULTS$gc.content          <- as.vector(gc.content)
  RESULTS$at.skew             <- as.vector(at.skew)
  RESULTS$gc.skew             <- as.vector(gc.skew)
  RESULTS$Fin.GC              <- Fin.GC            # the G+C content at the end of the simulation
  RESULTS$length.genome       <- length.genome
  RESULTS$length.out          <- length.out

  # record each billion year time step to observe evolution
  if(sim.evol){
    # 1bn yrs
    RESULTS$one.bn              <- NULL # 1bn years
    RESULTS$one.bn$at.content   <- at.content[when.one.bn]
    RESULTS$one.bn$gc.content   <- gc.content[when.one.bn]
    RESULTS$one.bn$at.skew      <- at.skew[when.one.bn]
    RESULTS$one.bn$gc.skew      <- gc.skew[when.one.bn]
    
    # 2bn yrs
    RESULTS$two.bn              <- NULL # 2bn years
    RESULTS$two.bn$at.content   <- at.content[when.two.bn]
    RESULTS$two.bn$gc.content   <- gc.content[when.two.bn]
    RESULTS$two.bn$at.skew      <- at.skew[when.two.bn]
    RESULTS$two.bn$gc.skew      <- gc.skew[when.two.bn]
    
    # 3bn yrs
    RESULTS$three.bn            <- NULL # 3bn years
    RESULTS$three.bn$at.content <- at.content[when.three.bn]
    RESULTS$three.bn$gc.content <- gc.content[when.three.bn]
    RESULTS$three.bn$at.skew    <- at.skew[when.three.bn]
    RESULTS$three.bn$gc.skew    <- gc.skew[when.three.bn]
    # 4bn yrs
    RESULTS$four.bn             <- NULL # 4bn years
    RESULTS$four.bn$at.content  <- at.content[when.four.bn]
    RESULTS$four.bn$gc.content  <- gc.content[when.four.bn]
    RESULTS$four.bn$at.skew     <- at.skew[when.four.bn]
    RESULTS$four.bn$gc.skew     <- gc.skew[when.four.bn]
  }

  # Genome Equilibrium tolerance
  if(EQtolerance != FALSE){
    dif <- abs(diff(out[,c("Ca","Cg","Ct","Cc")]))
    dif.max.min <- apply(dif, 1, function(x){max(x) - min(x)})
    Eq.time <- times[which(dif.max.min <= EQtolerance)[1]]

    # Chargaff compliance
    at.tol <- CHtolerance[CHtolerance$metadata == "AT_skew", "st.dev"]
    gc.tol <- CHtolerance[CHtolerance$metadata == "GC_skew", "st.dev"]
    Ch.compliance <- which(abs(gc.skew-0) <= gc.tol & abs(at.skew-0) <= at.tol)

    # if differences between compliant cases is 1, take first index
    # otherwise, take first index after the biggest difference
    if(length(table(diff(Ch.compliance))) == 1){
      Ch.time <- times[Ch.compliance[which(diff(Ch.compliance) == 1)[1]]]
    } else {
      Ch.time <- times[Ch.compliance[tail(which(diff(Ch.compliance) != 1), n = 1)+1]]
    }

    if(length(Ch.time) == 0){
      Ch.time <- NA
    }

    # time difference in equilibration and chargaff compliance
    Eq.Ch.time.diff <- Eq.time-Ch.time

    RESULTS$inp$CH.AT.tolerance <- at.tol
    RESULTS$inp$CH.GC.tolerance <- gc.tol
    RESULTS$inp$EQtolerance     <- EQtolerance
    RESULTS$Eq.time             <- Eq.time           # time to reach the equilibration tolerance limit
    RESULTS$Ch.time             <- Ch.time           # time to reach Chargaff's tolerance limit
    RESULTS$Eq.Ch.time.diff     <- Eq.Ch.time.diff   # time difference in equilibration and chargaff compliance
  }
  return(RESULTS)
}
################################################################################