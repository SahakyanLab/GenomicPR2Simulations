################################################################################
solsym <- function(
  Acont      = 0.25,
  Gcont      = 0.25,
  Ccont      = 0.25,
  span       = 4.28,
  step       = 0.001,
  max.runs   = 5, 
  muttype    = "Strand_Symmetric",
  dist       = "normal", 
  species    = "NONE", 
  scale.fac  = 1,
  tolerance  = FALSE,
  tol.return = "NONE",
  sim.evol   = FALSE, 
  NCPU       = 6,
  seed       = 1 
  ){

  # Function to generate systems for the simulation of genome dynamics evolution
  
  # Flag      Format       Description
  # Acont     <numeric>    Amount of adenine content
  # Gcont     <numeric>    Amount of guanine content
  # Ccont     <numeric>    Amount of cysotine content
  # span      <numeric>    Maximum time period for the time evolution
  #                        of genomic base content
  # step      <numeric>    Time step in evaluating the kinetic
  #                        mutation rate equations
  # max.runs  <numeric>    The number of systems to generate in 
  #                        the simulation. Value needs to be above zero
  # muttype   <character>  Mutation rate constants with symmetry constrained
  #                        ("Strand_Symmetric") or without symmetry constrained
  #                        ("Non_Symmetric")
  # dist      <character>  Type of distribution from which the mutation
  #                        rate constants are randomly drawn. Choose from
  #                        "uniform" or "normal" distributions
  # species   <character>  If the simulation is for the time periods 
  #                        to reach chargaff compliance and genome equilibration,
  #                        select species that were reported in the paper from 
  #                        Michael Lynch. Choose from: "H.sapiens", "D.melanogaster"
  #                        "C.elegans", "A.thaliana", "S.cerevisiae", "E.coli". If 
  #                        general simulation, use default parameter of: "NONE"
  # scale.fac  <numeric>   Scaling factor for the standard deviation of the 
  #                        random drawing of mutation rate constants
  # tolerance  <character> Chargaff tolerance for the simulation. 
  #                        If general simulation, use default parameter of FALSE. If
  #                        simulation for time periods to reach chargaff compliance
  #                        and genome equilibration, use TRUE boolean
  # tol.return <character> Return intermediate results for the chargaff tolerance. 
  #                        If general simulation, use default parameter of "NONE". If
  #                        simulation for time periods to reach chargaff compliance
  #                        and genome equilibration, use "equil_time" to obtain
  #                        mean and standard deviation fluctuation of the 
  #                        absolute difference in the base content. Use "fluctuation"
  #                        when results for "equil_time" are obtained
  # sim.evol   <character> Capture the 1 billion year time evolution of the change
  #                        of base content and meta-data. 
  # NCPU       <numeric>   Number of cores to use for parallel execution. 
  # seed       <numeric>   Random number generator, useful for reproducible random objects

  # check input classes
  inputchecking(Acont, Gcont, Ccont, span, step, max.runs, 
                muttype, dist, species, scale.fac, tolerance, 
                tol.return, sim.evol, NCPU, seed)
  
  # obtain states 
  if(tolerance){
    # obtain states 
    all.states <- states(iterations = max.runs, seed = seed, NCPU = NCPU)
  } else {
    Tcont <- 1 - Acont - Gcont - Ccont
    state <- c(Ca=Acont, Cg=Gcont, Ct=Tcont, Cc=Ccont) 
  }
  
  # set-up cluster for parallel computation
  cl <- makeCluster(NCPU)
  registerDoParallel(cl)
  cat(paste("Running with ", NCPU, " cores...", sep="", "\n"))
  
  # Declare that parallel RNG should be used for in a parallel foreach() call.
  # %dorng% will still result in parallel processing; uses %dopar% internally.
  registerDoRNG(seed = seed)

  # checks before simulation commences
  if((muttype == "Strand_Symmetric") & (dist == "uniform")){
    cat("Symmetry-constrained mutation rates only works with normal distribution. Automatically changing...", "\n")
    dist = "normal"
    cat("Running strand symmetric model with normal distribution...", "\n")
  }
  
  sim.run <- foreach(i = 1:max.runs, 
                     .combine="rbind",
                     .export=ls(globalenv()),
                     .packages=c("foreach", "truncnorm", "dplyr"),
                     .inorder=FALSE)%dopar%{
    if(tolerance){
      state <- all.states[i,]
    }

    if((muttype == "Non_Symmetric") & (dist == "uniform")){
     rates <- runif(n=12, min=0, max=note.three$MAXMEAN+note.three$MAXSD*scale.fac)
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
    
    if((muttype == "Non_Symmetric") & (dist == "normal")){
     kag = rtruncnorm(n=1, mean=note.one[note.one$MUT=="AG","MEAN"], 
                       sd=note.one[note.one$MUT=="AG","SD"]*scale.fac, a=0, b=Inf)
     
     kat = rtruncnorm(n=1, mean=note.one[note.one$MUT=="AT","MEAN"], 
                       sd=note.one[note.one$MUT=="AT","SD"]*scale.fac, a=0, b=Inf)
     
     kac = rtruncnorm(n=1, mean=note.one[note.one$MUT=="AC","MEAN"], 
                       sd=note.one[note.one$MUT=="AC","SD"]*scale.fac, a=0, b=Inf)
     
     kga = rtruncnorm(n=1, mean=note.one[note.one$MUT=="GA","MEAN"], 
                       sd=note.one[note.one$MUT=="GA","SD"]*scale.fac, a=0, b=Inf)
     
     kgt = rtruncnorm(n=1, mean=note.one[note.one$MUT=="GT","MEAN"], 
                       sd=note.one[note.one$MUT=="GT","SD"]*scale.fac, a=0, b=Inf)
     
     kgc = rtruncnorm(n=1, mean=note.one[note.one$MUT=="GC","MEAN"], 
                       sd=note.one[note.one$MUT=="GC","SD"]*scale.fac, a=0, b=Inf)
     
     kta = rtruncnorm(n=1, mean=note.one[note.one$MUT=="TA","MEAN"], 
                       sd=note.one[note.one$MUT=="TA","SD"]*scale.fac, a=0, b=Inf)
     
     ktg = rtruncnorm(n=1, mean=note.one[note.one$MUT=="TG","MEAN"], 
                       sd=note.one[note.one$MUT=="TG","SD"]*scale.fac, a=0, b=Inf)
     
     ktc = rtruncnorm(n=1, mean=note.one[note.one$MUT=="TC","MEAN"], 
                       sd=note.one[note.one$MUT=="TC","SD"]*scale.fac, a=0, b=Inf)
     
     kca = rtruncnorm(n=1, mean=note.one[note.one$MUT=="CA","MEAN"], 
                       sd=note.one[note.one$MUT=="CA","SD"]*scale.fac, a=0, b=Inf)
     
     kcg = rtruncnorm(n=1, mean=note.one[note.one$MUT=="CG","MEAN"], 
                       sd=note.one[note.one$MUT=="CG","SD"]*scale.fac, a=0, b=Inf)
     
     kct = rtruncnorm(n=1, mean=note.one[note.one$MUT=="CT","MEAN"], 
                       sd=note.one[note.one$MUT=="CT","SD"]*scale.fac, a=0, b=Inf)
    }
    
    if((muttype == "Strand_Symmetric") & (dist == "normal")){
      if(species == "NONE"){
        kag = ktc = rtruncnorm(n=1, mean=note.two[note.two$MUT=="AG|TC","MEAN"], 
                                sd=note.two[note.two$MUT=="AG|TC","SD"]*scale.fac, a=0, b=Inf)
        
        kat = kta = rtruncnorm(n=1, mean=note.two[note.two$MUT=="AT|TA","MEAN"], 
                                sd=note.two[note.two$MUT=="AT|TA","SD"]*scale.fac, a=0, b=Inf)
        
        kac = ktg = rtruncnorm(n=1, mean=note.two[note.two$MUT=="AC|TG","MEAN"], 
                                sd=note.two[note.two$MUT=="AC|TG","SD"]*scale.fac, a=0, b=Inf)
        
        kct = kga = rtruncnorm(n=1, mean=note.two[note.two$MUT=="CT|GA","MEAN"], 
                                sd=note.two[note.two$MUT=="CT|GA","SD"]*scale.fac, a=0, b=Inf)
        
        kca = kgt = rtruncnorm(n=1, mean=note.two[note.two$MUT=="CA|GT","MEAN"], 
                                sd=note.two[note.two$MUT=="CA|GT","SD"]*scale.fac, a=0, b=Inf)
        
        kcg = kgc = rtruncnorm(n=1, mean=note.two[note.two$MUT=="CG|GC","MEAN"], 
                                sd=note.two[note.two$MUT=="CG|GC","SD"]*scale.fac, a=0, b=Inf)
      } else {
        kag = ktc = rtruncnorm(n=1, mean=note.two[note.two$MUT=="AG|TC",species], 
                               sd=note.two[note.two$MUT=="AG|TC","SD"]*scale.fac, a=0, b=Inf)
   
        kat = kta = rtruncnorm(n=1, mean=note.two[note.two$MUT=="AT|TA",species], 
                               sd=note.two[note.two$MUT=="AT|TA","SD"]*scale.fac, a=0, b=Inf)
        
        kac = ktg = rtruncnorm(n=1, mean=note.two[note.two$MUT=="AC|TG",species], 
                               sd=note.two[note.two$MUT=="AC|TG","SD"]*scale.fac, a=0, b=Inf)
        
        kct = kga = rtruncnorm(n=1, mean=note.two[note.two$MUT=="CT|GA",species], 
                               sd=note.two[note.two$MUT=="CT|GA","SD"]*scale.fac, a=0, b=Inf)
        
        kca = kgt = rtruncnorm(n=1, mean=note.two[note.two$MUT=="CA|GT",species], 
                               sd=note.two[note.two$MUT=="CA|GT","SD"]*scale.fac, a=0, b=Inf)
        
        kcg = kgc = rtruncnorm(n=1, mean=note.two[note.two$MUT=="CG|GC",species], 
                               sd=note.two[note.two$MUT=="CG|GC","SD"]*scale.fac, a=0, b=Inf)
      }
    }   
    
    parameters <- c(kag=kag, kat=kat, kac=kac,
                    kga=kga, kgt=kgt, kgc=kgc,
                    kta=kta, ktg=ktg, ktc=ktc,
                    kca=kca, kcg=kcg, kct=kct) #mut/byr
    
    if(tolerance){
      if (tol.return == "equil_time"){
        atgc <- SolveATGC(parameters = parameters, state = state, 
                          step = step, span = span,
                          EQtolerance = EQtolerance, CHtolerance = CHtolerance)
      } else if (tol.return == "fluctuation") {
        atgc <- SolveATGC(parameters = parameters, state = state, 
                          step = step, span = span,
                          EQtolerance = FALSE, CHtolerance = FALSE)
      }
    } else {
      atgc <- SolveATGC(parameters = parameters, state = state, 
                        step = step, span = span,
                        EQtolerance = FALSE, CHtolerance = FALSE)
    }
    
    if(tolerance){
      if (tol.return == "equil_time"){
        return(data.frame(Chargaff   = atgc$Ch.time,
                          Equil      = atgc$Eq.time,
                          Difference = atgc$Eq.Ch.time.diff))
      } else if (tol.return == "fluctuation") {
        return(data.frame(Mean = atgc$fluctuation.mean,
                          Sd   = atgc$fluctuation.sd))
      }
    } else {
      if(sim.evol){
        return(data.frame(AT=atgc$at.content[atgc$length.out],
                          GC=atgc$gc.content[atgc$length.out],
                          ATratio=atgc$at.ratio[atgc$length.out],
                          GCratio=atgc$gc.ratio[atgc$length.out],
                          ATskew=atgc$at.skew[atgc$length.out],
                          GCskew=atgc$gc.skew[atgc$length.out],
                          one.bn.ATskew=atgc$one.bn$at.skew,
                          one.bn.GCskew=atgc$one.bn$gc.skew,
                          one.bn.AT=atgc$one.bn$at.content,
                          one.bn.GC=atgc$one.bn$gc.content,
                          two.bn.ATskew=atgc$two.bn$at.skew,
                          two.bn.GCskew=atgc$two.bn$gc.skew,
                          two.bn.AT=atgc$two.bn$at.content,
                          two.bn.GC=atgc$two.bn$gc.content,
                          three.bn.ATskew=atgc$three.bn$at.skew,
                          three.bn.GCskew=atgc$three.bn$gc.skew,
                          three.bn.AT=atgc$three.bn$at.content,
                          three.bn.GC=atgc$three.bn$gc.content,
                          four.bn.ATskew=atgc$four.bn$at.skew,
                          four.bn.GCskew=atgc$four.bn$gc.skew,
                          four.bn.AT=atgc$four.bn$at.content,
                          four.bn.GC=atgc$four.bn$gc.content,
                          kag=kag, kat=kat, kac=kac,
                          kga=kga, kgt=kgt, kgc=kgc,
                          kta=kta, ktg=ktg, ktc=ktc,
                          kca=kca, kcg=kcg, kct=kct))  
      } else { 
        return(data.frame(AT=atgc$at.content[atgc$length.out],
                          GC=atgc$gc.content[atgc$length.out],
                          ATratio=atgc$at.ratio[atgc$length.out],
                          GCratio=atgc$gc.ratio[atgc$length.out],
                          ATskew=atgc$at.skew[atgc$length.out],
                          GCskew=atgc$gc.skew[atgc$length.out],
                          kag=kag, kat=kat, kac=kac,
                          kga=kga, kgt=kgt, kgc=kgc,
                          kta=kta, ktg=ktg, ktc=ktc,
                          kca=kca, kcg=kcg, kct=kct))
      }
    }
  } %>% 
  suppressWarnings() # ignore 'already exporting variables' warning from foreach 
  ### end of foreach ######################################
  if(tolerance == FALSE){
    save(sim.run, file=paste0("../../data/Main_Simulation/",muttype,"-",dist,"/",
                              muttype,"-",dist,"-scaling-",scale.fac,".Rdata"))
  } else {
    return(sim.run)
  }
  stopImplicitCluster()
  stopCluster(cl)
}
################################################################################