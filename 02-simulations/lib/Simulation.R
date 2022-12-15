Simulation <- R6::R6Class(
    classname = "Simulation",
    public = list(
        initialize = function(Acont, Gcont, Ccont, span, step, max_runs,
                              muttype, distribution, species, scale_fac,
                              tolerance, tol_return, sim_evol,
                              sy_reg_run, NCPU, seed){
            if(!missing(Acont)) private$Acont <- Acont
            if(!missing(Gcont)) private$Gcont <- Gcont
            if(!missing(Ccont)) private$Ccont <- Ccont
            if(!missing(step)) private$step <- step
            if(!missing(span)) private$span <- span
            if(!missing(max_runs)) private$max_runs <- max_runs
            if(!missing(muttype)) private$muttype <- muttype
            if(!missing(distribution)) private$distribution <- distribution
            if(!missing(species)) private$species <- species 
            if(!missing(scale_fac)) private$scale_fac <- scale_fac 
            if(!missing(tolerance)) private$tolerance <- tolerance 
            if(!missing(tol_return)) private$tol_return <- tol_return
            if(!missing(sim_evol)) private$sim_evol <- sim_evol
            if(!missing(sy_reg_run)) private$sy_reg_run <- sy_reg_run
            if(!missing(NCPU)) private$NCPU <- NCPU 
            if(!missing(seed)) private$seed <- seed

            # checks
            if((!Acont<1 | !Gcont<1 | !Ccont<1))
                stop("Base content need to be numeric values below one.")
            if(sum(Acont, Gcont, Ccont)>=1)
                stop("Base contents need to sum to less than one so Tcont is non-zero.")
            if(max_runs%%1 != 0 | max_runs<=0)
                stop("max_runs needs to be an integer of positive value.")
            if(!is.character(muttype))
                stop("Muttype must be of character class.")
            if(!is.character(distribution))
                stop("Muttype must be of character class.")
            if(!is.character(species))
                stop("Muttype must be of character class.")
            if(scale_fac%%1 != 0 | scale_fac<=0)
                stop("scale_fac needs to be of positive integer value.")
            if(!scale_fac %in% c(0,1,2,5,10))
                stop("To reproduce the results in the paper, please use scale_fac = c(0,1,2,5,10).")
            if(!is.logical(tolerance))
                stop("Tolerance can only be TRUE or FALSE") 
            if(!is.character(tol_return))
                stop("tol_return must be of character class.")
            if(!is.logical(sim_evol))
                stop("sim_evol can only be TRUE or FALSE") 
            if(NCPU%%1 != 0 | NCPU<=0)
                stop("NCPU must be a positive integer values for parallel execution.")
            if(seed%%1 != 0 | seed<=0)
                stop("seed must be a positive integer values for producing random objects.")

            # get prokaryote and eukaryote genomic information
            private$get_species_df()

            # get pre-calculated mutation rates from trek paper
            private$get_mut_rates()
        },

        #' @description
        #' Generate systems by solving kinetic equations.
        #' @return None or Data.frame.
        run_simulation = function(){
            private$to_plot <- FALSE
            private$generate_systems()
        },

        #' @description
        #' Run one iteration of the genome dynamics simulation and generate plot.
        #' @param time.unit Character vector of the time unit used in the simulation.
        #' @param xlim Numeric vector of the X-axis limits in the plot.
        #' @param ylim Numeric vector of the Y-axis limits in the plot.
        #' @return None.
        get_example_plot = function(time.unit = "byr", xlim = c(0,5), ylim = c(0,60)){
            private$to_plot <- TRUE
            original.sim.runs <- private$max_runs
            private$max_runs <- 1
            sim.output <- private$generate_systems()
            private$plotATGC(atgc = sim.output, time.unit = time.unit, xlim = xlim, ylim = ylim)
            private$max_runs <- original.sim.runs
        }
    ),
    private = list(
        #' @field Acont Numeric vector of the amount of adenine starting content.
        Acont = 0.25,

        #' @field Gcont Numeric vector of the amount of guanine starting content.
        Gcont = 0.25,

        #' @field Ccont Numeric vector of the amount of cytosine starting content.
        Ccont = 0.25,

        #' @field Tcont Numeric vector of the amount of thymine starting content.
        Tcont = NULL,

        #' @field step Numeric vector of the time step in evaluating the kinetic equations.
        step = 0.001,

        #' @field span Numeric vector of the simulation end time in billion years.
        span = 4.28,

        #' @field  Numeric vector of the number of systems to generate
        #'  in the simulation. Value needs to be above zero.
        max_runs = 5,

        #' @field muttype Character vector of c("Strand_Symmetric", "Non_Symmetric") for
        #'  with and without symmetry-constrained mutation rate constants, respectively. 
        muttype = "Non_Symmetric",

        #' @field distribution Character vector of c("uniform", "normal"). 
        #'  The type of distribution from which the mutation rate constants to randomly draw from. 
        distribution = "uniform",

        #' @field species Character vector. If the simulation is for the time periods 
        #'  to reach chargaff compliance and genome equilibration, select species that
        #'  are reported in the paper from Michael Lynch. Choose from: 
        #'  c("H.sapiens", "D.melanogaster", "C.elegans", "A.thaliana", "S.cerevisiae", "E.coli"). 
        #'  If running a general simulation, use default parameter of: "NONE"
        species = "NONE", 

        #' @field lynch_rates Data.frame of mutation rate constants from 
        #'  real life species (Lynch paper).
        lynch_rates = NULL,
        
        #' @field Numeric vector of c(1,2,5,10), indicating the standard deviation of the
        #'  random drawing of mutation rate constants.
        scale_fac = 1,
        
        #' @field Boolean of the Charaff tolerance for the simulation. 
        #'  If FALSE, general simulation run (default). 
        #'  If TRUE, simulation for time periods to reach chargaff compliance
        #'  and genome equilibration.
        tolerance = FALSE,
        
        #' @field Character vector of c("NONE", "equil_time", "fluctuation"). 
        #'  Return intermediate results for the chargaff tolerance. 
        #'  Use "NONE" for a general simulation run. 
        #'  Use "equil_time" for simlation for time periods to reach 
        #'  chargaff compliance and equilibration and you will obtain the mean and st.dev fluctuation 
        #'  of the absolute difference in the base content. 
        #'  Use "fluctuation" when results for "equil_time" are obtained.
        tol_return = "NONE",
        
        #' @field sim_evol Boolean. If TRUE, capture the 1 billion year time evolution of the 
        #'  change of base content and meta-data.
        sim_evol = FALSE, 
        
        #' @field sy_reg_run Boolean. If TRUE, run simulation for obtaining 
        #'  the test results for the symbolic regression run. 
        sy_reg_run = FALSE,

        #' @field EQtolerance Numeric vector of the tolerance value for reaching genome equilibration
        EQtolerance = NULL,

        #' @field CHtolerance Numeric vector of the toelrance value for reahcing Chargaff compliance.
        CHtolerance = NULL,
        
        #' @field NCPU Numeric vector of the number of cores to use 
        #'  for parallel execution.
        NCPU = 4,

        #' @field seed Numeric vector to set the seed for the random number 
        #'  generator, useful for reproducible random objects.
        seed = 1,

        #' @field to_plot Boolean. If TRUE, will run the genome dynamics simulation with one
        #'  iteration and plot that result.
        to_plot = FALSE,

        #' @field prokaryote_df Data.frame of genomic information of all prokaryotes.
        prokaryote_df = NULL,

        #' @field eukaryote_df Data.frame of genomic information of all eukaryotes.
        eukaryote_df = NULL,

        #' @field notes Data.frame of mutation rate constants from trek paper.
        notes = NULL,

        #' @description
        #' Import genomic information for all species of the prokaryote and eukaryote kingdoms.
        #' @return None.
        get_species_df = function(){
            private$prokaryote_df <- read.csv(
                file = paste0("../../01-genome_composition/data/",
                              "01-Prokaryotes/All/all_filtered_dataframe.csv"),
                header = TRUE
            )

            private$eukaryote_df <- read.csv(
                file = paste0("../../01-genome_composition/data/",
                              "02-Eukaryotes/All/all_filtered_dataframe.csv"),
                header = TRUE
            )
        },

        #' @description
        #' Import the pre-calculated mutation rate constants from trek paper.
        #' @return None.
        get_mut_rates = function(){
            private$notes$note_one <- read.csv(
                file = "../data/Raw/Trek-paper-Note-1-mutation-rates.csv",
                header = TRUE
            )

            private$notes$note_two <- read.csv(
                file = "../data/Raw/Trek-paper-Note-2-mutation-rates.csv",
                header = TRUE
            )
            if(private$species != "NONE"){
                private$lynch_rates <- read.csv(
                    file = paste0("../data/Raw/Michael_Lynch/", 
                                  "Lynch-2010-converted-mutation-rates.csv"),
                    header = TRUE
                )
                private$notes$note_two <- private$notes$note_two %>% 
                    dplyr::mutate(lynch.rates[,5:length(lynch.rates)])
            }

            private$notes$note_three <- read.csv(
                file = "../data/Raw/Trek-paper-Note-3-mutation-rates.csv",
                header = TRUE
            )

            private$CHtolerance <- read.csv(
                file = paste0("../../01-genome_composition/data/01-Prokaryotes/", 
                              "PR2_compliance/PR2_fluctuations.csv"),
                header = TRUE
            )
        },

        #' @description
        #' Generate systems for the simulation of genome dynamics evolution
        #' @return None.
        generate_systems = function(){  
            start.time <- Sys.time()
            cur.msg <- "Running evolutionary genome dynamics simulations"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            # time taken for full processing for this experiment
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cur.msg <- paste0("\n",
                "Max. time period (byr):  ", private$span, "\n",
                "Iterations:              ", private$max_runs, "\n",
                "Scaling factor:          ", private$scale_fac, "\n",
                "Distribution:            ", private$distribution, "\n",
                "Mutation type:           ", private$muttype, "\n",
                "CPUs:                    ", private$NCPU
            )
            cat("Quick infosheet:", "\n", cur.msg, "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))

            # obtain states 
            if(private$tolerance){
                all.states <- private$get_states()
            } else {
                private$Tcont <- 1-private$Acont-private$Gcont-private$Ccont
                state <- c(
                    Ca=private$Acont, 
                    Cg=private$Gcont, 
                    Ct=private$Tcont, 
                    Cc=private$Ccont
                ) 
            }

            if(private$NCPU > 1){
                `%op%` <- `%dopar%`
            } else {
                `%op%` <- `%do%`
            }            

            # set-up cluster for parallel computation
            cl <- makeCluster(private$NCPU)
            registerDoParallel(cl)

            # Declare that parallel RNG should be used for in a parallel foreach() call.
            # %dorng% will still result in parallel processing; uses %dopar% internally.
            registerDoRNG(seed = private$seed)

            # checks before simulation commences
            if((private$muttype == "Strand_Symmetric") & (private$distribution == "uniform")){
                # Symmetry-constrained mutation rates only works with normal distribution.
                # Hence, running strand symmetric model with normal distribution.
                private$distribution = "normal"
            }

            sim.run <- foreach(i = 1:private$max_runs, 
                                .combine="rbind",
                                .export=c(ls(globalenv()), "private", "self"),
                                .packages=c("foreach", "truncnorm", "dplyr", "deSolve", "R6"),
                                .inorder=FALSE)%op%{
                Simulation$parent_env <- environment()
                if(private$tolerance){
                    state <- all.states[i,]
                }

                if((private$muttype == "Non_Symmetric") & (private$distribution == "uniform")){
                    rates <- runif(
                        n=12, 
                        min=0, 
                        max=private$notes$note_three$MAXMEAN+private$notes$note_three$MAXSD*private$scale_fac
                    )
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
                
                if((private$muttype == "Non_Symmetric") & (private$distribution == "normal")){
                    kag = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="AG","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="AG","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kat = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="AT","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="AT","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kac = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="AC","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="AC","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kga = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="GA","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="GA","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kgt = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="GT","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="GT","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kgc = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="GC","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="GC","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kta = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="TA","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="TA","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    ktg = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="TG","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="TG","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    ktc = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="TC","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="TC","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kca = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="CA","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="CA","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kcg = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="CG","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="CG","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                    
                    kct = rtruncnorm(
                        n=1, 
                        mean=private$notes$note_one[private$notes$note_one$MUT=="CT","MEAN"], 
                        sd=private$notes$note_one[private$notes$note_one$MUT=="CT","SD"]*private$scale_fac, 
                        a=0, b=Inf
                    )
                }
                
                if((private$muttype == "Strand_Symmetric") & (private$distribution == "normal")){
                    if(private$species == "NONE"){
                        kag = ktc = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="AG|TC","MEAN"], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="AG|TC","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kat = kta = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="AT|TA","MEAN"], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="AT|TA","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kac = ktg = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="AC|TG","MEAN"], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="AC|TG","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kct = kga = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="CT|GA","MEAN"], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="CT|GA","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kca = kgt = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="CA|GT","MEAN"], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="CA|GT","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kcg = kgc = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="CG|GC","MEAN"], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="CG|GC","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                    } else {
                        kag = ktc = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="AG|TC",private$species], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="AG|TC","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                
                        kat = kta = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="AT|TA",private$species], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="AT|TA","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kac = ktg = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="AC|TG",private$species], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="AC|TG","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kct = kga = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="CT|GA",private$species], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="CT|GA","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kca = kgt = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="CA|GT",private$species], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="CA|GT","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                        
                        kcg = kgc = rtruncnorm(
                            n=1, 
                            mean=private$notes$note_two[private$notes$note_two$MUT=="CG|GC",private$species], 
                            sd=private$notes$note_two[private$notes$note_two$MUT=="CG|GC","SD"]*private$scale_fac, 
                            a=0, b=Inf
                        )
                    }
                }   
                
                parameters <- c(kag=kag, kat=kat, kac=kac,
                                kga=kga, kgt=kgt, kgc=kgc,
                                kta=kta, ktg=ktg, ktc=ktc,
                                kca=kca, kcg=kcg, kct=kct) #mut/byr
                
                if(private$tolerance){
                    if (private$tol_return == "equil_time"){
                        atgc <- private$solveATGC(parameters = parameters, state = state)
                    } else if (tol_return == "fluctuation") {
                        private$EQtolerance = private$CHtolerance = FALSE
                        atgc <- private$solveATGC(parameters = parameters, state = state)
                    }
                } else {
                    private$EQtolerance = private$CHtolerance = FALSE
                    atgc <- private$solveATGC(parameters = parameters, state = state)
                }

                if(private$to_plot) return(atgc)

                if(private$tolerance){
                    if (private$tol_return == "equil_time"){
                        return(data.frame(
                            Chargaff   = atgc$Ch.time,
                            Equil      = atgc$Eq.time,
                            Difference = atgc$Eq.Ch.time.diff
                        ))
                    } else if(private$tol_return == "fluctuation") {
                        return(data.frame(
                            Mean = atgc$fluctuation.mean,
                            Sd   = atgc$fluctuation.sd
                        ))
                    }
                } else {
                    if(private$sim_evol){
                        return(data.frame(
                            AT=atgc$at.content[atgc$length.out],
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
                            kca=kca, kcg=kcg, kct=kct
                        ))
                    } else { 
                        return(data.frame(
                            AT=atgc$at.content[atgc$length.out],
                            GC=atgc$gc.content[atgc$length.out],
                            ATratio=atgc$at.ratio[atgc$length.out],
                            GCratio=atgc$gc.ratio[atgc$length.out],
                            ATskew=atgc$at.skew[atgc$length.out],
                            GCskew=atgc$gc.skew[atgc$length.out],
                            kag=kag, kat=kat, kac=kac,
                            kga=kga, kgt=kgt, kgc=kgc,
                            kta=kta, ktg=ktg, ktc=ktc,
                            kca=kca, kcg=kcg, kct=kct
                        ))
                    }
                }
            } %>% 
            suppressWarnings() # ignore 'already exporting variables' warning from foreach

            # return or save results
            if(private$tolerance | private$to_plot){
                return(sim.run)
            } else {
                if(private$sy_reg_run){
                return(sim.run)
                } else {
                    dir.create(
                        path = paste0("../data/Main_Simulation/",private$muttype,
                                      "-",private$distribution,"/"),
                        showWarnings = FALSE,
                        recursive = TRUE
                    )
                    saveRDS(
                        sim.run,
                        file=paste0("../data/Main_Simulation/",private$muttype,
                                    "-",private$distribution,"/", private$muttype,"-",
                                    private$distribution,"-scaling-",private$scale_fac,".Rdata")
                    )
                }
            }
            stopImplicitCluster()
            stopCluster(cl)

            # time taken for full processing for this experiment
            final.t <- Sys.time() - start.time
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },

        #' @description
        #' Obtains min and max single nucleotide content from real-life species.
        #' @return None.
        get_states = function(){
            df <- c(private$eukaryote_df$G_content,
                    private$eukaryote_df$C_content,
                    private$eukaryote_df$A_content,
                    private$eukaryote_df$T_content,
                    private$prokaryote_df$G_content,
                    private$prokaryote_df$C_content,
                    private$prokaryote_df$A_content,
                    private$prokaryote_df$T_content)
            
            # randomly select the order of single nucleotide bases
            rng <- set.seed(private$seed)
            
            random.state <- sapply(1:private$max_runs, function(x){
                base <- c(Ca = 0, Cg = 0, Ct = 0, Cc = 0)
                base.order <- sample(x = 4, size = 4, replace = FALSE)
                
                # first base replacement
                base[base.order[1]] <- runif(n = 1, min = range(df)[1], max = range(df)[2])
                
                # second base replacement
                base[base.order[2]] <- runif(n = 1, min = range(df)[1], max = range(df)[2])
                
                # third base replacement
                new.max <- 100-(sum(base)*100)
                base[base.order[3]] <- runif(n = 1, min = range(df)[1], max = new.max/100)
                
                # final base replacement
                remainder <- 100-(sum(base)*100)
                base[base.order[4]] <- remainder/100
                
                return(base)
            })
            return(t(random.state))
        },

        #' @description
        #' Numerically solves the kinetic rate equations for the simulation 
        #' @param parameters Numeric vector of all the mutation rate constants.
        #' @param state Numeric vector of the states of all four base contents. 
        #'  of the genome dynamics evolution
        #' @return List of info on system. 
        solveATGC = function(parameters, state){
            times <- seq(0, private$span, by=private$step) # byr

            #' @description
            #' solves system of ODEs
            #' @param t time sequence for which output is wanted.
            #' @param state the initial state values for the ODE system.
            #' @param parameters Numeric vector of all the mutation rate constants.
            #' @return matrix of class deSolve.
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
            RESULTS$inp$step            <- private$step
            RESULTS$inp$span            <- private$span
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
            if(private$sim_evol){
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
            if(private$tolerance){
                dif <- abs(diff(out[,c("Ca","Cg","Ct","Cc")]))
                dif.max.min <- apply(dif, 1, function(x){max(x) - min(x)})
                Eq.time <- times[which(dif.max.min <= private$EQtolerance)[1]]

                # Chargaff compliance
                at.tol <- private$CHtolerance[private$CHtolerance$metadata == "AT_skew", "st.dev"]
                gc.tol <- private$CHtolerance[private$CHtolerance$metadata == "GC_skew", "st.dev"]
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
                RESULTS$inp$EQtolerance     <- private$EQtolerance
                RESULTS$Eq.time             <- Eq.time           # time to reach the equilibration tolerance limit
                RESULTS$Ch.time             <- Ch.time           # time to reach Chargaff's tolerance limit
                RESULTS$Eq.Ch.time.diff     <- Eq.Ch.time.diff   # time difference in equilibration and chargaff compliance
            }
            return(RESULTS)
        },

        #' @description
        #' Plots the results of the numerically solved kinetic rate equations for the 
        #' simulation of the genome dynamics evolution.
        #' @param atgc List of all the numerically solved kinetic mutation rate
        #'  equations of one generation/iteration.
        #' @param time.unit Character vector of the time unit used in the simulation.
        #' @param xlim Numeric vector of the X-axis limits in the plot.
        #' @param ylim Numeric vector of the Y-axis limits in the plot.
        #' @return ggplot object.
        plotATGC = function(atgc, time.unit, xlim, ylim){
            p1 <- as_tibble(atgc$out) %>%
                dplyr::mutate(across(.cols = c(Ca, Cg, Ct, Cc), ~.*100)) %>%
                ggplot(aes(x = time,
                        y = Ca)) + 
                geom_line(aes(y = Ca), linewidth = 1.2, col = "forestgreen") + 
                geom_line(aes(y = Cg), linewidth = 1.2, col = "orange") + 
                geom_line(aes(y = Ct), linewidth = 1.2, col = "red") + 
                geom_line(aes(y = Cc), linewidth = 1.2, col = "blue")

            if(!is.null(atgc$Ch.time) & !is.null(atgc$Eq.time)){
                p1 <- p1 +
                    geom_vline(xintercept = atgc$Ch.time, linetype = "dashed") + 
                    geom_text(data = data.frame(xpos = atgc$Ch.time, 
                                                ypos =  min(ylim),
                                                annotateText = "Chargaff eq. reached",
                                                hjustvar = 0, vjustvar = 1.1), 
                            aes(x = xpos, 
                                y = ypos, 
                                hjust = hjustvar, 
                                vjust = vjustvar, 
                                label = annotateText,
                                angle = 90),
                            fontface = "bold", linewidth = 3) + 
                    geom_vline(xintercept = atgc$Eq.time) +
                    geom_text(data = data.frame(xpos = atgc$Eq.time, 
                                                ypos =  min(ylim),
                                                annotateText = "Genome eq. reached",
                                                hjustvar = 0, vjustvar = 1.1), 
                            aes(x = xpos, 
                                y = ypos, 
                                hjust = hjustvar, 
                                vjust = vjustvar, 
                                label = annotateText,
                                angle = 90),
                            fontface = "bold", linewidth = 3)
            }

            p1 <- p1 + 
                coord_cartesian(xlim = xlim, ylim = ylim) + 
                labs(x = "Time, byr",
                    y = "Base content, %") + 
                geom_text(data = data.frame(
                xpos = 0,
                ypos = max(ylim),
                annotateText = paste(
                    "t = 0 ",time.unit,"\n",
                    "nG (orange) = ",signif(atgc$inp$state["Cg"], 2),"\n",
                    "nC (blue) = ",signif(atgc$inp$state["Cc"], 2),"\n",
                    "nA (green) = ",signif(atgc$inp$state["Ca"], 2),"\n",
                    "nT (red) = ",signif(atgc$inp$state["Ct"], 2),"\n",
                    "G+C = ",signif(100*(atgc$inp$state["Cg"]+atgc$inp$state["Cc"])/
                                    atgc$length.genome,2),"%","\n",
                    "G/C = ",signif(atgc$inp$state["Cg"]/atgc$inp$state["Cc"],2),"\n",   
                    "A/T = ",signif(atgc$inp$state["Ca"]/atgc$inp$state["Ct"],2), sep=""),
                hjustvar = 0, vjustvar = 1), 
                        aes(x = xpos, 
                            y = ypos, 
                            hjust = hjustvar, 
                            vjust = vjustvar, 
                            label = annotateText,
                            angle = 0)) + 
                geom_text(data = data.frame(
                xpos = max(xlim)-1,
                ypos = max(ylim), 
                annotateText = paste(
                    "t = ",signif(atgc$inp$span,2)," ",time.unit,"\n",
                    "nG = ",signif(as.vector(atgc$out[atgc$length.out,"Cg"]),2),"\n",
                    "nC = ",signif(as.vector(atgc$out[atgc$length.out,"Cc"]),2),"\n",
                    "nA = ",signif(as.vector(atgc$out[atgc$length.out,"Ca"]),2),"\n",
                    "nT = ",signif(as.vector(atgc$out[atgc$length.out,"Ct"]),2),"\n",
                    "G+C = ",signif(atgc$gc.content[atgc$length.out],2),"%","\n",
                    "G/C = ",signif(atgc$out[atgc$length.out,"Cg"]/atgc$out[atgc$length.out,"Cc"],2),"\n",   
                    "A/T = ",signif(atgc$out[atgc$length.out,"Ca"]/atgc$out[atgc$length.out,"Ct"],2), sep=""),
                hjustvar = 0, vjustvar = 1), 
                        aes(x = xpos, 
                            y = ypos, 
                            hjust = hjustvar, 
                            vjust = vjustvar, 
                            label = annotateText,
                            angle = 0)) + 
                theme_bw() %>% 
                suppressWarnings()
            return(p1)
        }
    )
)