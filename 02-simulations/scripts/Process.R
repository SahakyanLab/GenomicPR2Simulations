# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
runs <- as.numeric(args[2])
ncpu <- as.numeric(args[3])
save_as <- as.character(args[4])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(deSolve)))
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(truncnorm)))
suppressPackageStartupMessages(suppressWarnings(library(doRNG)))
suppressPackageStartupMessages(suppressWarnings(library(doParallel)))
suppressPackageStartupMessages(suppressWarnings(library(foreach)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(gridExtra)))
suppressPackageStartupMessages(suppressWarnings(library(geneplotter)))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))

my.path="/Users/paddy/Documents/DPhil/01-Chargaff/02-simulations/scripts/"
ncpu=1
save_as="png"

setwd(my.path)
source("../lib/Simulation.R")
source("../lib/Plots.R")
source("../lib/Rates.R")

####################################################################
# Equilibration simulation runs
####################################################################
run_sims <- Simulation$new(
    Acont        = 0.25, # %
    Gcont        = 0.25, # %
    Ccont        = 0.25, # %
    span         = 10, # byr
    step         = 0.001, # byr
    max_runs     = 100000, 
    muttype      = "Strand_Symmetric",
    distribution = "normal",
    species      = "NONE",
    scale_fac    = 1,
    tolerance    = TRUE,
    tol_return   = "fluctuation",
    sim_evol     = FALSE,
    NCPU         = ncpu, 
    seed         = 1    
)
sim.results <- run_sims$run_simulation()
dir.create(
    path = "../data/Chargaff_Equilibrium/",
    showWarnings = FALSE,
    recursive = TRUE
)
saveRDS(
  sim.results, 
  file = "../data/Chargaff_Equilibrium/ChargaffEquilibrium.Rdata"
)
EQtolerance <- mean(sim.results$Mean, na.rm = TRUE)*((1/100)*25)

# run 1M simulations to obtain Chargaff tolerance and genome equilibrium
scaling <- c("zero" = 0, "one" = 1, "two" = 2, "five" = 5, "ten" = 10)
for(x in 1:length(scaling)){
    run_sims <- Simulation$new(
        Acont        = 0.25, # %
        Gcont        = 0.25, # %
        Ccont        = 0.25, # %
        span         = 10, # byr
        step         = 0.001, # byr
        max_runs     = 1000000, 
        muttype      = "Strand_Symmetric",
        distribution = "normal",
        species      = "NONE",
        scale_fac    = unname(scaling[x]),
        tolerance    = TRUE,
        tol_return   = "equil_time",
        sim_evol     = FALSE,
        NCPU         = ncpu, 
        seed         = 2021
    )
    sim.results <- run_sims$run_simulation()

    to.remove <- apply(sim.results, 1, function(x){any(is.na(x))})
    if(length(which(to.remove))>0){
        sim.results <- sim.results[!to.remove,]
        rownames(sim.results) <- NULL
    }

    saveRDS(
        sim.results, 
        file = paste0("../data/Chargaff_Equilibrium/", 
                      "ChargaffEquilibriumDistribution_scaling_", 
                      names(scaling[x]), ".Rdata")
    )
}

# generate plots
gen_plots <- Plots$new(equilibration = TRUE)
gen_plots$generate_equil_plots(save_as = save_as)

####################################################################
# Main simulation runs
####################################################################
# independent mutation rate constants
for(distr in c("uniform","normal")){
    for(scaling in c(1,2,5,10)){
        # run simulation
        run_sims <- Simulation$new(
            Acont        = 0.25, # %
            Gcont        = 0.25, # %
            Ccont        = 0.25, # %
            span         = 4.28, # byr
            step         = 0.001, # byr
            max_runs     = runs, 
            muttype      = "Non_Symmetric",
            distribution = distr,
            species      = "NONE",
            scale_fac    = scaling,
            tolerance    = FALSE,
            tol_return   = "NONE",
            sim_evol     = FALSE,
            NCPU         = ncpu, 
            seed         = 1    
        )
        run_sims$run_simulation()

        # generate plots
        gen_plots <- Plots$new(
            muttype       = "Non_Symmetric",
            distribution  = distr,
            scale_fac     = scaling,
            equilibration = FALSE
        )
        gen_plots$generate_sim_plots(save_as = save_as)
    }
}

# strand-symmetric mutation rate constants
for(scaling in c(1,2,5,10)){
    # run simulation
    run_sims <- Simulation$new(
        Acont        = 0.25, # %
        Gcont        = 0.25, # %
        Ccont        = 0.25, # %
        span         = 4.28, # byr
        step         = 0.001, # byr
        max_runs     = runs, 
        muttype      = "Strand_Symmetric",
        distribution = "normal",
        species      = "NONE",
        scale_fac    = scaling,
        tolerance    = FALSE,
        tol_return   = "NONE",
        sim_evol     = TRUE,
        NCPU         = ncpu, 
        seed         = 1    
    )
    run_sims$run_simulation()

    # generate plots
    gen_plots <- Plots$new(
        muttype       = "Strand_Symmetric",
        distribution  = "normal",
        scale_fac     = scaling,
        equilibration = FALSE
    )
    gen_plots$generate_sim_plots(save_as = save_as)
}

####################################################################
# Obtain mutation rate constants from real life species
####################################################################
source("../lib/Rates.R")
get_rates <- Rates$new(trek_scale = FALSE)
get_rates$process_rates(save_as = save_as)

# equilibrium simulation with these rate constants
lynch.rates <- read.csv(
    file = paste0("../data/Raw/Michael_Lynch/Lynch-", 
                  "2010-converted-mutation-rates.csv"),
    header = TRUE
)
species <- names(lynch.rates[,6:length(lynch.rates)])
for(i in 1:length(species)){
    # run simulation
    run_sims <- Simulation$new(
        Acont        = 0.25, # %
        Gcont        = 0.25, # %
        Ccont        = 0.25, # %
        span         = 10, # byr
        step         = 0.001, # byr
        max_runs     = 1000000, 
        muttype      = "Strand_Symmetric",
        distribution = "normal",
        species      = species[i],
        scale_fac    = 0,
        tolerance    = TRUE,
        tol_return   = "equil_time",
        sim_evol     = FALSE,
        NCPU         = ncpu, 
        seed         = 2021    
    )
    sim.results <- run_sims$run_simulation()

    to.remove <- apply(sim.results, 1, function(x){any(is.na(x))})
    if(length(which(to.remove))>0){
        sim.results <- sim.results[!to.remove,]
        rownames(sim.results) <- NULL
    }

    species.name <- gsub(pattern = "\\.", replacement = "_", x = species[i])
    saveRDS(
        sim.results, 
        file = paste0("../data/Chargaff_Equilibrium/", 
                      "ChargaffEquilibriumDistribution_",
                      species.name, ".Rdata")
    )
}

# generate plots
gen_plots <- Plots$new(equilibration = TRUE)
gen_plots$generate_equil_plots(
    save_as = save_as, 
    other_species = TRUE
)

####################################################################
# Simulation for test set for symbolic regression
####################################################################
# run simulation
run_sims <- Simulation$new(
    Acont        = 0.25, # %
    Gcont        = 0.25, # %
    Ccont        = 0.25, # %
    span         = 10, # byr
    step         = 0.001, # byr
    max_runs     = 10000000, 
    muttype      = "Non_Symmetric",
    distribution = "uniform",
    species      = "NONE",
    scale_fac    = 1,
    tolerance    = FALSE,
    tol_return   = "NONE",
    sim_evol     = FALSE,
    sy.reg.run   = TRUE,
    NCPU         = ncpu, 
    seed         = 2021    
)
sim.results <- run_sims$run_simulation()

# save output
dir.create(
    path = "../../04-symbolic_regression/data/",
    showWarnings = FALSE,
    recursive = TRUE
)
saveRDS(
    sim.results, 
    file = "../data/Simulation/simulation_batches.Rdata"
)