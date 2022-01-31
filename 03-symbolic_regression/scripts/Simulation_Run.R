args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[2])
runs <- as.numeric(args[3])
setwd(my_path)

# Load required supplementary functions and packages
suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(truncnorm))
suppressPackageStartupMessages(library(doRNG))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

# source dependencies
source("../../02-simulations/lib/SolveATGC.R")
source("../../02-simulations/lib/Simulation.R")
source("../../02-simulations/lib/CheckInput.R")
source("../../02-simulations/lib/States.R")

# Import calculated mutation rates from trek paper
note.one   <- read.csv("../../02-simulations/data/Raw/Trek-paper-Note-1-mutation-rates.csv", 
                       header = TRUE)
note.two   <- read.csv("../../02-simulations/data/Raw/Trek-paper-Note-2-mutation-rates.csv", 
                       header = TRUE)
note.three <- read.csv("../../02-simulations/data/Raw/Trek-paper-Note-3-mutation-rates.csv", 
                       header = TRUE)

# run simulation to obtain test set for Eureqa evaluation...
cat("Running simulation to obtain test set for Eureqa evaluation...", "\n")
sim.results <- Simulation(
    Acont      = 0.25, # %
    Gcont      = 0.25, # %
    Ccont      = 0.25, # %
    span       = 4.28, # byr
    step       = 0.001, # byr
    max.runs   = runs, # number of iterations
    muttype    = "Non_Symmetric",
    dist       = "uniform",
    species    = "NONE",
    scale.fac  = 1,
    tolerance  = FALSE,
    tol.return = "NONE",
    sim.evol   = FALSE,
    sy.reg.run = TRUE,
    NCPU       = ncpu, 
    seed       = 2021
)

# save output
saveRDS(sim.results, file = "../data/Simulation/simulation_batches.Rdata")
cat("Done!","\n")