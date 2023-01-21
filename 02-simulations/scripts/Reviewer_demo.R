# # read arguments from job submission
# args <- commandArgs(trailingOnly = TRUE)

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

# run simulation where the initial base content is random.
run_sims <- Simulation$new(
    span              = 4.28, # byr
    step              = 0.001, # byr
    max_runs          = 100000,
    muttype           = "Non_Symmetric",
    distribution      = "uniform",
    species           = "NONE",
    scale_fac         = 1,
    tolerance         = FALSE,
    random_init_bases = TRUE,
    tol_return        = "NONE",
    sim_evol          = FALSE,
    NCPU              = ncpu, 
    seed              = 1
)
run_sims$run_simulation()

# generate plots
gen_plots <- Plots$new(
    muttype           = "Non_Symmetric",
    distribution      = "uniform",
    scale_fac         = 1,
    equilibration     = FALSE,
    random_init_bases = TRUE
)
gen_plots$generate_sim_plots(save_as = save_as)

# run simulation where initial genome content is such that 
# each generation is any of the following at skew, gc skew 
# edge values of (1,1), (-1,1), (1,-1), (-1,-1).
run_sims <- Simulation$new(
    span              = 4.28, # byr
    step              = 0.001, # byr
    max_runs          = 100000,
    muttype           = "Non_Symmetric",
    distribution      = "uniform",
    species           = "NONE",
    scale_fac         = 1,
    tolerance         = FALSE,
    edge_skew_cases   = TRUE,
    tol_return        = "NONE",
    sim_evol          = FALSE,
    NCPU              = ncpu, 
    seed              = 1
)
run_sims$run_simulation()

# generate plots
gen_plots <- Plots$new(
    muttype           = "Non_Symmetric",
    distribution      = "uniform",
    scale_fac         = 1,
    equilibration     = FALSE,
    edge_skew_cases   = TRUE
)
gen_plots$generate_sim_plots(save_as = save_as)