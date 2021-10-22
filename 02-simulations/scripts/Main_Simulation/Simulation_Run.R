# Load required supplementary functions and packages
suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(truncnorm))
suppressPackageStartupMessages(library(doRNG))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[2])
runs <- as.numeric(args[3])
muttype <- as.character(args[4])
dist <- as.character(args[5])
scale.fac <- as.numeric(args[6])
sim.evol <- as.character(args[7])
sim.evol <- as.logical(sim.evol)
setwd(my_path)

# source dependencies
source("../../lib/SolveATGC.R")
source("../../lib/InputChecking.R")
source("../../lib/Simulation.R")

# Import calculated mutation rates from trek paper
note.one   <- read.csv("../../data/Raw/Trek-paper-Note-1-mutation-rates.csv", 
                       header = TRUE)
note.two   <- read.csv("../../data/Raw/Trek-paper-Note-2-mutation-rates.csv", 
                       header = TRUE)
note.three <- read.csv("../../data/Raw/Trek-paper-Note-3-mutation-rates.csv", 
                       header = TRUE)

solsym(Acont      = 0.25, # %
       Gcont      = 0.25, # %
       Ccont      = 0.25, # %
       span       = 4.28, # byr
       step       = 0.001, # byr
       max.runs   = runs, # number of iterations
       muttype    = muttype,
       dist       = dist,
       species    = "NONE",
       scale.fac  = scale.fac,
       tolerance  = FALSE,
       tol.return = "NONE",
       sim.evol   = sim.evol,
       NCPU       = ncpu, 
       seed       = 1)