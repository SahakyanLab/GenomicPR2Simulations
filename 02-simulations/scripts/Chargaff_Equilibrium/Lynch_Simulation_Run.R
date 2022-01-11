args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[2])
setwd(my_path)

# Load required supplementary functions and packages
suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(truncnorm))
suppressPackageStartupMessages(library(doRNG))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))

source("../../lib/SolveATGC.R")
source("../../lib/Simulation.R")
source("../../lib/States.R")
source("../../lib/CheckInput.R")

prokaryotes.df <- read.csv(file = "../../../01-genome_composition/data/01-Prokaryotes/All/all_filtered_dataframe.csv", 
                           header=TRUE)
eukaryotes.df  <- read.csv(file = "../../../01-genome_composition/data/02-Eukaryotes/All/all_filtered_dataframe.csv", 
                           header=TRUE)

# mutation rates from trek paper
note.one   <- read.csv("../../data/Raw/Trek-paper-Note-1-mutation-rates.csv", 
                       header = TRUE)
note.two   <- read.csv("../../data/Raw/Trek-paper-Note-2-mutation-rates.csv", 
                       header = TRUE)
note.three <- read.csv("../../data/Raw/Trek-paper-Note-3-mutation-rates.csv", 
                       header = TRUE)
CHtolerance  <- read.csv(file = "../../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv", 
                       header = TRUE)

# mutation rates from lynch paper
lynch.rates <- read.csv(file = "../../data/Raw/Michael_Lynch/Lynch-2010-converted-mutation-rates.csv",
                        header = TRUE)

note.two <- note.two %>% 
  mutate(lynch.rates[,5:length(lynch.rates)])

# Load Equilibrium fluctuation tolerance level
Fluc.tol <- readRDS(file = "../../data/Chargaff_Equilibrium/ChargaffEquilibrium.Rdata")
EQtolerance <- mean(Fluc.tol$Mean)*((1/100)*25)

# Run 1M simulations to obtain Chargaff tolerance and genome equilibration
species <- names(lynch.rates[,6:length(lynch.rates)])
for(i in 1:species){
  cat(paste0("Obtaining Chargaff tolerance and genome equilibration for species: ", species[i],"..."))
  sim.results <- Simulation(
    Acont      = 0.25, # %
    Gcont      = 0.25, # %
    Ccont      = 0.25, # %
    span       = 10, # byr
    step       = 0.001, # byr
    max.runs   = 1000000, # number of iterations
    muttype    = "Strand_Symmetric",
    dist       = "normal",
    species    = species[i],
    scale.fac  = 0,
    tolerance  = TRUE,
    tol.return = "equil_time",
    sim.evol   = FALSE,
    NCPU       = ncpu, 
    seed       = 2021
  )

  to.remove <- apply(sim.results, 1, function(x){any(is.na(x))})
  if(length(which(to.remove))>0){
    sim.results <- sim.results[!to.remove,]
    rownames(sim.results) <- NULL
  }

  species.name <- gsub(pattern = "\\.", replacement = "_", x = species[i])
  saveRDS(sim.results, 
          file = paste0("../../data/Chargaff_Equilibrium/ChargaffEquilibriumDistribution_",
                        species.name, ".Rdata"))
}