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

cat("Obtaining average fluctuation tolerance data...", "\n")
sim.results <- Simulation(
  Acont      = 0.25, # %
  Gcont      = 0.25, # %
  Ccont      = 0.25, # %
  span       = 10, # byr
  step       = 0.001, # byr
  max.runs   = 100000, # number of iterations
  muttype    = "Strand_Symmetric",
  dist       = "normal",
  species    = "NONE",
  scale.fac  = 1,
  tolerance  = TRUE,
  tol.return = "fluctuation",
  sim.evol   = FALSE,
  NCPU       = ncpu, 
  seed       = 1
)

saveRDS(sim.results, 
file = "../../data/Chargaff_Equilibrium/ChargaffEquilibrium.Rdata")
cat("Done!", "\n")

Fluc.tol <- readRDS(file = "../../data/Chargaff_Equilibrium/ChargaffEquilibrium.Rdata")
EQtolerance <- mean(Fluc.tol$Mean)*((1/100)*25)

# Run 1M simulations to obtain Chargaff tolerance and genome equilibration
scaling.factor <- c(0,1,2,5,10)
scaling.factor.name <- c("zero", "one", "two", "five", "ten")
for(i in 1:length(scaling.factor)){
  cat("Obtaining Chargaff tolerance and genome equilibration for scaling", scaling.factor[i],"...", "\n")
  sim.results <- Simulation(
    Acont      = 0.25, # %
    Gcont      = 0.25, # %
    Ccont      = 0.25, # %
    span       = 10, # byr
    step       = 0.001, # byr
    max.runs   = 1000000, # number of iterations
    muttype    = "Strand_Symmetric",
    dist       = "normal",
    species    = "NONE",
    scale.fac  = scaling.factor[i],
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

  saveRDS(sim.results, 
  file = paste0("../../data/Chargaff_Equilibrium/ChargaffEquilibriumDistribution_scaling_", 
  scaling.factor.name[i],".Rdata"))
}