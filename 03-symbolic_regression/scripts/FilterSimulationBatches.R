args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
setwd(my_path)

suppressPackageStartupMessages(library(dplyr))
source("../../02-simulations/lib/LoadData.R")

filter.df <- function(dataset, species = "eukaryotes"){

  # Apply the PR-2 tolerance to the compliance region of the datasets

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # species  <character> Apply the PR-2 tolerance zone of 
  #                      "prokaryotes", "eukaryotes" or "viruses" organisms

  if(species == "prokaryotes"){
    file.name <-  "../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv"
  } else if(species == "eukaryotes"){
    file.name <- "../../01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv"
  } else if(species == "viruses"){
    file.name <- "../../01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv"
  }
  
  if(file.exists(file.name)){
    # import chargaff compliancy region
    fluc.species <- read.csv(file = file.name, header = TRUE)
    
    # categorise data as compliant and non-compliant cases
    y <- dataset %>%
      as_tibble() %>%
      mutate(
        Label = ifelse(((ATskew >= fluc.species[fluc.species$metadata == "AT_skew","mean"]-
                           fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                          (ATskew <= fluc.species[fluc.species$metadata == "AT_skew","mean"]+
                             fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                          (GCskew >= fluc.species[fluc.species$metadata == "GC_skew","mean"]-
                             fluc.species[fluc.species$metadata == "GC_skew","st.dev"]) & 
                          (GCskew <= fluc.species[fluc.species$metadata == "GC_skew","mean"]+
                             fluc.species[fluc.species$metadata == "GC_skew","st.dev"])), 
                       1, 0)) %>%
      relocate(Label, .before = kag) %>%
      mutate(Label = as.factor(as.numeric(Label)),
             Label = forcats::fct_recode(Label, 
                                "NO" = "0", 
                                "YES" = "1"))
    
    # select only columns of Label and all rate constants
    sim.run.rates <- y %>%
      select(7:(dim(dataset)[2]-5)) %>%
      as_tibble() %>%
      na.omit()
    
    # filter out all compliant cases
    df.compliant <- sim.run.rates %>%
      filter(Label == "YES")
    
    # return compliant cases
    return(list(sim.run.rates, df.compliant))
  } else {
    stop("File does not exist!")
  }
}

sim.run <- LoadData(file.path = "../data/Simulation/simulation_batches.Rdata",
                    scaling = 1)

df <- filter.df(dataset = sim.run, species = "eukaryotes")

df[[2]] %>%
  write.csv(file = "../data/Test/simulation_batches_compliant.csv",
            row.names = FALSE)