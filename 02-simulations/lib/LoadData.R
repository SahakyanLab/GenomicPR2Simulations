################################################################################
LoadData <- function(file.path, scaling = 1){

  # Function to load the dataset as tibble class

  # Flag       Format        Description
  # file.path  <character>   Path to file for the dataset
  # scaling    <numeric>     Scaling factor for the standard deviation of the 
  #                          random drawing of mutation rate constants

  if(!is.character(file.path)){
    stop("file.path needs to be a character vector.")
  }

  if(!scaling %in% c(0,1,2,5,10)){
    stop("To reproduce the results in the paper, please use scale.fac = c(0,1,2,5,10).")
  }
  
  if(file.exists(file.path)){
    cat("Loading scaling ",scaling, "...")

    if(grepl(pattern = ".csv", x = file.path, fixed = TRUE)){
      sim.run <- read.csv(file = file.path, header = TRUE)
    } else if (grepl(pattern = ".Rdata", x = file.path, fixed = TRUE)){
      sim.run <- readRDS(file = file.path)
    }
    
    sim.run$nC <- (sim.run$GC/100)/(1+sim.run$GCratio)
    sim.run$nG <- (sim.run$GC/100)-sim.run$nC
    sim.run$nT <- (1-(sim.run$GC/100)) / (1+sim.run$ATratio)
    sim.run$nA <- (1-(sim.run$GC/100)) - sim.run$nT
    
    # obtain k-values from rhombus plot that fall within the 
    # "tolerance" values obtained from experimental values
    sim.run$A_minus_T <- sim.run$nA-sim.run$nT
    sim.run$G_minus_C <- sim.run$nG-sim.run$nC
    
    # return data
    cat(" Done!", "\n")
    return(sim.run)
  } else {
    stop("File does not exist!")
  }
}
################################################################################