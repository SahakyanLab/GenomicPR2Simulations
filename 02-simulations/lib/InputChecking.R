################################################################################
inputchecking <- function(Acont, Gcont, Ccont, span, step, 
                          max.runs, muttype, dist, species, 
                          scale.fac, tolerance, tol.return, 
                          sim.evol, NCPU, seed){

  # Function to check for the correct input parameters of the Simulation.R script

  # Dependencies
  #     Simulation.R

  # check if classes are correct
  if((!is.numeric(Acont) | !is.numeric(Gcont) | !is.numeric(Ccont))){
    stop("Base content need to be numeric values below zero.")
  }
  # check that base contents are fractions and all sum to one
  if((Acont%%1 == 0 | Gcont%%1 == 0 | Ccont%%1 == 0)){
    stop("Base contents need to be a fraction below zero.")
  }

  # run checks for max.runs
  if(max.runs < 1){
    stop("Max.runs needs to be a positive integer above zero.")
  } else if(max.runs%%1 != 0){
    stop("Max.runs needs to be an integer of positive value.")
  }

  # check if classes are correct
  if(!is.character(muttype)){
    stop("Muttype must be of character class.")
  }

  # check if classes are correct
  if(!is.character(dist)){
    stop("Muttype must be of character class.")
  }

  # check if classes are correct
  if(!is.character(species)){
    stop("Muttype must be of character class.")
  }

  # check if scale fac is a positive integer value
  if(scale.fac%%1 != 0){
    stop("Scale.fac needs to be of positive integer value.")
  }

  # remind user that scale.fac should be specific values for reproducing the paper
  if(!scale.fac %in% c(0,1,2,5,10)){
    cat("To reproduce the results in the paper, please use scale.fac = c(0,1,2,5,10).", "\n")
  }

  # check if classes are correct
  if(!is.logical(tolerance)){
    stop("Tolerance can only be TRUE or FALSE")
  } 

  # check if classes are correct
  if(!is.character(tol.return)){
    stop("tol.return must be of character class.")
  }

  # check if classes are correct
  if(!is.logical(sim.evol)){
    stop("sim.evol can only be TRUE or FALSE")
  } 

  # check if ncpu is a positive integer
  if(NCPU%%1 != 0){
    stop("NCPU must be a positive integer values for parallel execution.")
  }

  # check if seed is a positive integer
  if(seed%%1 != 0){
    stop("seed must be a positive integer values for producing random objects.")
  }
}
################################################################################