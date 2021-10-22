################################################################################
states <- function(iterations = 10, seed = 2021, NCPU = 6){
  # Function to obtain min and max single nucleotide content from real-life species

  # Dependencies
  #     Simulation.R

  # Flag       Format     Description
  # iterations <integer>  Number of systems to generate
  #                       in the simulation. Value needs to
  #                       be above zero.
  # seed       <integer>  Random number generator, useful for
  #                       reproducible random objects
  # NCPU       <integer>  Number of cores to use for parallel
  #                       execution. 

  df <- c(eukaryotes.df$G_content,
          eukaryotes.df$C_content,
          eukaryotes.df$A_content,
          eukaryotes.df$T_content,
          prokaryotes.df$G_content,
          prokaryotes.df$C_content,
          prokaryotes.df$A_content,
          prokaryotes.df$T_content)
  
  # randomly select the order of single nucleotide bases
  rng <- set.seed(seed)
  
  random.state <- sapply(1:iterations, function(x){
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
}
################################################################################