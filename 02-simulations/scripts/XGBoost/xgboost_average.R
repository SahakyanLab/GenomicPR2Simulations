suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(pbapply))
pbo = pboptions(type="txt")

args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[3])
setwd(my_path)

# load data as tibble
load.data <- function(file.path, scaling = 1){

  # Function to load the dataset as tibble class

  # Flag       Format        Description
  # file.path  <character>   Path to file for the dataset
  # scaling    <numeric>     Scaling factor for the standard deviation of the 
  #                          random drawing of mutation rate constants

  print(paste0("Loading scaling ",scaling,"..."), quote = FALSE)
  if(file.exists(file.path)){
    # load file
    if(grepl(pattern = ".csv", x = file.path, fixed = TRUE)){
      sim.run <- read.csv(file = file.path, header = TRUE)
    } else if (grepl(pattern = ".Rdata", x = file.path, fixed = TRUE)){
      load(file.path)
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
    return(sim.run)
  } else {
    stop("File does not exist!")
  }
}

filter.df <- function(dataset, species = "eukaryotes"){

  # Apply the PR-2 tolerance to the compliance region of the datasets

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # species  <character> Apply the PR-2 tolerance zone of 
  #                      "prokaryotes", "eukaryotes" or "viruses" organisms

  if(species == "prokaryotes"){
    file.name <- "../../../01-genome_composition/data/01-Prokaryotes/PR_compliance/PR2_fluctuations.csv"
  } else if(species == "eukaryotes"){
    file.name <- "../../../01-genome_composition/data/02-Eukaryotes/data/PR_compliance/PR2_fluctuations.csv"
  } else if(species == "viruses"){
    file.name <- "../../../01-genome_composition/data/03-Viruses/data/PR_compliance/PR2_fluctuations.csv"
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
             Label = fct_recode(Label, 
                                "NO" = "0", 
                                "YES" = "1"))
    
    # select only columns of Label and all rate constants
    sim.run.rates <- y %>%
      select(7:(dim(dataset)[2]-5)) %>%
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

non.compliant.sampling <- function(dataset, seed = 123){

  # Train the xgboost model

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # seed     <numeric>   Random number generator, useful for reproducible random objects

  # separate input variables
  sim.run.df <- dataset[[1]]
  df.compliant <- dataset[[2]]
  
  # random sampling
  df.non.compliant <- sim.run.df %>% 
    filter(Label == "NO")
  
  set.seed(seed = seed)
  to.sample <- sample(x = nrow(df.non.compliant),
                      size = nrow(df.compliant),
                      replace = FALSE)
  
  df.non.compliant <- sim.run.df %>%
    filter(Label == "NO") %>%
    slice(to.sample)
  
  # combine into df
  sim.run.df <- rbind(df.compliant, df.non.compliant)
  
  # Randomly shuffle rows of tbl
  df.train <- sim.run.df[sample(nrow(sim.run.df)),]
  
  # return data
  return(df.train)
}

xgb.train <- function(dataset, ncpu = 4, cv = 6, cvrep = 1, seed = 123, parallel = FALSE){

  # Train the xgboost model

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # cv       <numeric>   Number of cross-validation processes
  # cvrep    <numeric>   Number of times to repeat the cross-validation processes
  # seed     <numeric>   Random number generator, useful for reproducible random objects
  # parallel <boolean>   Run xgboost training process in parallel execution 

  # check if classes are correct
  if(ncpu%%1!=0){
    stop("Number of cores must be a positive integer.")
  }
  if(cv%%1!=0){
    stop("Cross-validation must be a positive integer.")
  }
  if(cvrep%%1!=0){
    stop("cvrep must be a positive integer.")
  }
  if(seed%%1!=0){
    stop("Seed must be a positive integer.")
  }

  # min / max normalisation
  normalize <- function(x, na.rm = TRUE) {
    return((x - min(x))/(max(x)-min(x)))
  }
  data <- dataset %>% 
    mutate(across(where(is.numeric), normalize))
  
  # set random seed
  set.seed(seed)
  
  # set up parallel processing
  if(parallel){
    cl = makeCluster(ncpu)
    registerDoParallel(cl)
    print(paste0("Running with ", ncpu, " cores..."), quote = F)
  }
  
  # Message output
  print(paste0("Validation scheme: ", cv, "-fold CV"), quote = F)
  
  # hyper-parameter search
  xgb.grid <- expand.grid(
    nrounds = c(1000), # boosting iterations
    max_depth = c(6), # max tree depth
    colsample_bytree = c(1), # subsample ratio of columns
    eta = c(0.1), # shrinkage
    gamma = c(0.1), # min loss reduction
    min_child_weight = c(1), # min sum of instance weight
    subsample = c(0.4) # subsample percentage
  )
  
  # define seeds for each CV process
  seeds <- lapply(1:((cv*cvrep)+1), FUN = function(x){
    sample.int(n = 50000, dim(xgb.grid)[1])
  })
  
  # create train control
  xgb.trcontrol <- trainControl(
    method = "repeatedcv",
    repeats = cvrep,
    number = cv,
    seeds = seeds,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    returnData = TRUE,
    verboseIter = TRUE,
    allowParallel = TRUE,
    savePredictions = TRUE
  )
  
  # train model
  print("Model fitting...", quote = F)
  xgb.model <- train(
    Label ~ ., 
    data = data,
    method = "xgbTree",
    trControl = xgb.trcontrol,
    tuneGrid = xgb.grid,
    metric = "ROC",
    verbose = FALSE
  )
  
  # stop cluster
  if(parallel){
    stopCluster(cl)
  }
  
  # best model
  best.model <- xgb.model$results %>% 
    filter(ROC == max(ROC))
  
  # return model
  return(list(xgb.model, best.model))
}

# Run multiple iterations to obtain average feature importance
# load data
sim.run <- load.data(file.path = "../../data/Main_Simulation/Non_Symmetric-uniform-scaling-1.Rdata", 
                    scaling = 1)

# obtain compliant cases
df <- filter.df(dataset = sim.run, 
                species = "viruses")
  
df[[2]] %>%
  write.csv(file = "../../../03-symbolic_regression/data/Training/simulation_batches_compliant.csv",
            row.names = FALSE)

# set up random number generator
set.seed(2021)
rng <- sample(x = 1:1000, size = 1000, replace = FALSE)

df.samples <- pblapply(rng, function(x){
  # random sampling of non compliant cases
  df.train <- non.compliant.sampling(dataset = df, 
                                     seed = x)
  
  # train model with optimised hyper-parameters
  model <- xgb.train(data = df.train, 
                     ncpu = 4, 
                     cv = 5, 
                     cvrep = 1, 
                     seed = 123, 
                     parallel = FALSE)
  
  # return feature importance as data frame
  return(varImp(model[[1]])$importance %>%
           rownames_to_column() %>%
           rename(rates = rowname,
                  importance = Overall) %>%
           arrange(rates) %>%
           select(importance))
})
featImp.all <- do.call(cbind, df.samples)

# rename columns
colnames(featImp.all) <- paste0(
  "output_", seq(from = 1, to = dim(featImp.all)[2], by = 1)
  )

# sort rate constants alphabetically
rates.sorted <- sort(colnames(df[[1]])[2:ncol(df[[1]])])

# transpose tibble
featImp.all <- featImp.all %>%
  as_tibble() %>%
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  pivot_wider(names_from = rowname,
              values_from = value) %>%
  select(-name) %>%
  rename_with(~ rates.sorted)

# statistical calculations
data <- rates.sorted %>%
  as_tibble() %>%
  mutate(mean = unname(colMeans(featImp.all)),
         sd   = unname(apply(featImp.all, 2, sd)),
         max  = mean+sd,
         min  = mean-sd) %>%
  rename(rate = value) 

# plot average feature importance
data.plot <- data %>%
  arrange(max) %>%
  mutate(rate = forcats::fct_inorder(rate)) %>%
  ggplot() + 
  geom_bar(aes(x = rate, y = mean),
           stat = "identity", fill = "skyblue", alpha = 0.7) + 
  geom_errorbar(aes(x = rate, ymin = min, ymax = max),
                width = 0.4, colour = "orange", alpha = 1, size = 1.5) +
  coord_flip()
  labs(x = "Importance",
       y = "Features",
       title = "Feature Importance averaged over 1000 iterations")

ggsave(filename = "../../figures/XGBoost/xgbtree_featimport_average_plots.pdf",
       plot = data.plot)
  