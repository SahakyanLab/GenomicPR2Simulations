args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[2])
setwd(my_path)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(pbapply))
pbo = pboptions(type="txt")
source("../../lib/LoadData.R")

filter.df <- function(dataset, species = "eukaryotes"){

  # Apply the PR-2 tolerance to the compliance region of the datasets

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # species  <character> Apply the PR-2 tolerance zone of 
  #                      "prokaryotes", "eukaryotes" or "viruses" organisms

  if(species == "prokaryotes"){
    file.name <- "../../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv"
  } else if(species == "eukaryotes"){
    file.name <- "../../../01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv"
  } else if(species == "viruses"){
    file.name <- "../../../01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv"
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
      dplyr::select(7:(dim(dataset)[2]-5)) %>%
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

xgb.train <- function(dataset, best.model, ncpu = 4, seed = 123){

  # Train the xgboost model

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate equations
  # seed     <numeric>   Random number generator, useful for reproducible random objects

  # check if classes are correct
  if(ncpu%%1!=0 | ncpu<=0) stop("Number of cores must be a positive integer.")
  if(seed%%1!=0 | seed<=0) stop("Seed must be a positive integer.")

  # min / max normalisation
  normalize <- function(x, na.rm = TRUE) {
    return((x - min(x))/(max(x)-min(x)))
  }

  data <- dataset %>% 
    mutate(across(where(is.numeric), normalize))
  
  set.seed(seed)
  
  # hyper-parameter search
  xgb.grid <- expand.grid(
    nrounds = best.model$nrounds, # boosting iterations
    max_depth = best.model$max_depth, # max tree depth
    colsample_bytree = best.model$colsample_bytree, # subsample ratio of columns
    eta = best.model$eta, # shrinkage
    gamma = best.model$gamma, # min loss reduction
    min_child_weight = best.model$min_child_weight, # min sum of instance weight
    subsample = best.model$subsample # subsample percentage
  )
  
  xgb.trcontrol <- trainControl(
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    returnData = TRUE,
    verboseIter = FALSE,
    allowParallel = TRUE,
    savePredictions = TRUE
  )
  
  # train model
  xgb.model <- train(
    Label ~ ., 
    data = data,
    method = "xgbTree",
    trControl = xgb.trcontrol,
    tuneGrid = xgb.grid,
    metric = "ROC",
    verbose = FALSE
  )
  
  # best model
  best.model <- xgb.model$results %>% 
    filter(ROC == max(ROC))
  
  # return model
  return(list(xgb.model, best.model))
}

# Run multiple iterations to obtain average feature importance
# load xgboost model file
xgbmodel <- readRDS(file = paste0("../../data/XGBoost/xgbtree.model"))
best.model <- xgbmodel$results %>% filter(ROC == max(ROC))
rm(xgbmodel)

# load data
sim.run <- LoadData(file.path = "../../data/Main_Simulation/Non_Symmetric-uniform/Non_Symmetric-uniform-scaling-1.Rdata", 
                    scaling = 1)

# obtain compliant cases
df <- filter.df(dataset = sim.run, species = "eukaryotes")

df[[2]] %>%
  select(-Label) %>%  
  write.csv(file = "../../../03-symbolic_regression/data/Train/simulation_batches_compliant.csv",
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
                     best.model = best.model, 
                     ncpu = ncpu, 
                     seed = 123)
  
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
  tidyr::pivot_longer(-rowname) %>%
  tidyr::pivot_wider(names_from = rowname,
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
  arrange(mean) %>%
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