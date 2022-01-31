args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[2])
setwd(my_path)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(MLeval))
source("../../lib/LoadData.R")

filtered.df <- function(dataset, seed = 123, species = "eukaryotes"){

  # Apply the PR-2 tolerance to the compliance region of the datasets

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # seed     <numeric>   Random number generator, useful for reproducible random objects
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
      as_tibble() %>%
      na.omit()
    
    # filter out all compliant cases
    df.compliant <- sim.run.rates %>%
      filter(Label == "YES")

    # random sampling
    df.non.compliant <- sim.run.rates %>% 
      filter(Label == "NO")

    set.seed(seed = seed)
    to.sample <- sample(x = nrow(df.non.compliant),
                        size = nrow(df.compliant),
                        replace = F)

    df.non.compliant <- sim.run.rates %>%
      filter(Label == "NO") %>%
      slice(to.sample)

    # combine into df
    sim.run.rates <- rbind(df.compliant, df.non.compliant)

    # Randomly shuffle rows of tbl
    df.train <- sim.run.rates[sample(nrow(sim.run.rates)),]
    
    # return data
    return(list(df.train, to.sample))
  } else {
    stop("File does not exist!")
  }
}

xgb.train <- function(dataset, ncpu = 4, cv = 6, cvrep = 1, seed = 123){

  # Train the xgboost model

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # cv       <numeric>   Number of cross-validation processes
  # cvrep    <numeric>   Number of times to repeat the cross-validation processes
  # seed     <numeric>   Random number generator, useful for reproducible random objects

  # check if classes are correct
  if(ncpu%%1!=0 | ncpu<=0) stop("Number of cores must be a positive integer.")
  if(cv%%1!=0 | cv<=0) stop("Cross-validation must be a positive integer.")
  if(cvrep%%1!=0 | cvrep<=0) stop("cvrep must be a positive integer.")
  if(seed%%1!=0 | seed<=0) stop("Seed must be a positive integer.")

  # extract only the sampled data set
  data <- dataset[[1]]
  
  # min / max normalisation
  normalize <- function(x, na.rm = TRUE) {
    return((x- min(x)) /(max(x)-min(x)))
  }
  
  data <- data %>% 
    mutate(across(where(is.numeric), normalize))
  
  set.seed(seed)
  cat("Validation scheme:", cv, "-fold CV", "\n")

  # hyper-parameter search
  xgb.grid <- expand.grid(
    nrounds = c(200, 500, 1000, 2000, 3000, 
                4000, 5000, 7000, 9000, 11000, 
                13000, 15000, 17000, 20000), # boosting iterations
    max_depth = c(5, 6, 8, 10), # max tree depth
    colsample_bytree = 1, # subsample ratio of columns
    eta = c(0.005, 0.01, 0.02, 0.1), # shrinkage
    gamma = c(0, 0.1, 0.2), # min loss reduction
    min_child_weight = c(1, 5, 10), # min sum of instance weight
    subsample = c(0.3, 0.4, 0.6, 0.8) # subsample percentage
  )
  
  # define seeds for each CV process
  seeds <- lapply(1:((cv*cvrep)+1), function(x){
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
  cat("Model fitting...", "\n")
  xgb.model <- train(
    Label ~ ., 
    data = data,
    method = "xgbTree",
    trControl = xgb.trcontrol,
    tuneGrid = xgb.grid,
    metric = "ROC",
    verbose = TRUE
  )
  
  # save AUC ROC curve
  res <- evalm(xgb.model)
  res.plot <- res$roc
  ggsave(paste0("../../figures/XGBoost/xgbtree_aucroc_plots.pdf"),
        device = "pdf",
        plot = res.plot)
  
  # feature importance
  feat.imp.plot <- varImp(xgb.model)$importance %>%
    rownames_to_column() %>%
    rename("variable" = rowname) %>%
    arrange(Overall) %>%
    filter(row_number() > max(row_number())-10) %>%
    mutate(variable = forcats::fct_inorder(variable)) %>%
    ggplot() +
    geom_col(aes(x = variable, y = Overall)) +
    coord_flip() +
    scale_fill_grey() +
    labs(x = "Importance",
         y = "Features")
  
  ggsave(paste0("../../figures/XGBoost/xgbtree_featimport_plots.pdf"),
        device = "pdf",
        plot = feat.imp.plot)

  # save best model
  saveRDS(xgb.model, file = paste0("../../data/XGBoost/xgbtree.model"))

  # best model
  best.model <- xgb.model$results %>% 
    filter(ROC == max(ROC))
  
  # return model
  return(list(xgb.model, best.model))
}

sim.run <- LoadData(file.path = "../../data/Main_Simulation/Non_Symmetric-uniform-scaling-1.Rdata", 
                    scaling = 1)
df.train <- filtered.df(dataset = sim.run, seed = 2022, species = "eukaryotes")
model <- xgb.train(dataset = df.train, ncpu = ncpu, cv = 6, cvrep = 1, seed = 1234)

# save model results plot
trellis.device(
  width = 14,
  height = 10,
  device = "pdf",
  file = "../../figures/XGBoost/xgbtree_roc-hyperparams_plots.pdf")
plot(model[[1]])
plot.save <- dev.off()