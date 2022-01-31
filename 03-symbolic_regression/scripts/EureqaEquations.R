args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
TEST <- as.logical(as.character(args[2]))

suppressPackageStartupMessages(library(dplyr))
source("../lib/Equations.R")

if(TEST){
  df <- read.csv("../data/Test/simulation_batches_compliant.csv")
  df <- df[2:length(df)]
} else {
  df <- read.csv("../data/Training/simulation_batches_compliant.csv")
}

rate.const <- colnames(df)

# obtain rate constant files
for(i in 1:length(rate.const)){
  write.table(x = Equations(data = df, var = rate.const[i]), 
              sep = ",",
              file = ifelse(TEST, 
                            paste0("../data/Test/EUREQA_prediction_", 
                            rate.const[i], ".csv"),
                            paste0("../data/Training/EUREQA_prediction_", 
                            rate.const[i], ".csv")))
}