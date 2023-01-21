# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
species <- as.character(args[2])
seed <- as.numeric(args[3])
ml_seed <- as.numeric(args[4])
avg_ml_seed <- as.numeric(args[5])
cv <- as.numeric(args[6])
cvrep <- as.numeric(args[7])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(caret)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(MLeval)))

setwd(my.path)
source("../lib/xgboostModel.R")

train_model <- xgboostModel$new(
    species = species,
    cv = cv, 
    cvrep = cvrep
)
train_model$find_best_model(
    seed = seed,
    ml_seed = ml_seed
)
train_model$get_avg_behaviour(
    seed = seed,
    avg_ml_seed = avg_ml_see
)