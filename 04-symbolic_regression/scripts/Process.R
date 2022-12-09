# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

setwd(my.path)
source("../lib/AnalyseSR.R")

sym_reg <- AnalyseSR$new(trek_scale = TRUE)
sym_reg$run_process()
