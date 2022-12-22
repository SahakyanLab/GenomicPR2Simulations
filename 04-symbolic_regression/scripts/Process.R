# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])

# import dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(stringr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))

my.path="/Users/paddy/Documents/DPhil/01-Chargaff/04-symbolic_regression/scripts"
setwd(my.path)
source("../lib/AnalyseSR.R")

for(trek_scale in c(FALSE, TRUE)){
    sym_reg <- AnalyseSR$new(trek_scale = trek_scale)
    sym_reg$run_process()
}