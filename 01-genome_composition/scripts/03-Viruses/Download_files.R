##############################################################################################
# This file downloads the text file of all bacteria species in the Ensembl Genome database,  #
# filters the data by omitting duplication to avoid skewing the analysis and figures,        #
# downloads each species one by one.                                                         #
# Calculations are made per species:                                                         #
#   G+C, G-C, A+T, A-T, GC skews, AT skews, genome length, base content of G, C, A and T.    #
#                                                                                            #
# Further information in the README.md file                                                  #
##############################################################################################
# read arguments from job submission
args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
species <- as.character(args[2])
setwd(paste0(my_path, "/", species, "/"))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))

#-----------------------------
# More memory efficient to request a small subset of sequences per iteration
fai <- fasta.index(paste0("../../data/", species,"/Raw/sequences.fasta"))
fai <- readDNAStringSet(fai)

# Initialise data frame for base calculations
genome_length <- width(fai)

# extract species name
Species_name <- names(fai)
Species_name <- str_extract(string = Species_name, pattern = "(?<=\\|).*(?=,)")

# extract base contents
all.letters <- letterFrequency(fai, letters="ACGT", OR=0)/genome_length

# obtain other metadata of each species
G_plus_C  <- all.letters[,"G"]+all.letters[,"C"]
G_minus_C <- all.letters[,"G"]-all.letters[,"C"]
A_plus_T  <- all.letters[,"A"]+all.letters[,"T"]
A_minus_T <- all.letters[,"A"]-all.letters[,"T"]
GC_ratio  <- all.letters[,"G"]/all.letters[,"C"]
AT_ratio  <- all.letters[,"A"]/all.letters[,"T"]
GC_skew   <- G_minus_C/G_plus_C
AT_skew   <- A_minus_T/A_plus_T

df <- data.frame("Species_name"  = Species_name,
                 "A_content"     = all.letters[,"A"],
                 "C_content"     = all.letters[,"C"],
                 "G_content"     = all.letters[,"G"],
                 "T_content"     = all.letters[,"T"],
                 "genome_length" = genome_length,
                 "G_plus_C"      = G_plus_C,
                 "G_minus_C"     = G_minus_C,
                 "A_plus_T"      = A_plus_T,
                 "A_minus_T"     = A_minus_T,
                 "GC_ratio"      = GC_ratio,
                 "AT_ratio"      = AT_ratio,
                 "GC_skew"       = GC_skew,
                 "AT_skew"       = AT_skew)

# remove any NAs if present
df <- df[complete.cases(df),]

write.csv(filtered.df, file = paste0("../../data/", species, "/All/all_filtered_dataframe.csv"), row.names = FALSE)
print("Dataframe of species names and associated meta-data saved!", quote = FALSE)