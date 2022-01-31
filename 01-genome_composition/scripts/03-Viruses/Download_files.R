args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
species <- as.character(args[2])
setwd(my_path)

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
source("../../lib/valid_url.R")

#-----------------------------
# download raw sequences in fasta format
# big files need more time to download
if(getOption('timeout') < 10000){
  options(timeout = 10000)
}

to.download <- "https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/?fq=%7B!tag=SeqType_s%7DSeqType_s:(%22Nucleotide%22)&fq=%7B!tag=VirusLineage_ss%7DVirusLineage_ss:(%22DNA%20viruses%22)&cmd=download&sort=SourceDB_s%20desc,CreateDate_dt%20desc,id%20asc&dlfmt=fasta&fl=AccVer_s,Definition_s,Nucleotide_seq"

cat("Downloading fasta files...")
# allow several attempts to donwload fasta files
# as internet disruptions may cause failure
attempt <- 1
while(attempt <= 3){
    url.test <- valid_url(to.download)
    cat("URL Test: ", url.test, "\n")

    if(isTRUE(url.test)){
        try(download.file(to.download, paste0("../../data/", species,"/Raw/sequences.fasta")))
        
        # If file doesn't exist, it'll likely be a server problem so allow it to retry up to 3 times
        # If it requires >3 times, likely problems are:
        #   (1) no internet connection
        #   (2) another pattern that's unaccounted for
        if(!file.exists(paste0("../../data/", species,"/Raw/sequences.fasta"))){
            Sys.sleep(3)
            attempt <- attempt + 1
        } else {
            break
        }
    } 
    if(isFALSE(url.test)){
        try(download.file(to.download, paste0("../../data/", species,"/Raw/sequences.fasta")))
        
        if(!file.exists(paste0("../../data/", species,"/Raw/sequences.fasta"))){
            Sys.sleep(3)
            attempt <- attempt + 1
        } else { 
            break
        }
    }
}
cat("Done!", "\n")

# unzip all files if needed
cat("Unzipping fasta file...")
system("gunzip ../../data/03-viruses/Raw/*")
cat("Done!", "\n")

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

df <- data.frame(
    "Species_name"  = Species_name,
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
    "AT_skew"       = AT_skew
)

# remove any NAs if present
df <- df[complete.cases(df),]

write.csv(
    df, 
    file = paste0("../../data/", species, "/All/all_download_dataframe.csv"), 
    row.names = FALSE
)
cat("Dataframe of species names and associated meta-data saved!", "\n")