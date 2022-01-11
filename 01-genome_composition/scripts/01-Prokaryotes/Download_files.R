##############################################################################################
# This file downloads the text file of all bacteria species in the Ensembl Genome database,  
# filters the data by omitting duplication to avoid skewing the analysis and figures,        
# downloads each species one by one.                                                         
# Calculations are made per species:                                                         
#   G+C, G-C, A+T, A-T, GC skews, AT skews, genome length, base content of G, C, A and T.    
#############################################################################################
args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
species <- as.character(args[2])
setwd(paste0(my_path, "/", species, "/"))

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(Biostrings))

# Download bacteria genome name file obtained from:
# ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48
if(!file.exists(paste0("../../data/", species, "/Raw/species_EnsemblBacteria.txt"))){
  print("Downloading text file of all Ensembl Bacteria names...", quote = F)
  
  download.file(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/",
                       "species_EnsemblBacteria.txt"),
                paste0("../../data/", species, "/Raw/species_EnsemblBacteria.txt"))
} else {
  print("species_EnsemblBacteria.txt file exists!", quote = F)
}

f <- read.delim(paste0("../../data/", species, "/Raw/species_EnsemblBacteria.txt"), sep = "\t", header = T)
write.csv(f, paste0("../../data/", species, "/Raw/species_EnsemblBacteria.csv"), row.names = TRUE)

#-----------------------------
# While files are systematically named, some files are stored using assembly instead of taxonomy_id
# and there seems to be no obvious difference. Possibly error on website,
# hence need to choose between taxonomy and assembly names later in the URL query
species.names <- character() 
taxonomy.id   <- character()
assembly      <- character()
address.all   <- character()
f.rownames    <- character()
for(i in 0:183){
  bacteria.number <- f[which(f$other_alignments==paste0("bacteria_",i,"_collection_core_48_101_1")),]
  address <- paste0("bacteria_",i,"_collection","/", bacteria.number$X.name,"/dna/")
  
  # URL address requires capitalized string unless it starts with an underscore
  bacteria.number$X.name <- str_replace(string = bacteria.number$X.name,
                                        pattern = "^(?!^_).*$", replacement = str_to_title)
  
  address.all   <- c(address.all, address)
  species.names <- c(species.names, bacteria.number$X.name)
  taxonomy.id   <- c(taxonomy.id,
                     gsub(pattern = "#|:|/| ", replacement = "_", bacteria.number$taxonomy_id))
  assembly      <- c(assembly,
                     gsub(pattern = "\\.", replacement = "_", bacteria.number$assembly))
  f.rownames    <- c(f.rownames, rownames(bacteria.number))
}

#-----------------------------
# filter the ~44k species names to omit duplicates
###########################################################################
###### CAUTION: THE BELOW IS NOT OPTIMISED - THIS TAKES HOURS TO RUN ######
###########################################################################
df <- data.frame("f.rownames"    = f.rownames,
                 "species.names" = species.names,
                 "taxonomy.id"   = taxonomy.id,
                 "assembly"      = assembly,
                 "address.all"   = address.all
                 )

print("Cleaning up all files...", quote = FALSE)
sorted.df <- sort(df$species.names)
sorted.df <- str_replace(string = sorted.df, pattern = "^(^_)", replacement = "")
filtered.species <- character()
pb  <- txtProgressBar(min = 1, max = length(species.names), style=3) 
i <- 1
while(i<=length(species.names)){
  split.string <- strsplit(sorted.df[i], split = "_")[[1]]
  
  if(length(split.string)>1){
    pos <- which(grepl(pattern = split.string[1], x = sorted.df) & 
                   grepl(pattern = split.string[2], x = sorted.df))
  } else {
    pos <- which(grepl(pattern = split.string[1], x = sorted.df))
  }
  compare.integer <- sort(table(sapply(pos, function(i) nchar(i))), decreasing = TRUE)
  
  final.pos <- tail(pos, n=1) # default is the final element of list 
  if(length(names(compare.integer))>1){
    if((max(as.integer(names(compare.integer)))-min(as.integer(names(compare.integer))))>1){
      # take most frequently occurring element
      freq.integer <- as.integer(names(compare.integer)[1])
      final.pos <- tail(which(nchar(pos)==freq.integer), n=1) 
    }
    if(max(diff(pos))-min(diff(pos))>100){
      # check that differences aren't big
      final.pos <- pos[which.max(diff(pos))]
    }
  }
  filtered.species <- c(filtered.species, sorted.df[pos[1]])
  i <- final.pos+1
  
  setTxtProgressBar(pb, i)
}
close(pb)

# re-insert underscore in front of strings starting with lower case letters
filtered.species <- gsub("^([a-z])", "\\_\\1", filtered.species,  perl = TRUE)

# re-define data frame
df <- df[match(filtered.species, df$species.names),]

# remove any possible duplicates
df <- df[-which(duplicated(df$species.names)==TRUE),]

# removing any remaining extra strings we missed
sorted.df        <- sort(df$species.names)
sorted.df        <- str_replace(string = sorted.df, pattern = "^(^_)", replacement = "")
split.string     <- strsplit(sorted.df, split = "_")
filtered.species <- character()

pb  <- txtProgressBar(min = 1, max = length(split.string), style=3) 
for(i in 1:length(split.string)){
  check <- agrep(pattern = split.string[i], split.string, max.distance = 0.1)
  filtered.species <- c(filtered.species, df$species.names[check])
  
  if(length(check)>1){
    split.string[check[2:length(check)]] <- NULL
  }
  setTxtProgressBar(pb, i)
}
close(pb)

###########################################################################
###### CAUTION: THE ABOVE IS NOT OPTIMISED - THIS TAKES HOURS TO RUN ######
###########################################################################

# re-define data frame
filtered.species <- gsub("^([a-z])", "\\_\\1", filtered.species,  perl = TRUE)
df <- df[match(filtered.species, df$species.names),]
df <- df[-which(duplicated(df$species.names)==TRUE),]

# removing special cases which occurred a few times 
df <- df[-grep("unidentified|uncultured|unknown|untyped|
                                 |cloning|clones|patent|Treatment|Oligonucleotide|
                                 |methods|involving|vectors|means|event|initiation|
                                 |improvements|relating|identification|highly|
                                 |procedure|treating|amplification|joining|
                                 |complementing|comprising|predict|
                                 |\\bleft\\b|\\bright\\b|\\bwith\\b|\\bspecific\\b|
                                 |\\.{3}|\\=", 
                                 df$Species_name, ignore.case = TRUE), ]

# remove demonstrative pronouns which occurred a few times
df <- df[-grep("this|that|these|those|such|here|there|thereof|
                                 |\\buse\\b|\\bit\\b", 
                                 df$Species_name, ignore.case = TRUE), ]

# removing any possible NAs
df <- df[-which(is.na(df$genome_length)),]
write.csv(df, file = paste0("../../data/", species, "/All/all_download_dataframe.csv"), row.names = FALSE)
print("Cleaned up all files!", quote = FALSE)

#-----------------------------
# Initialise content vectors
G_content     <- numeric()
C_content     <- numeric()
A_content     <- numeric()
T_content     <- numeric()
G_plus_C      <- numeric()
G_minus_C     <- numeric()
A_plus_T      <- numeric()
A_minus_T     <- numeric()
genome_length <- numeric()

# Function to check if a given URL exists or not
valid_url <- function(url_in, t = 300){

  # Function to check if a given URL exists or not

  # Flag      Format       Description
  # url_in   <character>   Character vector of the URL to download
  # t        <numeric>     Maximum time until timeout reached

  con <- url(url_in)
  check <- suppressWarnings(try(open.connection(con, open = "rt", timeout = t), silent = T)[1])
  suppressWarnings(try(close.connection(con), silent = T))
  ifelse(is.null(check), TRUE, FALSE)
}

# big files need more time to download
if(getOption('timeout') < 1000){
  options(timeout = 2000)
}

for(i in 1:length(df$species.names)){
  # ----------------------------------------------
  print("Loading the toplevel assembly from the following address:", quote=F)
  print(df$address.all[i], quote = F)
  
  # ----------------------------------------------
  # First, check if a given URL exists or not;
  # Sometimes a server error comes up so my solution is to put the try block 
  # in a loop and allow it to retry in case it generates an error
  attempt <- 1
  while(attempt <= 3){
    url.test <- valid_url(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/fasta/",
                                  df$address.all[i],
                                  df$species.names[i],
                                  ".",
                                  df$taxonomy.id[i],
                                  ".dna.toplevel.fa.gz"))
    print(paste0("URL Test: ", url.test), quote = F)
    
    if(isTRUE(url.test)){
      try(download.file(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/fasta/",
                                  df$address.all[i],
                                  df$species.names[i],
                                  ".",
                                  df$taxonomy.id[i],
                                  ".dna.toplevel.fa.gz"), 
                        paste0("./", df$species.names[i], ".dna.toplevel.fa.gz")))
      
      # If file doesn't exist, it'll likely be a server problem so allow it to retry up to 3 times
      # If it requires >3 times, likely problems are:
      #   (1) no internet connection
      #   (2) another pattern that's unaccounted for
      if(!file.exists(paste0(df$species.names[i], ".dna.toplevel.fa.gz"))){
        Sys.sleep(3)
        attempt <- attempt + 1
      } else {
        break
      }
    } 
    if(isFALSE(url.test)){
        try(download.file(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/fasta/",
                                    df$address.all[i],
                                    df$species.names[i],
                                    ".",
                                    df$assembly[i],
                                    ".dna.toplevel.fa.gz"), 
                          paste0("./", df$species.names[i], ".dna.toplevel.fa.gz")))
      
      if(!file.exists(paste0(df$species.names[i], ".dna.toplevel.fa.gz"))){
        Sys.sleep(3)
        attempt <- attempt + 1
        
        # If still nothing downloaded until this point, then 
        # likely a pattern that hasn't been accounted for;
        # Web scrape the url address and download
        if(attempt==3){
          new.url <- webscrape.filename(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/fasta/",
                                      df$address.all[i]))
          
          try(download.file(paste0("ftp://ftp.ensemblgenomes.org/pub/bacteria/release-48/fasta/",
                                  df$address.all[i],
                                  new.url,
                                  ".dna.toplevel.fa.gz"),
                        paste0(new.url, ".dna.toplevel.fa.gz")))
          
          df$species.names[i] <- new.url
          break
        }
      } else { 
        break
      }
    }
  }

  # ----------------------------------------------
  print("Unpacking the toplevel combined fasta file...", quote=F)
  system(paste0("gunzip ", df$species.names[i], ".dna.toplevel.fa.gz"))

  # ----------------------------------------------
  print("File unpacked. Now calculating...", quote=F)
  
  # More memory efficient to request a small subset of sequences per iteration
  fai <- fasta.index(paste0(df$species.names[i], ".dna.toplevel.fa"))
  fai <- readDNAStringSet(fai)
  
  # Initialise data frame for base calculations
  vec.length <- width(fai)
  
  # extract base contents
  all.letters <- letterFrequency(fai, letters="ACGT", OR=0)
  
  # obtain other metadata of each species
  G_content     <- c(G_content, sum(all.letters[,"G"])/sum(vec.length))
  C_content     <- c(C_content, sum(all.letters[,"C"])/sum(vec.length))
  A_content     <- c(A_content, sum(all.letters[,"A"])/sum(vec.length))
  T_content     <- c(T_content, sum(all.letters[,"T"])/sum(vec.length))
  G_plus_C      <- c(G_plus_C, sum(all.letters[,"G"]+all.letters[,"C"])/sum(vec.length))
  G_minus_C     <- c(G_minus_C, sum(all.letters[,"G"]-all.letters[,"C"])/sum(vec.length))
  A_plus_T      <- c(A_plus_T, sum(all.letters[,"A"]+all.letters[,"T"])/sum(vec.length))
  A_minus_T     <- c(A_minus_T, sum(all.letters[,"A"]-all.letters[,"T"])/sum(vec.length))
  genome_length <- c(genome_length, sum(vec.length))
  
  print(paste0("G+C content:   ", G_plus_C[i]), quote = F)
  print(paste0("A+T content:   ", A_plus_T[i]), quote = F)
  print(paste0("Genome length: ", genome_length[i]), quote = F)
  
  # ----------------------------------------------
  # Save space on the active memory and hard disk
  rm(seq.string)
  rm(all.letters)
  file.remove(paste0(df$species.names[i], ".dna.toplevel.fa"))

  # ----------------------------------------------
  cat(paste0("\n\033[0;", 32, "m", "Loaded ", i, "/",
             length(df$species.names), " prokaryotes","\033[0m","\n\n"))

  # ----------------------------------------------
  # Prevent excessive server querying breakdown
  Sys.sleep(3)
}
print("All fasta files downloaded and meta-data obtained!", quote = FALSE)

#-----------------------------
# Calculate skews and ratios
GC_skew  <- G_minus_C/G_plus_C
AT_skew  <- A_minus_T/A_plus_T
GC_ratio <- G_content/C_content
AT_ratio <- A_content/T_content

# Combine all vectors into a new data frame
filtered.df <- data.frame("rowname"           = df$f.rownames,
                          "Species_name"      = df$species.names,
                          "G_content"         = G_content,
                          "C_content"         = C_content,
                          "A_content"         = A_content,
                          "T_content"         = T_content,
                          "G_plus_C"          = G_plus_C,
                          "G_minus_C"         = G_minus_C,
                          "GC_skew"           = GC_skew,
                          "A_plus_T"          = A_plus_T,
                          "A_minus_T"         = A_minus_T,
                          "AT_skew"           = AT_skew,
                          "GC_ratio"          = GC_ratio,
                          "AT_ratio"          = AT_ratio,
                          "genome_length"     = genome_length
                          )

# If any NAs in the vectors, indices to download saved in vector
# input vector in the for-loop in "download_files.R
if(dim(filtered.df[rowSums(is.na(filtered.df))!=0,])[1]!=0){
  to.download <- as.numeric(row.names(filtered.df[rowSums(is.na(filtered.df))!=0,]))
} else {
  print("No NAs in the data frame!", quote = F)
}

write.csv(filtered.df, file = paste0("../../data/", species, "/All/all_filtered_dataframe.csv"), row.names = FALSE)
print("Dataframe of species names and associated meta-data saved!", quote = FALSE)