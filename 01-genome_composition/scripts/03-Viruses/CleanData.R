args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
species <- as.character(args[2])
setwd(my_path)

suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
suppressPackageStartupMessages(library(dplyr))

df <- read.csv(file = paste0("../../data/", species, "/All/all_download_dataframe.csv"), header=TRUE)

#-----------------------------
# conduct same filtering process as with viruses in case of any duplication
cat("Cleaning up all files...", "\n")
sorted.df <- sort(df$Species_name)
sorted.df <- sorted.df[grep(pattern = "^[a-zA-Z]", x = sorted.df)]
sorted.df <- sorted.df[which(str_locate(sorted.df, " ")[,1]!=2)]
filtered.species <- rep(NA, nrow(df))

pb  <- txtProgressBar(min = 1, max = nrow(df), style=3) 
i <- 1
j <- 1

while(i <= length(df$Species_name)){
  split.string <- strsplit(sorted.df[i], split = " ")[[1]]
  split.ind    <- setdiff(1:length(split.string), grep(pattern = "([()])|&)", x = split.string))
  
  possibleError <- tryCatch(
    which(
      grepl(pattern = split.string[split.ind[1]], x = sorted.df) & 
      grepl(pattern = split.string[split.ind[2]], x = sorted.df)),
    error = function(e){e}
  )
  
  if(!inherits(possibleError, "error")){
    # First filtering
    if(length(split.string)>1){
      pos <- which(
        grepl(pattern = split.string[split.ind[1]], x = sorted.df) & 
        grepl(pattern = split.string[split.ind[2]], x = sorted.df)
        )
    } else {
      pos <- which(grepl(pattern = split.string[split.ind[1]], x = sorted.df))
    }
    
    if(length(split.string)>1){
      pos <- which(
        grepl(pattern = split.string[split.ind[1]], x = sorted.df) & 
        grepl(pattern = split.string[split.ind[2]], x = sorted.df)
        )
    } else {
      pos <- which(grepl(pattern = split.string[split.ind[1]], x = sorted.df))
    }    
    
    # Second filtering
    concat.string <- paste0("^", split.string[1],"+[[:space:]]+")
    pos           <- pos[grep(concat.string, sorted.df[pos])]
    
    # Third filtering
    pos <- pos[which(unlist(lapply(strsplit(sorted.df[pos], split = " "), "[[", 2))==split.string[2])] 
  } else {
    filtered.species[j] <- sorted.df[i]
    i<-i+1
    next
  }

  if(!(length(pos)>0)){
    filtered.species[j] <- sorted.df[pos[1]]
    i <- i+1
    next
  }
  
  compare.integer <- sort(table(nchar(pos)), decreasing = TRUE)
  compare.integer <- as.integer(names(compare.integer))
  
  final.pos <- tail(pos, n=1) # default is the final element of list 

  if(length(compare.integer)>1){
    if(diff(compare.integer)>1){
      # take most frequently occurring element
      final.pos <- tail(which(nchar(pos)==compare.integer[[1]]), n=1) 
    }

    if(max(diff(pos))-min(diff(pos))>100){
      # check that differences aren't big
      final.pos <- pos[which.max(diff(pos))]
    }
  }

  filtered.species[j] <- sorted.df[pos[1]]
  i <- final.pos+1
  j <- j+1
  
  setTxtProgressBar(pb, i)
}
close(pb)

# Remove all caps strings
check.allcaps <- grep("^[[:upper:][:space:]]+$", filtered.species)
if(length(check.allcaps) > 0){
  filtered.species <- filtered.species[-check.allcaps]
}

write.csv(
  filtered.species, 
  file = paste0("../../data/", species,"/All/filtered_species_backup.csv"), 
  row.names = FALSE
)

# Remove any possible duplication
string.split.one <- unlist(lapply(strsplit(filtered.species, split = " "), "[", 1))
string.split.two <- unlist(lapply(strsplit(filtered.species, split = " "), "[", 2))
patterns         <- paste(string.split.one, string.split.two)
filtered.species <- filtered.species[unique(match(patterns, patterns))]

# Sort and re-define data frame
filtered.species <- sort(filtered.species)
df <- df[match(filtered.species, df$Species_name),]

# remove any lower-cased letters
check.lowercase <- grep("(^[[:lower:]])", df$Species_name)
if(length(check.lowercase) > 0){
  df <- df[-grep("(^[[:lower:]])", df$Species_name), ]
}

# remove words whose 1st word contains special characters
first.word <- str_extract(df$Species_name, "([^\\s]+\\s)")
df <- df[is.na(str_extract(first.word, "([-?.,;:'_+=()!@#$%^&*|~`{}])")), ]

# remove all-caps strings
check.allcaps <- grep("^[^a-z]*$", df$Species_name)
if(length(check.allcaps) > 0){
  df <- df[-grep("^[^a-z]*$", df$Species_name), ]
}

# removing any remaining extra strings we missed
sorted.df <- sort(df$Species_name)
sorted.df <- str_replace(string = sorted.df, pattern = "^(^_)", replacement = "")

# allow dynamic-changing max value of progress bar
progressbar <- function(vector.len){
  txtProgressBar(min = 1, max = length(vector.len), style = 3)
}

for(i in 1:length(sorted.df)){
  pb <- progressbar(sorted.df)
  matches <- agrep(pattern = sorted.df[i], x = sorted.df, max.distance = 1)
  
  if(matches[1]==i & length(matches)>1){
    sorted.df <- sorted.df[-matches[2:length(matches)]]
  }
  setTxtProgressBar(pb, i)
}
close(pb)

df <- df[match(sorted.df, df$Species_name),]

# removing special cases which occurred a few times
check.case <- grep("^WO|JP|KR|Sequence", df$Species_name)
if(length(check.case) > 0){
  df <- df[-check.case, ]
}

check.case <- grep(
  "unidentified|uncultured|unknown|untyped|
  |cloning|clones|patent|Treatment|Oligonucleotide|
  |methods|involving|vectors|means|event|initiation|
  |improvements|relating|identification|highly|
  |procedure|treating|amplification|joining|
  |complementing|comprising|predict|comparison|
  |delivery|test|usable|detecting|
  |\\bleft\\b|\\bright\\b|\\bwith\\b|\\bspecific\\b|
  |\\bcryoem\\b|\\bcrystal\\b|\\bproduction\\b|
  |\\btrapping\\b|\\bmodified\\b|\\bloaded\\b|
  |attachment|
  |\\bon\\b|\\busing\\b|\\.{3}|\\=", 
  df$Species_name, ignore.case = TRUE)

if(length(check.case) > 0){
  df <- df[-check.case, ]
}

# remove demonstrative pronouns which occurred a few times
check.case <- grep(
  "this|that|these|those|such|here|there|thereof|
  |\\buse\\b|\\bit\\b", df$Species_name, ignore.case = TRUE)

if(length(check.case) > 0){
  df <- df[-check.case, ]
}

# removing individual cases
df <- df[-which(df$Species_name %in% head(df$Species_name[order(df$AT_ratio, decreasing = TRUE)], n=2)),]

# remove any possible NAs in the remaining dataframe
df <- df[complete.cases(df),]

# additional removals
df <- as_tibble(df) %>%
  mutate(id = seq(from = 1, to = dim(df)[1], by = 1))

to.remove <- c(11, 102, 587, 909, 1688, 2220, 2281, 2334, 2339, 3252, 3865, 3915)
df <- df[-to.remove,]

write.csv(
  df, 
  file = paste0("../../data/", species, "/All/all_filtered_dataframe.csv"), 
  row.names = FALSE
)
cat("Cleaned up all files!", "\n")