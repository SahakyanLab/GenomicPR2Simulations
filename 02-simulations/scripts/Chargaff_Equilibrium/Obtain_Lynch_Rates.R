args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
TREK.SCALE <- as.logical(as.character(args[2]))
setwd(my_path)

# Load required supplementary functions and packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(suppressWarnings(library(Biostrings)))
source("../../../01-genome_composition/lib/valid_url.R")

# obtain mutation rates from Table 1 in the following paper
# https://www.pnas.org/content/107/3/961

# Import calculated mutation rates from trek paper
note.two <- read.csv("../../data/Raw/Trek-paper-Note-2-mutation-rates.csv",
                     header = TRUE)
lynch.rates <- read.csv("../../data/Raw/Michael_Lynch/Lynch-2010-mutation-rates_FREQUENCY.csv",
                        header = TRUE)

lynch.rates <- lynch.rates[match(note.two$MUT, lynch.rates$MUT),]

# match rate constants from lynch with trek
lynch.rates <- lynch.rates %>%
  mutate(MEAN = note.two$MEAN, .after = 1)

lynch.plot <- lynch.rates %>%
  ggplot(aes(x = MEAN)) +
  geom_point(aes(y = H.sapiens, color = "H.sapiens")) +
  geom_line(aes(y = H.sapiens, color = "H.sapiens")) +
  geom_point(aes(y = D.melanogaster, color = "D.melanogaster")) +
  geom_point(aes(y = C.elegans, color = "C.elegans")) +
  geom_point(aes(y = A.thaliana, color = "A.thaliana")) +
  geom_point(aes(y = S.cerevisiae, color = "S.cerevisiae")) +
  geom_point(aes(y = E.coli, color = "E.coli")) +
  xlim(0, 1.2) + 
  ylim(0, 0.6) + 
  labs(x = "Trek (Mean), byr",
       y = "Lynch, 2010")

ggsave("../../figures/Chargaff_Equilibrium/LynchRateFreq.pdf", plot = lynch.plot)

# conversion of rate frequencies to time bound rate constants
Coef <- lynch.rates %>%
  lm(formula = MEAN ~ H.sapiens-1) %>%
  coefficients

conversion.formula <- function(x){
  return(Coef*x)
}

final.rates <- lynch.rates %>%
  select(1:3) %>%
  cbind(apply(lynch.rates[4:length(lynch.rates)], 2, conversion.formula)) %>%
  mutate(
    MEDIAN = note.two$MEDIAN,
    SD = note.two$SD,
    .after = MEAN
  )

# save converted rate constants as csv
write.table(x = final.rates, 
            file = "../../data/Raw/Michael_Lynch/Lynch-2010-converted-mutation-rates.csv", 
            sep = ",", row.names = FALSE)
cat("Converted mutation rates to time domain for simulation!", "\n")

# fasta files to download 
species.names <- c(
  "Homo_sapiens", 
  "Arabidopsis_thaliana", 
  "Caenorhabditis_elegans", 
  "Drosophila_melanogaster",
  "Escherichia_coli", 
  "Saccharomyces_cerevisiae",
  "Photorhabdus_luminescens",
  "Teredinibacter_turnerae",
  "Mycobacterium_smegmatis",
  "Pseudomonas_fluorescens",
  "Rhodosporidium_toruloides",
  "Mus_musculus",
  "Daphnia_pulex",
  "Pristionchus_pacificus",
  "Daphnia_magna",
  "Pan_troglodytes",
  "Aotus_nancymaae"
)

download.files.url <- c(
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/115/GCF_000157115.1_ASM15711v1/GCF_000157115.1_ASM15711v1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/756/235/GCA_000756235.1_ASM75623v1/GCA_000756235.1_ASM75623v1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/102/985/GCF_900102985.1_IMG-taxon_2597490348_annotated_assembly/GCF_900102985.1_IMG-taxon_2597490348_annotated_assembly_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/964/255/GCF_000964255.1_ASM96425v1/GCF_000964255.1_ASM96425v1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/349/145/GCF_013349145.1_ASM1334914v1/GCF_013349145.1_ASM1334914v1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/414/285/GCF_001414285.1_ATCC948-1/GCF_001414285.1_ATCC948-1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/320/785/GCF_000320785.1_RHOziaDV1.0/GCF_000320785.1_RHOziaDV1.0_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/911/175/335/GCA_911175335.1_PA42_4.2/GCA_911175335.1_PA42_4.2_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/180/635/GCA_000180635.4_El_Paco_v._4/GCA_000180635.4_El_Paco_v._4_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/631/705/GCF_020631705.1_ASM2063170v1.1/GCF_020631705.1_ASM2063170v1.1_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/880/755/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_genomic.fna.gz",
  "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/952/055/GCF_000952055.2_Anan_2.0/GCF_000952055.2_Anan_2.0_genomic.fna.gz"
)

# big files need more time to download
if(getOption('timeout') < 1000){
  options(timeout = 1000)
}

cat("Downloading fasta files...")
# download files
for(i in 1:length(species.names)){
  # allow several attempts to donwload fasta files
  # as internet disruptions may cause failure
  print(paste0("Downloading fasta file of species: ", species.names[i],"..."), quote = FALSE)
  attempt <- 1
  while(attempt <= 3){
    url.test <- valid_url(download.files.url[i])
    print(paste0("URL Test: ", url.test), quote = F)
    
    if(isTRUE(url.test)){
      try(download.file(download.files.url[i], 
                        paste0("../../data/Raw/Michael_Lynch/", i, "-", species.names[i], ".dna.fna.gz")))
      
      # If file doesn't exist, it'll likely be a server problem so allow it to retry up to 3 times
      # If it requires >3 times, likely problems are:
      #   (1) no internet connection
      #   (2) another pattern that's unaccounted for
      if(!file.exists(paste0("../../data/Raw/Michael_Lynch/", i, "-", species.names[i], ".dna.fna.gz"))){
        Sys.sleep(3)
        attempt <- attempt + 1
      } else {
        break
      }
    } 
    if(isFALSE(url.test)){
        try(download.file(download.files.url[i], 
                          paste0("../../data/Raw/Michael_Lynch/", i, "-", species.names[i], ".dna.fna.gz")))
      
      if(!file.exists(paste0("../../data/Raw/Michael_Lynch/", i, "-", species.names[i], ".dna.fna.gz"))){
        Sys.sleep(3)
        attempt <- attempt + 1
      } else { 
        break
      }
    }
  }
}
cat("Done!", "\n")

# unzip all files if needed
cat("Unzipping all fasta files...")
system("gunzip ../../data/Raw/Michael_Lynch/*")
cat("Done!", "\n")

# read fasta files
files <- list.files(path = "../../data/Raw/Michael_Lynch/", 
                    pattern = ".fna$")
files <- str_sort(files, numeric = TRUE)

# obtain GC content from species
cat("Obtaining GC content and GC skews from species...")
base.values <- lapply(files, function(x){
  fai <- readDNAStringSet(paste0("../../data/Raw/Michael_Lynch/", x))
  
  # Initialise data frame for base calculations
  genome_length <- width(fai)
  
  # extract base contents
  all.letters <- letterFrequency(fai, letters="ACGT", OR=0)
  
  # obtain GC average
  all.letters <- cbind(all.letters,genome_length)
  all.letters <- colSums(all.letters)
  G_plus_C  <- (all.letters["G"]+all.letters["C"])/all.letters["genome_length"]
  G_minus_C <- (all.letters["G"]-all.letters["C"])/all.letters["genome_length"]
  GC_skew <- G_minus_C/G_plus_C

  A_plus_T  <- (all.letters["A"]+all.letters["T"])/all.letters["genome_length"]
  A_minus_T <- (all.letters["A"]-all.letters["T"])/all.letters["genome_length"]
  AT_skew <- A_minus_T/A_plus_T

  return(
    list(G_plus_C, GC_skew, AT_skew)
  )
})

GC.content <- sapply(base.values, `[[`, 1)
GC.skew <- sapply(base.values, `[[`, 2)
AT.skew <- sapply(base.values, `[[`, 3)

other.species.gc.avg <- data.frame(
  species = files,
  GC.content = unlist(GC.content),
  GC.skew = unlist(GC.skew),
  AT.skew = unlist(AT.skew)
)

other.species.gc.avg %>% 
  dplyr::select(-AT.skew) %>% 
  write.csv(file = "../../data/Raw/Michael_Lynch/GC_values.csv",
            row.names = TRUE)
cat("Done!", "\n")

if(TREK.SCALE){
  lynch.rates <- lynch.rates %>% 
    dplyr::select(1:3) %>% 
    cbind(apply(lynch.rates[4:length(lynch.rates)], 2, conversion.formula)) %>% 
    as_tibble() %>% 
    mutate(
      RATES = c("j", "n", "l", "i", "k", "m"),
      .after = MUT
    )
} else {
  lynch.rates <- lynch.rates %>% 
    as_tibble() %>% 
    mutate(
      RATES = c("j", "n", "l", "i", "k", "m"),
      .after = MUT
    )
}

# calculate rate constant ratios
ind.n <- which(lynch.rates$RATES == "n")
ind.j <- which(lynch.rates$RATES == "j")
ind.m <- which(lynch.rates$RATES == "m")
ind.i <- which(lynch.rates$RATES == "i")

# assign new GC contents
rate.constant.ratios <- sapply(4:length(lynch.rates), function(x){
  sum(lynch.rates[ind.n, x],lynch.rates[ind.j, x])/
    sum(lynch.rates[ind.i, x],lynch.rates[ind.m, x])
})

# assign ratios to new column 
other.species.gc.avg <- data.frame(
  Species = colnames(lynch.rates)[4:length(lynch.rates)],
  Rates = rate.constant.ratios,
  GC.average = other.species.gc.avg$GC.content*100,
  GC.skew = other.species.gc.avg$GC.skew,
  AT.skew = other.species.gc.avg$AT.skew
)

# save new data frame as csv
other.species.gc.avg %>%
  dplyr::select(-c(GC.skew, AT.skew)) %>% 
  write.csv(
    file = paste0("../../data/Raw/Michael_Lynch", 
    ifelse(TREK.SCALE, "/Trek_scale_", "/"), 
    "GC_vs_Rates.csv"),
    row.names = FALSE
  )

other.species.gc.avg %>%
  dplyr::select(-GC.average) %>% 
  write.csv(
    file = paste0("../../data/Raw/Michael_Lynch", 
    ifelse(TREK.SCALE, "/Trek_scale_", "/"), 
    "GC_AT_skew_vs_Rates.csv"),
    row.names = FALSE
  )
  
# save new data frame as csv
if(TREK.SCALE){
  lynch.rates %>%
  select(-c(2:3)) %>%
    write.csv(
      file = "../../data/Raw/Michael_Lynch/Trek_scale_Lynch-2010-mutation-rates_FREQUENCY.csv", 
      row.names = FALSE
    )
}