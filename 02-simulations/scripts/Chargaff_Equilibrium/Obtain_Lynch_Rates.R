args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
setwd(my_path)

# Load required supplementary functions and packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Biostrings))

# obtain mutation rates from Table 1 in the following paper
# https://www.pnas.org/content/107/3/961

# Import calculated mutation rates from trek paper
note.two <- read.csv("../../data/Raw/Trek-paper-Note-2-mutation-rates.csv",
                     header = TRUE)
if(!file.exists("../../data/Raw/Michael_Lynch/Lynch-2010-mutation-rates_FREQUENCY.csv")){
  stop("File does not exist. Please obtain from the paper.")
} 

lynch.rates <- read.csv("../../data/Raw/Michael_Lynch/Lynch-2010-mutation-rates_FREQUENCY.csv",
                        header = TRUE)

lynch.rates <- lynch.rates[match(note.two$MUT, lynch.rates$MUT),]

# match rate constants from lynch with trek
lynch.rates <- lynch.rates %>%
  mutate(MEAN = note.two$MEAN) %>%
  relocate(MEAN, .after = 1)

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

ggsave("../../figures/Chargaff_Equilibrium/Michael_Lynch/LynchRateFreq.pdf",
       plot = lynch.plot)

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
  mutate(MEDIAN = note.two$MEDIAN,
         SD = note.two$SD) %>%
  relocate(c(MEDIAN, SD), .after = MEAN)

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

# obtain mean GC content from species
cat("Obtaining mean GC content from species...")
GC.average <- sapply(files, function(x){
  fai <- readDNAStringSet(paste0("../../data/Raw/Michael_Lynch/", x))
  
  # Initialise data frame for base calculations
  genome_length <- width(fai)
  
  # extract base contents
  all.letters <- letterFrequency(fai, letters="ACGT", OR=0)
  
  # obtain GC average
  all.letters <- cbind(all.letters,genome_length)
  all.letters <- colSums(all.letters)
  G_plus_C  <- (all.letters["G"]+all.letters["C"])/all.letters["genome_length"]

  return(G_plus_C)
})

other.species.gc.avg <- as.data.frame(GC.average)
other.species.gc.avg %>% 
  write.csv(file = "../../data/Raw/Michael_Lynch/GC_average.csv",
            row.names = TRUE)
cat("Done!", "\n")

lynch.rates <- lynch.rates %>%
  as_tibble() %>%
  mutate(RATES = c("n", "m", "l", "i", "j", "k")) %>%
  relocate(RATES, .after = MUT)

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
  GC.average = other.species.gc.avg$GC.average*100
)

# save new data frame as csv
other.species.gc.avg %>%
  write.csv(file = "../../data/Raw/Michael_Lynch/GC_vs_Rates.csv",
            row.names = FALSE)