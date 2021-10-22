suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(suppressWarnings(library(geneplotter)))
colfun = colorRampPalette(c("white","blue","skyblue",
                            "chartreuse3","green","yellow",
                            "orange","red","darkred"))

args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save.as <- as.character(args[2])
scaling <- as.numeric(args[3])
Tolerance <- as.character(args[4])
Tolerance <- as.logical(Tolerance)
setwd(my_path)

load(paste0("../../data/Main_Simulation/Strand_Symmetric-normal/Strand_Symmetric-normal-scaling-",scaling,".Rdata"))
sim.run$nC <- (sim.run$GC/100)/(1+sim.run$GCratio)
sim.run$nG <- (sim.run$GC/100)-sim.run$nC
sim.run$nT <- (1-(sim.run$GC/100)) / (1+sim.run$ATratio)
sim.run$nA <- (1-(sim.run$GC/100)) - sim.run$nT

to.keep <- which(is.na(str_extract(string = names(sim.run), pattern = ".bn.")))
sim.run <- sim.run[,to.keep]

#----------------------
# G+C content
GC_hist <- function(dataset, scaling){

  # Plot of the GC content as a function of the mutation rate constants 

  # Flag      Format     Description
  # dataset   <Rdata>     Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # scaling   <numeric>   Scaling factor of the simulation results

  Plot <- dataset %>%
    as_tibble() %>%
    ggplot(aes(x = GC, 
               y = ..density..)) + 
    geom_histogram(fill = "skyblue",
                   color = "black",
                   alpha = 1,
                   breaks = seq(min(dataset$GC),
                                max(dataset$GC),
                                length.out = 51)) + 
    labs(x = paste0("GC hist (%)", "\n",
                    "Mean = ", format(round(mean(dataset$GC), 2), 
                                      nsmall = 2), " ",
                    "St.Dev = ", format(round(sd(dataset$GC), 2), 
                                        nsmall = 2)),
         y = "Density",
         title = scaling)
  return(Plot)
}

ggsave(width=15, height=8, 
       filename = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                         scaling,"/GC_hist.", save.as),
       plot = GC_hist(sim.run, paste0("Scaling = ", scaling)))

#----------------------
# A+T content
AT_hist <- function(dataset, scaling){

  # Plot the AT content as histograms

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results

  Plot <- dataset %>%
    as_tibble() %>%
    ggplot(aes(x = AT, 
               y = ..density..)) + 
    geom_histogram(fill = "skyblue",
                   color = "black",
                   alpha = 1,
                   breaks = seq(min(dataset$AT),
                                max(dataset$AT),
                                length.out = 51)) + 
    labs(x = paste0("AT hist (%)", "\n",
                    "Mean = ", format(round(mean(dataset$AT), 2), 
                                      nsmall = 2), " ",
                    "St.Dev = ", format(round(sd(dataset$AT), 2), 
                                        nsmall = 2)),
         y = "Density",
         title = scaling)
  return(Plot)
}

ggsave(width=15, height=8, 
       filename = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                         scaling,"/AT_hist.", save.as),
       plot = AT_hist(sim.run, paste0("Scaling = ", scaling)))

#----------------------
# G-C vs. A-T 
difference.plots <- function(dataframe, scaling){

  # Plot the G-C vs. A-T difference plots

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results
  
  x = dataframe$nA-dataframe$nT
  y = dataframe$nG-dataframe$nC
  
  smoothScatter(y=y, x=x,
                nrpoints=100, nbin=1000,
                bandwidth=c(diff(range(x))/500, diff(range(y))/500),
                # xlim=c(-1,1),ylim=c(-1,1),
                xlab="A-T", ylab="G-C", main=scaling,
                colramp=colfun,
                cex.axis = 1.2, cex.lab = 1.1)
}

if(save.as == "pdf"){
  pdf(width=15, height=8, 
      paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
             scaling,"/G-C_vs_A-T.", save.as))
} else {
  png(width=600, height=600, 
      paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
             scaling,"/G-C_vs_A-T.", save.as))
}
difference.plots(sim.run, paste0("Scaling = ", scaling))
pic.saved <- dev.off()

#-----------------------------
# GC-skew vs. AT-skew
skew.plot <- function(dataframe, scaling){

  # Plot the GC vs. AT skew plots

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results

  smoothScatter(y=dataframe$GCskew, 
                x=dataframe$ATskew,
                nrpoints=100, nbin=1000,
                xlab="AT skew", ylab="GC skew", main=scaling,
                colramp=colfun,
                xlim=c(-1e-15,1e-15),
                ylim=c(-1e-15,1e-15),
                cex.axis = 1.2, cex.lab = 1.1)
  
  abline(v = 0, lty = 2, lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=15, height=8, 
      paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
             scaling,"/GC-skew_vs_AT-skew.", save.as))  
} else {
  png(width=600, height=600, 
      paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
             scaling,"/GC-skew_vs_AT-skew.", save.as))  
}
skew.plot(sim.run, paste0("Scaling = ", scaling))
pic.saved <- dev.off()

#----------------------
# G+C content vs. rates
GC.content.rates <- function(
  dataframe, 
  scaling,
  species = "Prokaryotes",
  Tolerance = FALSE 
){

  # Plot of the GC content as a function of the mutation rate constants 

  # Flag      Format     Description
  # dataset   <Rdata>     Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # scaling   <numeric>   Scaling factor of the simulation results

  # Chargaff tolerance
  sim.run <- dataframe
  if(Tolerance){
    if(species == "Prokaryotes"){
      tolerance  <- read.csv(file = "../../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv", 
                             header = TRUE)
    } else if(species == "Eukaryotes"){
      tolerance  <- read.csv(file = "../../../01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv", 
                             header = TRUE)
    } else if(species == "Viruses"){
      tolerance  <- read.csv(file = "../../../01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv", 
                             header = TRUE)
    }
    at.tol <- tolerance[tolerance$metadata == "AT_skew", "st.dev"]
    gc.tol <- tolerance[tolerance$metadata == "GC_skew", "st.dev"]
    ind <- which(abs(sim.run$GCskew-0) <= gc.tol & abs(sim.run$ATskew-0) <= at.tol)
    rates <- c("kag", "kac", "kct", "kca", "kat", "kcg")
    
    plots <- lapply(1:length(rates), function(x){
      x.values <- sim.run[ind,rates[x]]
      y.values <- sim.run$GC[ind]
      return(smoothScatter(y=y.values, x=x.values,
                           nrpoints=100, 
                           nbin=400,
                           bandwidth=c(diff(range(x.values))/200, 
                                       diff(range(sim.run$GC))/200),
                           xlab=paste(rates[x], ", byr-1"),
                           ylab="G+C content, %",
                           colramp=colfun, 
                           ylim=c(0,100),
                           cex.axis = 1.4, cex.lab = 1.4))
    })
  } else {
    rates <- c("kag", "kac", "kct", "kca", "kat", "kcg") 
    
    plots <- lapply(1:length(rates), function(x){
      x.values <- sim.run[,rates[x]]
      return(smoothScatter(y=sim.run$GC, x=x.values,
                           nrpoints=100, 
                           nbin=400,
                           bandwidth=c(diff(range(x.values))/200, 
                                       diff(range(sim.run$GC))/200),
                           xlab=paste(rates[x], ", byr-1"),
                           ylab="G+C content, %",
                           colramp=colfun, 
                           ylim=c(0,100),
                           cex.axis = 1.4, cex.lab = 1.4))
    })
  }
  return(plots)
}

if(Tolerance){
  if(save.as == "pdf"){
    pdf(width=6, height=15, 
        file = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                      scaling,"/tolerance-GC-content_vs_rates-all.", save.as))
  } else {
    png(width = 1200, height = 1100,
        file = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                      scaling,"/tolerance-GC-content_vs_rates-all.", save.as))
  }
  par(mfcol=c(6,3))
  GC.content.rates(dataframe = sim.run, species = "Eukaryotes", 
                   scaling = paste0("Scaling = ", scaling), Tolerance = Tolerance)
  GC.content.rates(dataframe = sim.run, species = "Prokaryotes", 
                   scaling = paste0("Scaling = ", scaling), Tolerance = Tolerance)
  GC.content.rates(dataframe = sim.run, species = "Viruses", 
                   scaling = paste0("Scaling = ", scaling), Tolerance = Tolerance)
  pic.saved <- dev.off() 
} else {
  if(save.as == "pdf"){
    pdf(width=6, height=15, 
        file = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                      scaling,"/GC-content_vs_rates-all.", save.as))
  } else {
    png(width = 600, height = 1100,
        file = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                      scaling,"/GC-content_vs_rates-all.", save.as))
  }
  par(mfcol=c(6,1))
  GC.content.rates(dataframe = sim.run, 
                   scaling = paste0("Scaling = ", scaling), Tolerance = Tolerance)
  pic.saved <- dev.off() 
}
