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
scaling <- as.numeric(args[2])
setwd(my_path)

# load simulation data set
load(paste0("../../data/Main_Simulation/Strand_Symmetric-normal/Strand_Symmetric-normal-scaling-",scaling,".Rdata"))
to.keep <- which(!is.na(str_extract(string = names(sim.run), pattern = "skew")))
sim.run <- sim.run[,to.keep]
sim.run <- sim.run[,c(3:10, 1:2)]

# extract year timelines
bn.years <- str_extract(string = names(sim.run), pattern = ".+(?=.bn.)")
bn.years <- unique(bn.years)
bn.years <- c(bn.years[1:(length(bn.years)-1)], "four_28")

#-----------------------------
# GC-skew vs. AT-skew
first.col <- seq(from = 1, to = length(sim.run), by = 2)
second.col <- seq(from = 2, to = length(sim.run), by = 2)

for(i in 1:length(bn.years)){
  skew.plot.evolution <- function(dataframe, scaling){
    
    # Plot the GC vs. AT skew plots
    
    # Flag      Format     Description
    # dataset  <Rdata>     Dataset of the equilibrium outputs from the
    #                      numerically solved kinetic mutation rate
    #                      equations.
    # scaling  <numeric>   Scaling factor of the simulation results
    
    smoothScatter(y=dataframe[,first.col[i]], 
                  x=dataframe[,second.col[i]],
                  nrpoints=100, nbin=1000,
                  xlab="AT skew", ylab="GC skew", main=scaling,
                  colramp=colfun,
                  xlim=c(-1e-15,1e-15),
                  ylim=c(-1e-15,1e-15),
                  cex.axis = 1.2, cex.lab = 1.1)
    
    abline(v = 0, lty = 2, lwd = 2)
    abline(h = 0, lty = 2, lwd = 2)
  }
  
  png(width=600, height=600, 
      paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
             scaling,"/evolution_GC-skew_vs_AT-skew_", bn.years[i], "bn.png"))  
  skew.plot.evolution(sim.run, paste0("Scaling = ", scaling))
  dev.off()
}


# make animation movie out of strand symmetry evolution images
suppressPackageStartupMessages(library(magick))
suppressPackageStartupMessages(library(animation))

# load photos for video clip
SLF.photos <- list.files(path = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                                       scaling, "/"),
                         pattern = "evolution", 
                         all.files = TRUE, ignore.case = TRUE)
SLF.photos <- SLF.photos[c(3,5,4,2,1)]

#This extracts the underlying height, width, and type of image.
img.height <- magick::image_info(
  image_read(paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                    scaling, "/", SLF.photos[1])))$height
img.width <- magick::image_info(
  image_read(paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                    scaling, "/", SLF.photos[1])))$width
img.type <- magick::image_info(
  image_read(paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                    scaling, "/", SLF.photos[1])))$format

# display each picture for 0.25 seconds
animation::ani.options(interval = 0.25,
                       ani.height = img.height,
                       ani.width = img.width,
                       ani.dev = tolower(img.type),
                       ani.type = tolower(img.type))

# increasing video dimensions for better image quality
opts <- paste("-s ", img.height * 1.5, "x", img.width * 1.5, sep = "")

animation::saveVideo(
  for(i in 1:length(SLF.photos)){
    SLF.image <- magick::image_read(
      paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                        scaling, "/", SLF.photos[i]))
    plot(SLF.image)
  },
  video.name = paste0(file="../../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                      scaling, "/SkewEvolutionAnimation.mp4"),
  other.opts = paste0("-pix_fmt yuv420p -b 300k ", opts)
)
