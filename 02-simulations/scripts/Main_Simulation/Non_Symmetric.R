args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save.as <- as.character(args[2])
scaling <- as.numeric(args[3])
dist <- as.character(args[4])
Tolerance <- as.logical(as.character(args[5]))
setwd(my_path)

if(!scaling %in% c(0,1,2,5,10)){
  stop("To reproduce the results in the paper, please use scaling = c(0,1,2,5,10).")
}

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(suppressWarnings(library(geneplotter)))

colfun = colorRampPalette(c("white","blue","skyblue",
                            "chartreuse3","green","yellow",
                            "orange","red","darkred"))

sim.run <- readRDS(file = paste0("../../data/Main_Simulation/Non_Symmetric-",dist,
"/Non_Symmetric-",dist,"-scaling-",scaling,".Rdata"))
sim.run$nC <- (sim.run$GC/100)/(1+sim.run$GCratio)
sim.run$nG <- (sim.run$GC/100)-sim.run$nC
sim.run$nT <- (1-(sim.run$GC/100)) / (1+sim.run$ATratio)
sim.run$nA <- (1-(sim.run$GC/100)) - sim.run$nT

#----------------------
# G+C content
GC_hist <- function(dataset, scaling){

  # Plot the GC content as histograms

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results

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
       filename = paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/Scaling_",
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
       filename = paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
                         scaling,"/AT_hist.", save.as),
       plot = AT_hist(sim.run, paste0("Scaling = ", scaling)))

#----------------------
# G-C vs. A-T
difference.plots <- function(dataset, scaling){

  # Plot the G-C vs. A-T difference plots

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results

  x = dataset$nA-dataset$nT
  y = dataset$nG-dataset$nC

  smoothScatter(y=y, x=x,
                nrpoints=100, nbin=1000,
                bandwidth=c(diff(range(x))/500, diff(range(y))/500),
                xlim=c(-1,1),ylim=c(-1,1),
                xlab="A-T", ylab="G-C", main=scaling,
                colramp=colfun,
                cex.axis = 1.2, cex.lab = 1.1)

  abline(v = 0, lty = 2, lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=6, height=6,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/G-C_vs_A-T.", save.as))
} else {
  png(width = 550, height = 550,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/G-C_vs_A-T.", save.as))
}
difference.plots(sim.run, paste0("Scaling = ", scaling))
pic.saved <- dev.off()

#-----------------------------
# GC-ratio vs. AT-ratio
ratio.plot <- function(dataset, scaling, xlim=c(0,100), ylim=c(0,100)){

  # Plot the GC vs. AT ratio plots

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results
  # xlim     <numeric>   X-axis limit for the plot of the simulation.
  # ylim     <numeric>   Y-axis limit for the plot of the simulation.

  smoothScatter(y=dataset$GCratio,
                x=dataset$ATratio,
                nrpoints=100, nbin=1000,
                xlim=xlim,ylim=ylim,
                xlab="AT-ratio", ylab="GC-ratio", main=scaling,
                colramp=colfun,
                cex.axis = 1.2, cex.lab = 1.1)

  abline(v = 1, lty = 2, lwd = 2)
  abline(h = 1, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=15, height=8,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-ratio_vs_AT-ratio.", save.as))
  ratio.plot(sim.run, paste0("Scaling = ", scaling))
  pic.saved <- dev.off()

  pdf(width=15, height=8,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-ratio_vs_AT-ratio.", save.as))
  ratio.plot(sim.run, paste0("Scaling = ", scaling), xlim=c(0,10), ylim=c(0,10))
  pic.saved <- dev.off()
} else {
  png(width = 850, height = 850,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-ratio_vs_AT-ratio.", save.as))
  ratio.plot(sim.run, paste0("Scaling = ", scaling))
  pic.saved <- dev.off()

  png(width = 550, height = 550,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-ratio_vs_AT-ratio_0-10-range.", save.as))
  ratio.plot(sim.run, paste0("Scaling = ", scaling), xlim=c(0,10), ylim=c(0,10))
  pic.saved <- dev.off()
}

#-----------------------------
# GC-skew vs. AT-skew
skew.plot <- function(dataset, scaling){

  # Plot the GC vs. AT skew plots

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # scaling  <numeric>   Scaling factor of the simulation results

  smoothScatter(y=dataset$GCskew,
                x=dataset$ATskew,
                nrpoints=100, nbin=1000,
                xlab="AT skew", ylab="GC skew", main=scaling,
                colramp=colfun,
                cex.axis = 1.2, cex.lab = 1.1)

  abline(v = 0, lty = 2, lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=6, height=6,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-skew_vs_AT-skew.", save.as))
} else {
  png(width = 550, height = 550,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-skew_vs_AT-skew.", save.as))
}
skew.plot(sim.run, paste0("Scaling = ", scaling))
pic.saved <- dev.off()

#----------------------
# Rate constants 
Fluc.tol <- readRDS(file = "../../data/Chargaff_Equilibrium/ChargaffEquilibrium.Rdata")
EQtolerance <- mean(Fluc.tol$Mean)*((1/100)*25)

rate.con.plots <- function(
  dataset = sim.run,
  species = "Prokaryotes", 
  Tolerance = FALSE
  ){

  # Plot of the mutation rate constants as ratio distribution plots

  # Flag      Format     Description
  # dataset   <Rdata>     Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # species   <character> Apply the PR-2 tolerance zone of 
  #                       "Prokaryotes", "Eukaryotes" or "Viruses" organisms
  # Tolerance <boolean>   Apply Chargaff tolerance to range of values if TRUE
  
  rates <- data.frame(col1 = c("kag", "kat", "kac", "kct", "kca", "kcg"),
                      col2 = c("ktc", "kta", "ktg", "kga", "kgt", "kgc"))
  length.out <- 10*5+1
  
  # Chargaff compliance
  if(Tolerance){
    if(species == "Prokaryotes"){
      tolerance  <- read.csv(
        file = "../../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv",
        header = TRUE
      )
    } else if(species == "Eukaryotes"){
      tolerance  <- read.csv(
        file = "../../../01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv", 
        header = TRUE
      )
    } else if(species == "Viruses"){
      tolerance  <- read.csv(
        file = "../../../01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv", 
        header = TRUE
      )
    }
    at.tol <- tolerance[tolerance$metadata == "AT_skew", "st.dev"]
    gc.tol <- tolerance[tolerance$metadata == "GC_skew", "st.dev"]
    ind <- which(abs(sim.run$GCskew-0) <= gc.tol & abs(sim.run$ATskew-0) <= at.tol)
    
    Plots <- lapply(1:6, function(x){
      to.plot <- sim.run %>%
      as_tibble() %>%
      dplyr::slice(ind) %>%
      dplyr::select(rates$col1[x], rates$col2[x]) %>%
      setNames(c('rate_1', 'rate_2')) %>%
      mutate(frac = rate_1/rate_2) %>%
      dplyr::select(frac) %>%
      dplyr::filter(frac <= 5)
    
    plots <- to.plot %>%
      ggplot(aes(x = frac, ..density..)) +
      geom_histogram(fill = "skyblue", 
                      color = "black",
                      alpha = 1,
                      breaks = seq(min(to.plot$frac), 
                                  max(to.plot$frac), 
                                  length.out = length.out)) +
      geom_vline(xintercept = 1,
                  linetype = "dashed",
                  size = 1) + 
      scale_fill_manual(values = c("#69b3a2")) + 
      {if(dist == "uniform")coord_cartesian(xlim = c(0, 5),
                                            ylim = c(0, 0.9))} +
      {if(dist == "normal")coord_cartesian(xlim = c(0, 5),
                                            ylim = c(0, 1.2))} +
      labs(x = "",
            y = "",
            title = paste0(rates$col1[x], "/", rates$col2[x]))
      
      return(plots)
    })
  } else {
    Plots <- lapply(1:6, function(x){
      to.plot <- sim.run %>%
        as_tibble() %>%
        dplyr::select(rates$col1[x], rates$col2[x]) %>%
        setNames(c('rate_1', 'rate_2')) %>%
        mutate(frac = rate_1/rate_2) %>%
        dplyr::select(frac) %>%
        dplyr::filter(frac <= 5)
      
      plots <- to.plot %>%
        ggplot(aes(x = frac, ..density..)) +
        geom_histogram(fill = "skyblue", 
                       color = "black",
                       alpha = 1,
                       breaks = seq(min(to.plot$frac), 
                                    max(to.plot$frac), 
                                    length.out = length.out)) +
        geom_vline(xintercept = 1,
                   linetype = "dashed",
                   size = 1) + 
        scale_fill_manual(values = c("#69b3a2")) + 
        {if(dist == "uniform")coord_cartesian(xlim = c(0, 5),
                                              ylim = c(0, 0.6))} +
        {if(dist == "normal")coord_cartesian(xlim = c(0, 5),
                                              ylim = c(0, 1))} +
        labs(x = "",
             y = "",
             title = paste0(rates$col1[x], "/", rates$col2[x]))
      return(plots)
    })
  }
  return(do.call(grid.arrange, c(Plots, nrow = 6)))
}

p1 <- rate.con.plots(dataset = sim.run, species = "Eukaryotes", Tolerance = Tolerance)
p2 <- rate.con.plots(dataset = sim.run, species = "Prokaryotes", Tolerance = Tolerance)
p3 <- rate.con.plots(dataset = sim.run, species = "Viruses", Tolerance = Tolerance)

grid.arrange(p1, p2, p3, ncol = 3) %>%
  ggsave(width = 18, height = 14, 
          filename = paste0(file ="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_", 
          scaling, ifelse(Tolerance, "/tolerance", "/allind"),
          "-rate_ratios.", save.as))

#----------------------
# G+C content vs. rates

GC.content.rates <- function(dataset, scaling){

  # Plot of the GC content as a function of the mutation rate constants 

  # Flag      Format     Description
  # dataset   <Rdata>     Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # scaling   <numeric>   Scaling factor of the simulation results

  sim.run <- dataset
  rates <- c("kag", "kat", "kac", "kga", "kgt", "kcg")
  plots <- lapply(1:length(rates), function(x){
    x.values <- sim.run[,rates[x]]
    return(smoothScatter(y=sim.run$GC, x=x.values,
                         nrpoints=100, 
                         nbin=400,
                         bandwidth=c(diff(range(x.values))/200, 
                                     diff(range(sim.run$GC))/200),
                         xlab=rates[x],
                         ylab="G+C content, %",
                         colramp=colfun, 
                         cex.axis = 1.5, cex.lab = 1.5))
  })
  return(plots)
}

if(save.as == "pdf"){
  pdf(width=6, height=15, 
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-content_vs_rates-all.", save.as))
} else {
  png(width = 500, height = 1300,
      paste0(file="../../figures/Main_Simulation/Non_Symmetric-",dist,"/scaling_",
             scaling,"/GC-content_vs_rates-all.", save.as))
}
par(mfcol=c(6,1))
GC.content.rates(dataset = sim.run, scaling = paste0("Scaling = ", scaling))
pic.saved <- dev.off()

# # ########################################################################################
# # import tolerance
# prok = read.csv("/Users/paddy/Documents/DPhil/01-Chargaff/01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv")
# euk = read.csv("/Users/paddy/Documents/DPhil/01-Chargaff/01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv")
# virus = read.csv("/Users/paddy/Documents/DPhil/01-Chargaff/01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv")

# other.species <- read.csv("/Users/paddy/Documents/DPhil/01-Chargaff/02-simulations/data/Raw/Michael_Lynch/Trek_scale_GC_AT_skew_vs_Rates.csv")
# other.species <- other.species[,-2]

# euk.species <- c(
#   "H.sapiens", "D.melanogaster",	"C.elegans",	"A.thaliana",	
#   "M.m.domesticus",	"D.pulex",	"P.pacificus",	
#   "D.magna",	"P.troglodytes",	"A.nancymaae"
# )

# prok.species <- c(
#   "E.coli",	"P.luminescens ATCC29999", "T.turnerae", 
#   "S.cerevisiae", "M.smegmatis",	"P.fluorescens ATCC948",	"R. toruloides"
# )

# other.species <- other.species %>% 
#   dplyr::mutate(
#     kingdom = ifelse(Species %in% euk.species, "eukaryotes", "prokaryotes"),
#     kingdom.col = ifelse(kingdom == "eukaryotes", "purple", "darkgreen")
#   )

# #helper function to make transparent ramps
# alpharamp<-function(alpha) {
#   stopifnot(alpha>=0 & alpha<=256)
#   function(n) paste(colorRampPalette(c("white","blue","skyblue",
#                             "chartreuse3","green","yellow",
#                             "orange","red","darkred"))(n), format(as.hexmode(alpha), upper.case=T), sep="")
# }

# skew.plot <- function(dataset, scaling){
#   smoothScatter(
#     y=dataset$GCskew,
#     x=dataset$ATskew,
#     nrpoints=100, nbin=1000,
#     xlab="AT skew", ylab="GC skew", main=scaling,
#     colramp= alpharamp(alpha=20),
#     cex.axis = 1.2, cex.lab = 1.1,
#     xlim = c(-0.016, 0.016),
#     ylim = c(-0.03, 0.03)
#   )

#   rect(
#     xleft = euk[euk[, "metadata"] == "AT_skew", "mean"]-euk[euk[, "metadata"] == "AT_skew", "st.dev"],
#     xright = euk[euk[, "metadata"] == "AT_skew", "mean"]+euk[euk[, "metadata"] == "AT_skew", "st.dev"],
#     ybottom = euk[euk[, "metadata"] == "GC_skew", "mean"]-euk[euk[, "metadata"] == "GC_skew", "st.dev"],
#     ytop = euk[euk[, "metadata"] == "GC_skew", "mean"]+euk[euk[, "metadata"] == "GC_skew", "st.dev"],
#     lty = 3,
#     lwd = 4,
#     border = "purple"
#   )

#   rect(
#     xleft = prok[prok[, "metadata"] == "AT_skew", "mean"]-prok[prok[, "metadata"] == "AT_skew", "st.dev"],
#     xright = prok[prok[, "metadata"] == "AT_skew", "mean"]+prok[prok[, "metadata"] == "AT_skew", "st.dev"],
#     ybottom = prok[prok[, "metadata"] == "GC_skew", "mean"]-prok[prok[, "metadata"] == "GC_skew", "st.dev"],
#     ytop = prok[prok[, "metadata"] == "GC_skew", "mean"]+prok[prok[, "metadata"] == "GC_skew", "st.dev"],
#     lty = 3,
#     lwd = 4,
#     border = "darkgreen"
#   )

#   par(new=TRUE)
#   plot(
#     x=other.species$AT.skew,
#     y=other.species$GC.skew,
#     axes=FALSE,
#     cex=2,
#     pch=20,
#     col=other.species$kingdom.col,
#     xlim = c(-0.016, 0.016),
#     ylim = c(-0.03, 0.03),
#     xlab="", ylab=""
#   )
# }

# png(width = 550, height = 550, file ="plot_otherspecies.png")
# skew.plot(dataset = sim.run, scaling = paste0("Scaling = ", scaling))
# pic.saved <- dev.off()

# as_tibble(other.species) %>% 
#   filter(kingdom == "prokaryotes") %>% 
#   arrange(desc(AT.skew))

# skew.plot <- function(dataset, scaling){

#   # Plot the GC vs. AT skew plots

#   # Flag      Format     Description
#   # dataset  <Rdata>     Dataset of the equilibrium outputs from the
#   #                      numerically solved kinetic mutation rate
#   #                      equations.
#   # scaling  <numeric>   Scaling factor of the simulation results

#   smoothScatter(
#     y=dataset$GCskew,
#     x=dataset$ATskew,
#     nrpoints=100, nbin=1000,
#     xlab="AT skew", ylab="GC skew", main=scaling,
#     colramp= alpharamp(alpha=30),
#     cex.axis = 1.2, cex.lab = 1.1,
#     xlim = c(-RANGE, RANGE),
#     ylim = c(-RANGE, RANGE)
#   )

#   # plot(
#   #   x = dataset$ATskew, y = dataset$GCskew, col = "white",
#   #   xlim = c(-RANGE, RANGE),
#   #   ylim = c(-RANGE, RANGE)
#   # )

#   # rect(
#   #   xleft = -RANGE,
#   #   xright = RANGE,
#   #   ybottom = -RANGE,
#   #   ytop = RANGE,
#   #   lty = 3,
#   #   lwd = 6,
#   #   border = "white"
#   # )

#   rect(
#     xleft   = prok[prok$metadata == "AT_skew", "mean"]-
#               prok[prok$metadata == "AT_skew", "st.dev"],
#     xright  = prok[prok$metadata == "AT_skew", "mean"]+
#               prok[prok$metadata == "AT_skew", "st.dev"],
#     ybottom = prok[prok$metadata == "GC_skew", "mean"]-
#               prok[prok$metadata == "GC_skew", "st.dev"],
#     ytop    = prok[prok$metadata == "GC_skew", "mean"]+
#               prok[prok$metadata == "GC_skew", "st.dev"],
#     lwd = 3,
#     border = "black"
#   )

#   rect(
#     xleft   = euk[euk$metadata == "AT_skew", "mean"]-
#               euk[euk$metadata == "AT_skew", "st.dev"],
#     xright  = euk[euk$metadata == "AT_skew", "mean"]+
#               euk[euk$metadata == "AT_skew", "st.dev"],
#     ybottom = euk[euk$metadata == "GC_skew", "mean"]-
#               euk[euk$metadata == "GC_skew", "st.dev"],
#     ytop    = euk[euk$metadata == "GC_skew", "mean"]+
#               euk[euk$metadata == "GC_skew", "st.dev"],
#     lwd = 3,
#     border = "black"
#   )

#   rect(
#     xleft   = virus[virus$metadata == "AT_skew", "mean"]-
#               virus[virus$metadata == "AT_skew", "st.dev"],
#     xright  = virus[virus$metadata == "AT_skew", "mean"]+
#               virus[virus$metadata == "AT_skew", "st.dev"],
#     ybottom = virus[virus$metadata == "GC_skew", "mean"]-
#               virus[virus$metadata == "GC_skew", "st.dev"],
#     ytop    = virus[virus$metadata == "GC_skew", "mean"]+
#               virus[virus$metadata == "GC_skew", "st.dev"],
#     lwd = 3,
#     border = "black"
#   )

#   xlim = c(-RANGE, RANGE)
#   ylim = c(-RANGE, RANGE)
# }

# png(width = 550, height = 550, file ="plot_zoomed.png")
# # png(width = 550, height = 550, file ="plot.png")
# skew.plot(dataset = sim.run, scaling = paste0("Scaling = ", scaling))
# pic.saved <- dev.off()
# ########################################################################################