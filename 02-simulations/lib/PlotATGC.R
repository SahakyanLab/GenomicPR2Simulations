################################################################################
PlotATGC <- function(atgc=atgc, time.unit="byr", xlim=c(0,5), ylim=c(0,60)){

  # Function to plot the results of the numerically solved kinetic rate 
  # equations for the simulation of the genome dynamics evolution

  # Dependencies
  #     SolveATGC.R, InputChecking.R, Simulation.R

  # Flag         Format       Description
  # atgc         <list>       Output of all the numerically solved kinetic 
  #                           mutation rate equations of one generation
  # time.unit    <character>  Time unit of the simulation.
  # xlim         <numeric>    X-axis limit for the plot of the simulation.
  # ylim         <numeric>    Y-axis limit for the plot of the simulation.

  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(dplyr))
  
  as_tibble(atgc$out) %>%
    mutate(across(.cols = c(Ca, Cg, Ct, Cc), ~.*100)) %>%
    ggplot(aes(x = time,
               y = Ca)) + 
    geom_line(aes(y = Ca), size = 1.2, col = "forestgreen") + 
    geom_line(aes(y = Cg), size = 1.2, col = "orange") + 
    geom_line(aes(y = Ct), size = 1.2, col = "red") + 
    geom_line(aes(y = Cc), size = 1.2, col = "blue") + 
    geom_vline(xintercept = atgc$Ch.time, linetype = "dashed") + 
    # geom_text(data = data.frame(xpos = atgc$Ch.time, 
    #                             ypos =  min(ylim),
    #                             annotateText = "Chargaff eq. reached",
    #                             hjustvar = 0, vjustvar = 1.1), 
    #           aes(x = xpos, 
    #               y = ypos, 
    #               hjust = hjustvar, 
    #               vjust = vjustvar, 
    #               label = annotateText,
    #               angle = 90),
    #           fontface = "bold", size = 3) + 
    # geom_vline(xintercept = atgc$Eq.time) +
    # geom_text(data = data.frame(xpos = atgc$Eq.time, 
    #                             ypos =  min(ylim),
    #                             annotateText = "Genome eq. reached",
    #                             hjustvar = 0, vjustvar = 1.1), 
    #           aes(x = xpos, 
    #               y = ypos, 
    #               hjust = hjustvar, 
    #               vjust = vjustvar, 
    #               label = annotateText,
    #               angle = 90),
    #           fontface = "bold", size = 3) + 
    coord_cartesian(xlim = xlim, ylim = ylim) + 
    labs(x = "Time, byr",
         y = "Base content, %") + 
    # geom_text(data = data.frame(
    #   xpos = 0,
    #   ypos = max(ylim),
    #   annotateText = paste(
    #     "t = 0 ",time.unit,"\n",
    #     "nG (orange) = ",round(atgc$inp$state["Cg"], 2),"\n",
    #     "nC (blue) = ",round(atgc$inp$state["Cc"], 2),"\n",
    #     "nA (green) = ",round(atgc$inp$state["Ca"], 2),"\n",
    #     "nT (red) = ",round(atgc$inp$state["Ct"], 2),"\n",
    #     "G+C = ",round(100*(atgc$inp$state["Cg"]+atgc$inp$state["Cc"])/
    #                      atgc$length.genome,2),"%","\n",
    #      "G/C = ",format(round(atgc$inp$state["Cg"]/atgc$inp$state["Cc"],2),nsmall=2),"\n",   
    #      "A/T = ",format(round(atgc$inp$state["Ca"]/atgc$inp$state["Ct"],2),nsmall=2), sep=""),
    #   hjustvar = 0, vjustvar = 1), 
    #           aes(x = xpos, 
    #               y = ypos, 
    #               hjust = hjustvar, 
    #               vjust = vjustvar, 
    #               label = annotateText,
    #               angle = 0), size = 2) + 
    # geom_text(data = data.frame(
    #   xpos = max(xlim)-1,
    #   ypos = max(ylim), 
    #   annotateText = paste(
    #     "t = ",round(atgc$inp$span,2)," ",time.unit,"\n",
    #     "nG = ",round(as.vector(atgc$out[atgc$length.out,"Cg"]),2),"\n",
    #     "nC = ",round(as.vector(atgc$out[atgc$length.out,"Cc"]),2),"\n",
    #     "nA = ",round(as.vector(atgc$out[atgc$length.out,"Ca"]),2),"\n",
    #     "nT = ",round(as.vector(atgc$out[atgc$length.out,"Ct"]),2),"\n",
    #     "G+C = ",round(atgc$gc.content[atgc$length.out],2),"%","\n",
    #     "G/C = ",format(round(atgc$out[atgc$length.out,"Cg"]/atgc$out[atgc$length.out,"Cc"],2),nsmall=2),"\n",   
    #     "A/T = ",format(round(atgc$out[atgc$length.out,"Ca"]/atgc$out[atgc$length.out,"Ct"],2),nsmall=2), sep=""),
    #   hjustvar = 0, vjustvar = 1), 
    #           aes(x = xpos, 
    #               y = ypos, 
    #               hjust = hjustvar, 
    #               vjust = vjustvar, 
    #               label = annotateText,
    #               angle = 0), size = 2) + 
    theme_bw()
}
################################################################################
setwd("/Users/paddy/Documents/DPhil/01-Chargaff/02-simulations/lib/")

# Load required supplementary functions and packages
suppressPackageStartupMessages(library(deSolve))
suppressPackageStartupMessages(library(dplyr))
source("SolveATGC.R")
source("States.R")

prokaryotes.df <- read.csv(file = "../../01-genome_composition/data/01-Prokaryotes/All/all_filtered_dataframe.csv", 
                           header=TRUE)
eukaryotes.df  <- read.csv(file = "../../01-genome_composition/data/02-Eukaryotes/All/all_filtered_dataframe.csv", 
                           header=TRUE)

# Import calculated mutation rates from trek paper
note.one   <- read.csv("../data/Raw/Trek-paper-Note-1-mutation-rates.csv", 
                       header = TRUE)
note.two   <- read.csv("../data/Raw/Trek-paper-Note-2-mutation-rates.csv", 
                       header = TRUE)
note.three <- read.csv("../data/Raw/Trek-paper-Note-3-mutation-rates.csv", 
                       header = TRUE)

Acont      = 0.278
Gcont      = 0.167
Ccont      = 0.278
span       = 4.28
step       = 0.001
max.runs   = 1
muttype    = "Non_Symmetric"
dist       = "uniform"
species    = "NONE"
scale.fac  = 1
tolerance  = FALSE
tol.return = "NONE"
sim.evol   = FALSE 
sy.reg.run = FALSE
NCPU       = 1
seed       = 1 

all.states <- States(iterations = max.runs, seed = seed, NCPU = NCPU)
state <- all.states[1,]

rates <- runif(n=12, min=0, max=note.three$MAXMEAN+note.three$MAXSD*scale.fac)
kag = rates[1]
kat = rates[2]
kac = rates[3]
kga = rates[4]
kgt = rates[5]
kgc = rates[6]
kta = rates[7]
ktg = rates[8]
ktc = rates[9]
kca = rates[10]
kcg = rates[11]
kct = rates[12]

parameters <- c(kag=kag, kat=kat, kac=kac,
                kga=kga, kgt=kgt, kgc=kgc,
                kta=kta, ktg=ktg, ktc=ktc,
                kca=kca, kcg=kcg, kct=kct) #mut/byr

atgc <- SolveATGC(parameters = parameters, state = state, 
                  step = step, span = span,
                  EQtolerance = FALSE, CHtolerance = FALSE)

p = PlotATGC(atgc = atgc, xlim = c(0, 4.5))
p

ggsave(plot = p, filename = "uniform_NS.pdf", width = 5, height = 4)
