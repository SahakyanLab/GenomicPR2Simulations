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
  
  as.data.frame(atgc$out) %>%
    as_tibble() %>%
    mutate(across(.cols = c(Ca, Cg, Ct, Cc), ~.*100)) %>%
    ggplot(aes(x = time,
               y = Ca)) + 
    geom_line(aes(y = Ca), size = 1.2, col = "forestgreen") + 
    geom_line(aes(y = Cg), size = 1.2, col = "orange") + 
    geom_line(aes(y = Ct), size = 1.2, col = "red") + 
    geom_line(aes(y = Cc), size = 1.2, col = "blue") + 
    geom_vline(xintercept = atgc$Ch.time, linetype = "dashed") + 
    geom_text(data = data.frame(xpos = atgc$Ch.time, 
                                ypos =  min(ylim),
                                annotateText = "Chargaff eq. reached",
                                hjustvar = 0, vjustvar = 1.1), 
              aes(x = xpos, 
                  y = ypos, 
                  hjust = hjustvar, 
                  vjust = vjustvar, 
                  label = annotateText,
                  angle = 90),
              fontface = "bold", size = 3) + 
    geom_vline(xintercept = atgc$Eq.time) +
    geom_text(data = data.frame(xpos = atgc$Eq.time, 
                                ypos =  min(ylim),
                                annotateText = "Genome eq. reached",
                                hjustvar = 0, vjustvar = 1.1), 
              aes(x = xpos, 
                  y = ypos, 
                  hjust = hjustvar, 
                  vjust = vjustvar, 
                  label = annotateText,
                  angle = 90),
              fontface = "bold", size = 3) + 
    xlim(xlim) +
    ylim(ylim) +
    labs(x = "Time, byr",
         y = "Base content, %") + 
    geom_text(data = data.frame(
      xpos = 0,
      ypos = max(ylim),
      annotateText = paste(
        "t = 0 ",time.unit,"\n",
        "nG (orange) = ",round(atgc$inp$state["Cg"], 2),"\n",
        "nC (blue) = ",round(atgc$inp$state["Cc"], 2),"\n",
        "nA (green) = ",round(atgc$inp$state["Ca"], 2),"\n",
        "nT (red) = ",round(atgc$inp$state["Ct"], 2),"\n",
        "G+C = ",round(100*(atgc$inp$state["Cg"]+atgc$inp$state["Cc"])/
                         atgc$length.genome,2),"%","\n",
         "G/C = ",format(round(atgc$inp$state["Cg"]/atgc$inp$state["Cc"],2),nsmall=2),"\n",   
         "A/T = ",format(round(atgc$inp$state["Ca"]/atgc$inp$state["Ct"],2),nsmall=2), sep=""),
      hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, 
                  y = ypos, 
                  hjust = hjustvar, 
                  vjust = vjustvar, 
                  label = annotateText,
                  angle = 0), size = 2) + 
    geom_text(data = data.frame(
      xpos = max(xlim)-2,
      ypos = max(ylim), 
      annotateText = paste(
        "t = ",round(atgc$inp$span,2)," ",time.unit,"\n",
        "nG = ",round(as.vector(atgc$out[atgc$length.out,"Cg"]),2),"\n",
        "nC = ",round(as.vector(atgc$out[atgc$length.out,"Cc"]),2),"\n",
        "nA = ",round(as.vector(atgc$out[atgc$length.out,"Ca"]),2),"\n",
        "nT = ",round(as.vector(atgc$out[atgc$length.out,"Ct"]),2),"\n",
        "G+C = ",round(atgc$gc.content[atgc$length.out],2),"%","\n",
        "G/C = ",format(round(atgc$out[atgc$length.out,"Cg"]/atgc$out[atgc$length.out,"Cc"],2),nsmall=2),"\n",   
        "A/T = ",format(round(atgc$out[atgc$length.out,"Ca"]/atgc$out[atgc$length.out,"Ct"],2),nsmall=2), sep=""),
      hjustvar = 0, vjustvar = 1), 
              aes(x = xpos, 
                  y = ypos, 
                  hjust = hjustvar, 
                  vjust = vjustvar, 
                  label = annotateText,
                  angle = 0), size = 2) + 
    theme_bw()
}
################################################################################