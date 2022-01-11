args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save.as <- as.character(args[2])
setwd(my_path)

# load dependencies
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(geneplotter))

colfun = colorRampPalette(c("white","blue","skyblue",
                            "chartreuse3","green","yellow",
                            "orange","red","darkred"))

# Import data frames of the three species types
prokaryotes.df <- read.csv(file = "../../data/01-Prokaryotes/All/all_filtered_dataframe.csv", header=TRUE)
eukaryotes.df  <- read.csv(file = "../../data/02-Eukaryotes/All/all_filtered_dataframe.csv", header=TRUE)
viruses.df     <- read.csv(file = "../../data/03-Viruses/All/all_filtered_dataframe.csv", header=TRUE)

#-----------------------------
# Create data frame and save plots 
# G+C content 
gc <- c(eukaryotes.df$G_plus_C*100, 
        prokaryotes.df$G_plus_C*100, 
        viruses.df$G_plus_C*100)

min(viruses.df$G_plus_C*100)
max(viruses.df$G_plus_C*100)

GC_hist <- function(dataset, title = ""){

  # Plot the GC content as histograms

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # title    <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  Plot <- dataset %>%
    as_tibble() %>%
    ggplot(aes(x = G_plus_C*100, 
               y = ..density..)) + 
    geom_histogram(fill = "skyblue",
                   color = "black",
                   alpha = 1,
                   breaks = seq(min(dataset$G_plus_C)*100,
                                max(dataset$G_plus_C)*100,
                                length.out = 51)) + 
    coord_cartesian(xlim = c(min(gc), max(gc))) + 
    labs(x = "",
         y = "",
         title = paste0(title, "\n", 
                        "Mean = ", format(round(mean(dataset$G_plus_C)*100, 2), 
                                          nsmall = 2), " ",
                        "St.Dev = ", format(round(sd(dataset$G_plus_C)*100, 2), 
                                            nsmall = 2)))
  return(Plot)
}

grid.arrange(GC_hist(eukaryotes.df, "Eukaryotes"),
             GC_hist(prokaryotes.df, "Prokaryotes"),
             GC_hist(viruses.df, "Viruses"),
             ncol = 3, nrow = 1) %>%
  ggsave(width=15, height=8, filename = paste0("../../figures/00-All_Species/GC_hist.", save.as))
print("GC histogram plots done!", quote = FALSE)

#-----------------------------
# A+T content 
AT_hist <- function(dataset, title = ""){

  # Plot the AT content as histograms

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # title    <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  Plot <- dataset %>%
    as_tibble() %>%
    ggplot(aes(x = A_plus_T*100, 
               y = ..density..)) + 
    geom_histogram(fill = "skyblue",
                   color = "black",
                   alpha = 1,
                   breaks = seq(min(dataset$A_plus_T)*100,
                                max(dataset$A_plus_T)*100,
                                length.out = 51)) + 
    labs(x = paste0("A+T content (%)", "\n",
                    "Mean = ", format(round(mean(dataset$A_plus_T)*100, 2), 
                                      nsmall = 2), " ",
                    "St.Dev = ", format(round(sd(dataset$A_plus_T)*100, 2), 
                                        nsmall = 2)),
         y = "Density",
         title = title)
  return(Plot)
}

grid.arrange(GC_hist(eukaryotes.df, "Eukaryotes"),
             GC_hist(prokaryotes.df, "Prokaryotes"),
             GC_hist(viruses.df, "Viruses"),
             ncol = 3, nrow = 1) %>%
  ggsave(width=15, height=8, filename = paste0("../../figures/00-All_Species/AT_hist.", save.as))
print("AT histogram plots done!", quote = FALSE)
#-----------------------------
# G-C vs. A-T 
difference.plots <- function(dataset, species){

  # Plot the G-C vs. A-T difference 

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # species  <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  x = dataset$A_minus_T
  y = dataset$G_minus_C
  
  smoothScatter(y=y, x=x,
                nrpoints=100, nbin=1000,
                bandwidth=c(diff(range(x))/500, diff(range(y))/500),
                xlab="A-T", ylab="G-C", main=species,
                colramp=colfun,
                cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  abline(v = 0, lty = 2, lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=15, height=8, file="../../figures/00-All_Species/G-C_vs_A-T.pdf")
} else {
  png(width=1200, height=500, file="../../figures/00-All_Species/G-C_vs_A-T.png")
}
par(mfrow=c(1,3))
difference.plots(eukaryotes.df, "Eukaryotes")
difference.plots(prokaryotes.df, "Prokaryotes")
difference.plots(viruses.df, "Viruses")
plot.save <- dev.off()
print("G-C vs. A-T plots done!", quote = FALSE)

#-----------------------------
# GC-ratio vs. AT-ratio
ratio.plot <- function(dataset, species){

  # Plot the GC vs. AT ratios

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # species  <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  smoothScatter(y=dataset$GC_ratio,
                x=dataset$AT_ratio,
                nrpoints=100, nbin=1000,
                xlab="AT ratio", ylab="GC-ratio", main=species,
                colramp=colfun,
                cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  abline(v = 1, lty = 2, lwd = 2)
  abline(h = 1, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=15, height=8, file="../../figures/00-All_Species/GC-ratio_vs_AT-ratio.pdf")
} else {
  png(width=1200, height=500, file="../../figures/00-All_Species/GC-ratio_vs_AT-ratio.png")
}
par(mfrow=c(1,3))
ratio.plot(eukaryotes.df, "Eukaryotes")
ratio.plot(prokaryotes.df, "Prokaryotes")
ratio.plot(viruses.df, "Viruses")
plot.save <- dev.off()
print("GC vs. AT ratio plots done!", quote = FALSE)

#-----------------------------
# GC-skew vs. AT-skew
skew.plot <- function(dataset, species){

  # Plot the GC vs. AT skews

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # species  <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  smoothScatter(y=dataset$GC_skew, 
                x=dataset$AT_skew,
                nrpoints=100, nbin=1000,
                xlab="AT skew", ylab="GC skew", main=species,
                colramp=colfun,
                cex.lab=1.5, cex.axis=1.5, cex.main=1.5)
  abline(v = 0, lty = 2, lwd = 2)
  abline(h = 0, lty = 2, lwd = 2)
}

if(save.as == "pdf"){
  pdf(width=15, height=8, file="../../figures/00-All_Species/GC-skew_vs_AT-skew.pdf")
} else {
  png(width=1200, height=500, file="../../figures/00-All_Species/GC-skew_vs_AT-skew.png")
}
par(mfrow=c(1,3))
skew.plot(eukaryotes.df, "Eukaryotes")
skew.plot(prokaryotes.df, "Prokaryotes")
skew.plot(viruses.df, "Viruses")
plot.save <- dev.off()
print("GC vs. AT skew done!", quote = FALSE)

#-----------------------------
# GC-skew
GC_skew <- function(dataset, species){

  # Plot the GC skew as a histogram

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # species  <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  Plot <- dataset %>%
    as_tibble() %>%
    ggplot(aes(x = GC_skew, 
               y = ..density..)) + 
    geom_histogram(fill = "skyblue",
                   color = "black",
                   alpha = 1,
                   breaks = seq(min(dataset$GC_skew),
                                max(dataset$GC_skew),
                                length.out = 51)) + 
    labs(x = paste0("GC skew (%)", "\n",
                    "Mean = ", format(round(mean(dataset$GC_skew)*100, 2), 
                                      nsmall = 2), " ",
                    "St.Dev = ", format(round(sd(dataset$GC_skew)*100, 2), 
                                        nsmall = 2)),
         y = "Density",
         title = species)
  return(Plot)
}

grid.arrange(GC_skew(eukaryotes.df, "Eukaryotes"),
             GC_skew(prokaryotes.df, "Prokaryotes"),
             GC_skew(viruses.df, "Viruses"),
             ncol = 3, nrow = 1) %>%
  ggsave(width=15, height=8, filename = paste0("../../figures/00-All_Species/GC_skew.", save.as))
print("GC skew plots done!", quote = FALSE)

#-----------------------------
# AT-skew
AT_skew <- function(dataset, species){

  # Plot the AT skew as a histogram

  # Flag      Format      Description
  # dataset  <Rdata>      Dataset of the equilibrium outputs from the
  #                       numerically solved kinetic mutation rate
  #                       equations.
  # species  <character>  Name of the species for which the PR-2
  #                       compliance is used to apply to the dataset

  Plot <- dataset %>%
    as_tibble() %>%
    ggplot(aes(x = AT_skew, 
               y = ..density..)) + 
    geom_histogram(fill = "skyblue",
                   color = "black",
                   alpha = 1,
                   breaks = seq(min(dataset$AT_skew),
                                max(dataset$AT_skew),
                                length.out = 51)) + 
    labs(x = paste0("AT skew (%)", "\n",
                    "Mean = ", format(round(mean(dataset$AT_skew)*100, 2), 
                                      nsmall = 2), " ",
                    "St.Dev = ", format(round(sd(dataset$AT_skew)*100, 2), 
                                        nsmall = 2)),
         y = "Density",
         title = species)
  return(Plot)
}

grid.arrange(AT_skew(eukaryotes.df, "Eukaryotes"),
             AT_skew(prokaryotes.df, "Prokaryotes"),
             AT_skew(viruses.df, "Viruses"),
             ncol = 3, nrow = 1) %>%
  ggsave(width=15, height=8, filename = paste0("../../figures/00-All_Species/AT_skew.", save.as))
print("AT skew plots done!", quote = FALSE)

#-----------------------------
# Base content compliance
# GC content
r.squared <- function(x, y){summary(lm(x~y))$r.squared}

rsq.plots <- function(dataset, 
                      species, 
                      content = "GC", 
                      title = ""){

  # Plot the GC vs. AT ratios

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # species  <character> Apply the PR-2 tolerance zone of 
  #                      "Prokaryotes", "Eukaryotes" or "Viruses" organisms
  # content  <character> Plot the G vs. C content or A vs. T content
  # title    <character> Name of the species for which the PR-2
  #                      compliance is used to apply to the dataset

  dataset <- dataset %>%
    as_tibble() %>%
    dplyr::select(G_content, 
                  C_content,
                  A_content,
                  T_content) %>%
    dplyr::mutate(across(.cols = everything(), ~.*100))
  
  ColMax <- max(dataset)
  
  if(content == "GC"){
    rsq <- r.squared(dataset$G_content, dataset$C_content)
    
    Plot <- dataset %>%
      as_tibble() %>%
      ggplot(aes(x = G_content,
                 y = C_content)) + 
      geom_point() + 
      geom_smooth(method = "lm", formula = y ~ x) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      coord_cartesian(xlim = c(0, ColMax),
                      ylim = c(0, ColMax)) + 
      labs(x = "",
           y = "",
           title = paste0(title, "\n",
                          "R = ", format(round(rsq, 3), 
                                         nsmall = 3), "\n",
                          "Nr. of Points = ", dim(dataset)[1]))
    
  } else {
    rsq <- r.squared(dataset$A_content, dataset$T_content)
    
    Plot <- dataset %>%
      as_tibble() %>%
      ggplot(aes(x = A_content,
                 y = T_content)) + 
      geom_point() + 
      geom_smooth(method = "lm", formula = y ~ x) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
      coord_cartesian(xlim = c(0, ColMax),
                      ylim = c(0, ColMax)) + 
      labs(x = "",
           y = "",
           title = paste0(title, "\n",
                          "R = ", format(round(rsq, 3), 
                                         nsmall = 3), "\n",
                          "Nr. of Points = ", dim(dataset)[1]))
  }
  return(Plot)
}

euk.p1 <- rsq.plots(eukaryotes.df, "Eukaryotes", "GC", "Eukaryotes")
euk.p2 <- rsq.plots(eukaryotes.df, "Eukaryotes", "AT")

prok.p1 <- rsq.plots(prokaryotes.df, "Prokaryotes", "GC", "Prokaryotes")
prok.p2 <- rsq.plots(prokaryotes.df, "Prokaryotes", "AT")

vir.p1 <- rsq.plots(viruses.df, "Viruses", "GC", "Viruses")
vir.p2 <- rsq.plots(viruses.df, "Viruses", "AT")

grid.arrange(euk.p1, prok.p1, vir.p1,
             euk.p2, prok.p2, vir.p2,
             ncol = 3, nrow = 2) %>%
  ggsave(width=10, height=7, 
         filename = paste0("../../figures/00-All_Species/content_comparisons.", save.as))
print("Base composition plots for all species done!", quote = FALSE)