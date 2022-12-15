args <- commandArgs(trailingOnly = TRUE)
my.path <- as.character(args[1])
species.name <- as.character(args[2])
setwd(my.path)

# load dependencies
suppressPackageStartupMessages(suppressWarnings(library(dplyr)))
suppressPackageStartupMessages(suppressWarnings(library(ggplot2)))
suppressPackageStartupMessages(suppressWarnings(library(tidyr)))
suppressPackageStartupMessages(suppressWarnings(library(broom)))

sim.run <- readRDS(file = paste0("../../02-simulations/data/Main_Simulation/", 
                                 "Non_Symmetric-uniform/Non_Symmetric-uniform", 
                                 "-scaling-1.Rdata"))
sim.run$nC <- (sim.run$GC/100)/(1+sim.run$GCratio)
sim.run$nG <- (sim.run$GC/100)-sim.run$nC
sim.run$nT <- (1-(sim.run$GC/100)) / (1+sim.run$ATratio)
sim.run$nA <- (1-(sim.run$GC/100)) - sim.run$nT
sim.run$A_minus_T <- sim.run$nA-sim.run$nT
sim.run$G_minus_C <- sim.run$nG-sim.run$nC

#----------------------
# obtain k-values from rhombus plot that fall within the 
# "tolerance" values obtained from experimental values

#' @description
#' Apply the PR-2 tolerance to the compliance region of the datasets
#' @param dataset Rdata. Dataset of the equilibrium outputs from the numerically
#'  solved kinetic mutation rate equations.
#' @param species Character vector. Apply the PR-2 tolerance zone of
#'  c("prokaryotes", "eukaryotes", "viruses") kingdoms.
#' @param analysis Character vector. 
#'  "compliance": apply PR-2 tolerance zone of the species
#'  "pca": for principal component analysis
filtered.df <- function(dataset, species, analysis = "compliance"){
  file.name <- switch(species,
    "prokaryotes" = "../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv",
    "eukaryotes" = "../../01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv",
    "viruses" = "../../01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv"
  )
  
  if(file.exists(file.name)){
    # import chargaff compliancy region
    fluc.species <- read.csv(file = file.name, header = TRUE)

    if(analysis == "pca"){
      y <- dataset %>%
        as_tibble() %>%
        dplyr::mutate(compliance = 
                ifelse(((ATskew >= fluc.species[fluc.species$metadata == "AT_skew","mean"]-
                    fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                    (ATskew <= fluc.species[fluc.species$metadata == "AT_skew","mean"]+
                      fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                    (GCskew >= fluc.species[fluc.species$metadata == "GC_skew","mean"]-
                      fluc.species[fluc.species$metadata == "GC_skew","st.dev"]) & 
                    (GCskew <= fluc.species[fluc.species$metadata == "GC_skew","mean"]+
                      fluc.species[fluc.species$metadata == "GC_skew","st.dev"])),
                    "Compliant", "Not Compliant"),
                    .after = kct)
    } else {
      y <- dataset %>%
        as_tibble() %>%
        dplyr::filter(ATskew >= fluc.species[fluc.species$metadata == "AT_skew","mean"]-
                fluc.species[fluc.species$metadata == "AT_skew","st.dev"] & 
                ATskew <= fluc.species[fluc.species$metadata == "AT_skew","mean"]+
                fluc.species[fluc.species$metadata == "AT_skew","st.dev"] & 
                GCskew >= fluc.species[fluc.species$metadata == "GC_skew","mean"]-
                fluc.species[fluc.species$metadata == "GC_skew","st.dev"] & 
                GCskew <= fluc.species[fluc.species$metadata == "GC_skew","mean"]+
                fluc.species[fluc.species$metadata == "GC_skew","st.dev"])
    }
  } else {
    stop("File does not exist!")
  }
  return(y)
}

#----------------------
# Principal Component Analysis
species.name = "prokaryotes"

print(paste0("Performing PCA on ", species.name, "..."))
sim.run.rates <- filtered.df(dataset = sim.run, 
                            species = species.name, 
                            analysis = "pca") %>%
  dplyr::select(7:(dim(sim.run)[2]-5))

df.compliant <- sim.run.rates %>% 
  dplyr::filter(compliance == "Compliant")

df.not.compliant <- sim.run.rates %>% 
  dplyr::filter(compliance == "Not Compliant") %>% 
  dplyr::slice_sample(n = 200000)
  
sim.run.rates <- rbind(df.compliant, df.not.compliant)
compliance.col <- dplyr::pull(sim.run.rates["compliance"])
cols.names <- colnames(sim.run.rates)

# normalise
norm.sim.run.rates <- apply(
  dplyr::select(sim.run.rates, -compliance),
  2, scale, scale = TRUE, center = TRUE
)
colnames(norm.sim.run.rates) <- cols.names[1:(length(cols.names)-1)]
rates.pca.fit <- prcomp(
  x = norm.sim.run.rates, 
  center = FALSE, 
  scale. = FALSE
)
pca.df <- data.frame(pc_1 = rates.pca.fit$x[,1],
                     pc_2 = rates.pca.fit$x[,2],
                     pc_3 = rates.pca.fit$x[,3],
                     compliance = compliance.col)

pca.df.compliant <- pca.df %>%
  dplyr::filter(compliance == "Compliant")

# 2-Dimensional plots
pca.plot <- pca.df %>%
  ggplot(aes(x = pc_1, 
             y = pc_2,
             color = compliance)) + 
  geom_point(alpha = 0.03) +
  geom_point(data = pca.df.compliant,
             aes(x = pc_1, 
                 y = pc_2)) + 
  labs(x = "PC1",
       y = "PC2")

dir.create(
  path = "../figures/PCA/",
  showWarnings = FALSE,
  recursive = TRUE
)
ggsave(
  filename = paste0("../figures/PCA/",species.name,"_PCA.png"),
  plot = pca.plot, 
  dpi = 150
)

# Variance explained by each principal component
pc_eigenvalues <- rates.pca.fit$sdev^2
pc_eigenvalues <- tibble(
  PC = factor(1:length(pc_eigenvalues)),
  variance = pc_eigenvalues) %>%
  mutate(
    pct = variance/sum(variance)*100,
    pct_cum = cumsum(pct)
  )

# Eigenvalues of all PCs
pc_eigenvalues.plot <- pc_eigenvalues %>%
  ggplot(aes(x = PC, y = variance)) + 
  geom_line(aes(y = variance, group = 1)) +
  geom_point(aes(y = variance)) + 
  labs(x = "Principal component",
       y = "Variances",
       title = "Eigenvalues of all principal components")

ggsave(
  filename = paste0("../figures/PCA/", species.name,"_Eigenvalues.png"),
  plot = pc_eigenvalues.plot, 
  dpi = 300
)
  
# Cumulative variance of all PCs
pc_eigenvalues.plot <- pc_eigenvalues %>%
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) + 
  labs(x = "Principal component",
       y = "Fraction variance explained",
       title = "Screeplot of cumulative variance of all PCs")

ggsave(
  filename = paste0("../figures/PCA/", species.name,"_Cumulative_variance.png"),
  plot = pc_eigenvalues.plot, 
  dpi = 300
)