suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(broom))

args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
ncpu <- as.numeric(args[3])
species.name <- as.character(args[4])
scaling = 1

setwd(my_path)

load("../../data/simulation/Non_Symmetric-uniform-scaling-1.Rdata")
sim.run$nC <- (sim.run$GC/100)/(1+sim.run$GCratio)
sim.run$nG <- (sim.run$GC/100)-sim.run$nC
sim.run$nT <- (1-(sim.run$GC/100)) / (1+sim.run$ATratio)
sim.run$nA <- (1-(sim.run$GC/100)) - sim.run$nT
sim.run$A_minus_T <- sim.run$nA-sim.run$nT
sim.run$G_minus_C <- sim.run$nG-sim.run$nC

#----------------------
# obtain k-values from rhombus plot that fall within the 
# "tolerance" values obtained from experimental values

filtered.df <- function(dataset, species, analysis = "compliance"){

  # Apply the PR-2 tolerance to the compliance region of the datasets

  # Flag      Format     Description
  # dataset  <Rdata>     Dataset of the equilibrium outputs from the
  #                      numerically solved kinetic mutation rate
  #                      equations.
  # species  <character> Apply the PR-2 tolerance zone of 
  #                      "Prokaryotes", "Eukaryotes" or "Viruses" organisms
  # analysis <character> If "compliance", apply the PR-2 tolerance zone of the species

  if(species == "prokaryotes"){
    file.name <- "../../../01-genome_composition/data/01-Prokaryotes/PR_compliance/PR2_fluctuations.csv"
  } else if(species == "eukaryotes"){
    file.name <- "../../../01-genome_composition/data/02-Eukaryotes/data/PR_compliance/PR2_fluctuations.csv"
  } else if(species == "viruses"){
    file.name <- "../../../01-genome_composition/data/03-Viruses/data/PR_compliance/PR2_fluctuations.csv"
  }
  
  if(file.exists(file.name)){
    # import chargaff compliancy region
    fluc.species <- read.csv(file = file.name, header = TRUE)

    if(analysis == "pca"){
      y <- dataset %>%
        as_tibble() %>%
        mutate(compliance = 
                ifelse(((ATskew >= fluc.species[fluc.species$metadata == "AT_skew","mean"]-
                    fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                    (ATskew <= fluc.species[fluc.species$metadata == "AT_skew","mean"]+
                      fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                    (GCskew >= fluc.species[fluc.species$metadata == "GC_skew","mean"]-
                      fluc.species[fluc.species$metadata == "GC_skew","st.dev"]) & 
                    (GCskew <= fluc.species[fluc.species$metadata == "GC_skew","mean"]+
                      fluc.species[fluc.species$metadata == "GC_skew","st.dev"])),
                    "Compliant", "Not Compliant")) %>%
        relocate(compliance, .after = kct)
    } else {
      y <- dataset %>%
        as_tibble() %>%
        filter(ATskew >= fluc.species[fluc.species$metadata == "AT_skew","mean"]-
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
print(paste0("Performing PCA on ", species.name, "..."))
sim.run.rates <- filtered.df(dataset = sim.run, 
                            species = species.name, 
                            analysis = "pca") %>%
  select(7:(dim(sim.run)[2]-5))

df.compliant <- sim.run.rates %>% 
  filter(compliance == "Compliant")

df.not.compliant <- sim.run.rates %>% 
  filter(compliance == "Not Compliant") %>% 
  slice_sample(n = 200000)
  
sim.run.rates <- rbind(df.compliant, df.not.compliant)

rates.pca.fit <- sim.run.rates %>%
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)

summary(rates.pca.fit)

pca.df <- data.frame(pc_1 = rates.pca.fit$x[,1],
                     pc_2 = rates.pca.fit$x[,2],
                     pc_3 = rates.pca.fit$x[,3],
                     compliance = sim.run.rates$compliance)

pca.df.compliant <- pca.df %>%
  filter(compliance == "Compliant")

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

ggsave(filename = paste0("../../figures/PCA/",species.name,"_PCA.jpeg"),
       plot = pca.plot, dpi = 150)

# Variance explained by each principal component
pc_eigenvalues <- rates.pca.fit$sdev^2
pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)),
                         variance = pc_eigenvalues) %>%
  mutate(pct = variance/sum(variance)*100) %>%
  mutate(pct_cum = cumsum(pct))

print(paste0("Obtaining eigenvalues of all PCs for ", species.name, "..."))
# Eigenvalues of all PCs
pc_eigenvalues.plot <- pc_eigenvalues %>%
  ggplot(aes(x = PC, y = variance)) + 
  geom_line(aes(y = variance, group = 1)) +
  geom_point(aes(y = variance)) + 
  labs(x = "Principal component",
       y = "Variances",
       title = "Eigenvalues of all principal components")

ggsave(filename = paste0("../../figures/PCA/", species.name,"_Eigenvalues.png"),
       plot = pc_eigenvalues.plot, dpi = 300)
  
print(paste0("Obtaining cumulative variance of all PCs for ", species.name, "..."))
# Cumulative variance of all PCs
pc_eigenvalues.plot <- pc_eigenvalues %>%
  ggplot(aes(x = PC)) +
  geom_col(aes(y = pct)) +
  geom_line(aes(y = pct_cum, group = 1)) + 
  geom_point(aes(y = pct_cum)) + 
  labs(x = "Principal component",
       y = "Fraction variance explained",
       title = "Screeplot of cumulative variance of all PCs")

ggsave(filename = paste0("../../figures/PCA/", species.name,"_Cumulative_variance.png"),
       plot = pc_eigenvalues.plot, dpi = 300)
