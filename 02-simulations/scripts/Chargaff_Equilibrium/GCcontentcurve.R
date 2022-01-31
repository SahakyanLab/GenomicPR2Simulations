args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save.as <- as.character(args[2])
setwd(my_path)

# Load required supplementary functions and packages
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))

# ----------------
# define helper functions
theoretical <- function(x){
  return(x/(1-x))
}

hist.count <- function(data){
  p <- hist(data[, "GC"], plot = FALSE)
  return(p$counts)
}

# ----------------

# load data 
prokaryotes.df <- read.csv(file = "../../../01-genome_composition/data/01-Prokaryotes/All/all_filtered_dataframe.csv", 
                           header=TRUE)
eukaryotes.df  <- read.csv(file = "../../../01-genome_composition/data/02-Eukaryotes/All/all_filtered_dataframe.csv", 
                           header=TRUE)
viruses.df     <- read.csv(file = "../../../01-genome_composition/data/03-Viruses//All/all_filtered_dataframe.csv", 
                           header=TRUE)
species.gc     <- read.csv(file = "../../data/Raw/Michael_Lynch/GC_vs_Rates.csv",
                           header = TRUE)

euk <- c(
  "H.sapiens", "D.melanogaster",	"C.elegans",	"A.thaliana",	
  "M.m.domesticus",	"D.pulex",	"P.pacificus",	
  "D.magna",	"P.troglodytes",	"A.nancymaae"
)

prok <- c(
  "E.coli",	"P.luminescens ATCC29999", "T.turnerae", "S.cerevisiae", 
  "M.smegmatis",	"P.fluorescens ATCC948",	"R. toruloides"
)

species.gc <- species.gc %>% 
  dplyr::mutate(
    kingdom = ifelse(Species %in% euk, "eukaryotes", "prokaryotes"),
    kingdom.col = ifelse(kingdom == "eukaryotes", "purple", "darkgreen")
  ) %>% 
  arrange(desc(kingdom))

prok.hist <- data.frame(GC = prokaryotes.df$G_plus_C)
euk.hist <- data.frame(GC = eukaryotes.df$G_plus_C)
vir.hist <- data.frame(GC = viruses.df$G_plus_C)

# theoretical gc curve
theo.line <- theoretical(x = seq(from = 0, to = 1, by = 0.001))
theo.line <- theo.line[is.finite(theo.line)]
time <- seq(from = 0, to = 100, length.out = length(theo.line))

theo.gc.line <- data.frame(
     Species = "Theoretical",
     GC.average = time, 
     Rates = theo.line,
     kingdom = "Theoretical",
     kingdom.col = "darkred"
)

df <- rbind(theo.gc.line, species.gc)

# euclidean distance from point to tangent
spl <- smooth.spline(theo.gc.line$GC.average, theo.gc.line$Rates, spar=0.3)
newx <- seq(min(theo.gc.line$GC.average), max(theo.gc.line$Rates), 0.1)
pred <- predict(spl, x=newx, deriv=0)

eucl.dist <- sapply(1:nrow(species.gc), function(x){
  newx  <- species.gc[x, "GC.average"]
  pred0 <- predict(spl, x=newx, deriv=0)
  pred1 <- predict(spl, x=newx, deriv=1)
  yint  <- pred0$y - (pred1$y*newx)
  xint  <- -yint/pred1$y

  tangent <- data.frame(
    x = df$GC.average,
    y = pred1$y*df$GC.average + yint) %>% 
    dplyr::arrange(x)

  distance <- tangent[match(newx, tangent$x), "y"]-species.gc[x, "Rates"]

  return(distance)
})

eucl.df <- cbind(species.gc, eucl.dist) 

p.den <- eucl.df %>% 
  ggplot(aes(x = eucl.dist, fill = kingdom)) + 
  geom_density(alpha = 0.5) +  
  scale_fill_manual(
    values = unique(eucl.df$kingdom.col),
    limits = unique(eucl.df$kingdom)
  ) + 
  theme_bw() + 
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) + 
  labs(
    x = "Euclidean Distance",
    y = "Density"
  )

ggsave(
  filename = paste0("../../figures/Chargaff_Equilibrium/Theoretical_GC/Euclidean_Dist.", save.as),
  plot = p.den,
  height = 7, width = 10
)

# label each species
p <- df %>% 
  ggplot() + 
  geom_line(
      data = subset(df, Species == "Theoretical"),
      aes(x = GC.average, y = Rates),
      colour = df$kingdom.col[df$Species == "Theoretical"],
      size = 1.2,
      show.legend = FALSE
  ) + 
  geom_point(
      data = subset(df, Species != "Theoretical"),
      aes(x = GC.average, y = Rates),
      colour = df$kingdom.col[df$Species != "Theoretical"],
  ) + 
  theme_bw() + 
  theme(
      axis.text.x = element_text(size = 15),
      axis.text.y = element_text(size = 15),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
  ) + 
  labs(
    x = "GC content (%)",
    y = "(n+j)/(i+m)"
  )

p1 <- p + ggrepel::geom_text_repel(
       data = subset(df, Species != "Theoretical"),
       aes(x = GC.average, y = Rates, label = Species),
       box.padding = 1, 
       max.overlaps = Inf,
       min.segment.length = unit(0.1, "lines"),
       size = 3
  ) + 
  coord_cartesian(
    xlim = c(0, 100), 
    ylim = c(-1.5, 12)
  )

ggsave(
     filename = paste0("../../figures/Chargaff_Equilibrium/Theoretical_GC/GCcontenthist_otherspecies_withlabels.", save.as),
     plot = p1,
     height = 8, width = 10
)

p2 <- p + 
  coord_cartesian(
    xlim = c(0, 100), 
    ylim = c(0, 12)
  )

ggsave(
     filename = paste0("../../figures/Chargaff_Equilibrium/Theoretical_GC/GCcontenthist_otherspecies.", save.as),
     plot = p2,
     height = 7, width = 10
)