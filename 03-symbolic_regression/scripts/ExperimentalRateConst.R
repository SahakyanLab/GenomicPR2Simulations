args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
TREK.SCALE <- as.logical(as.character(args[2]))
setwd(my_path)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
source("../lib/Equations.R")

#--------------------
# helper function

r.squared <- function(x, y){
  return(
    formatC(
      signif(summary(lm(x~y))$r.squared, digits=3),
      digits=3,
      format="fg",
      flag="#"
    )
  )
}

#--------------------

# load experimental rate constants
rate.const <- read.csv(
  file = paste0(
    "../../02-simulations/data/Raw/Michael_Lynch", 
    ifelse(TREK.SCALE, "/Trek_scale_", "/"), 
    "Lynch-2010-mutation-rates_FREQUENCY.csv"),
  header = TRUE
)

rate.const <- as_tibble(rate.const) %>%
    mutate(
        MUT = plyr::mapvalues(
            x = rate.const$MUT,
            from = c("AG|TC", "CT|GA", "AT|TA", "CA|GT", "AC|TG", "CG|GC"),
            to = c("kag|ktc", "kct|kga", "kat|kta", "kca|kgt", "kac|ktg", "kcg|kgc")
        )
    ) %>% 
    tidyr::separate(MUT, c("MUT_1", "MUT_2"))

# generate plots
rates <- c("kag", "ktc", "kct", "kga", "kat", "kta", "kca", "kgt", "kac", "ktg", "kcg", "kgc")
rate.const.columns <-  colnames(rate.const[3:length(rate.const)])

df <- rate.const %>% 
  tidyr::pivot_longer(1:2) %>% 
  dplyr::select(-name) %>% 
  dplyr::rename(rates = value) %>% 
  tidyr::pivot_longer(-rates) %>% 
  tidyr::pivot_wider(names_from = rates, values_from = value)

results <- lapply(1:length(rate.const.columns), function(cols){
  out <- sapply(rates, function(i){
    return(Equations(data = df[cols, ], var = i))
  })

  filtered.df <- as.data.frame(t(out)) %>% 
    dplyr::mutate(across(where(is.list), unlist)) %>% 
    dplyr::rename_with(~c("y_true", "y_pred")) %>% 
    dplyr::mutate(
      species = rate.const.columns[cols],
      rates = colnames(out)
    )

  return(filtered.df)
})

results <- do.call(rbind, results)
rownames(results) <- NULL

euk <- c(
  "H.sapiens", "D.melanogaster",	"C.elegans",	"A.thaliana",	
  "S.cerevisiae", "M.m.domesticus",	"D.pulex",	"P.pacificus",	
  "D.magna",	"P.troglodytes",	"A.nancymaae"
)
prok <- c(
  "E.coli",	"P.luminescens ATCC29999", "T.turnerae", 
  "M.smegmatis",	"P.fluorescens ATCC948",	"R. toruloides"
)

results <- results %>% 
  dplyr::mutate(
    kingdom = ifelse(species %in% euk, "eukaryotes", "prokaryotes"),
    kingdom.col = ifelse(kingdom == "eukaryotes", "purple", "darkgreen")
  )
  
# one plot with 12 rate constants for each species
p <- results %>% 
  group_by(species) %>% 
  mutate(R2 = paste0("R = ", r.squared(y_true, y_pred))) %>% 
  ggplot(aes(x = y_true, y = y_pred, col = rates)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  geom_point() + 
  facet_wrap(~species, ncol = 6) + 
  coord_cartesian(xlim = c(0, 1.7), ylim = c(0, 1.7)) +
  labs(x = "Y_True", y = "Y_Pred") + 
  theme_bw() + 
  theme(
    strip.text.x = element_text(size = 9),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
ggsave(
  filename = paste0(
    "../figures", 
    ifelse(TREK.SCALE, "/Trek_scale_", "/"), 
    "ExperimentalRateConstants_by_species.pdf"), 
  plot = p,
  height = 6, width = 11
)

# one plot with all species for each rate constant
p <- results %>% 
  group_by(rates) %>% 
  mutate(R2 = paste0("R = ", r.squared(y_true, y_pred))) %>% 
  ggplot(aes(x = y_true, y = y_pred, col = species)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  geom_point() + 
  facet_wrap(~rates, ncol = 6) + 
  # coord_cartesian(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1.5)) +
  labs(x = "Y_True", y = "Y_Pred") + 
  theme_bw() + 
  theme(
    strip.text.x = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
ggsave(
  filename = paste0(
    "../figures", 
    ifelse(TREK.SCALE, "/Trek_scale_", "/"), 
    "ExperimentalRateConstants_by_rates.pdf"), 
  plot = p,
  height = 4, width = 12
)

# col by kingdoms
p <- results %>% 
  group_by(rates) %>% 
  mutate(R2 = paste0("R = ", r.squared(y_true, y_pred))) %>% 
  ggplot(aes(x = y_true, y = y_pred)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") + 
  geom_point(colour = results$kingdom.col) + 
  facet_wrap(~rates, ncol = 6) + 
  geom_text(
    aes(x = 1.15, y = 0.05, label = R2),
    show.legend = FALSE
  ) + 
  coord_cartesian(xlim = c(0, 1.5), ylim = c(0, 1.5)) +
  labs(x = "Y_True", y = "Y_Pred") + 
  theme_bw() + 
  theme(
    strip.text.x = element_text(size = 12),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )
  
ggsave(
  filename = paste0(
    "../figures", 
    ifelse(TREK.SCALE, "/Trek_scale_", "/"), 
    "ExperimentalRateConstants_by_rates_bykingdomcol.pdf"), 
  plot = p,
  height = 4, width = 10
)
