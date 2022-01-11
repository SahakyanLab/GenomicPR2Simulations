args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save.as <- as.character(args[2])
species <- as.character(args[3])
setwd(my_path)

# load dependencies
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(geneplotter))

colfun = colorRampPalette(c("white","blue","skyblue",
                            "chartreuse3","green","yellow",
                            "orange","red","darkred"))

# Import data frame from downloaded data
filtered.df <- read.csv(file = paste0("../data/", species, "/All/all_filtered_dataframe.csv"), header=TRUE)

#-----------------------------
# Mean/St.dev values for fluctuations from 2nd parity rule 
mean.sd.df = data.frame("metadata" = c("A_minus_T", "G_minus_C",
                                       "AT_skew", "GC_skew",
                                       "AT_ratio", "GC_ratio"),
                        "mean"     = c(mean(filtered.df$A_minus_T),
                                       mean(filtered.df$G_minus_C),
                                       mean(filtered.df$AT_skew),
                                       mean(filtered.df$GC_skew),
                                       mean(filtered.df$AT_ratio),
                                       mean(filtered.df$GC_ratio)),
                        "st.dev"   = c(sd(filtered.df$A_minus_T),
                                       sd(filtered.df$G_minus_C),
                                       sd(filtered.df$AT_skew),
                                       sd(filtered.df$GC_skew),
                                       sd(filtered.df$AT_ratio),
                                       sd(filtered.df$GC_ratio)))

write.csv(mean.sd.df, file = paste0("../data/", species, "/All/PR2_compliance/PR2_fluctuations.csv"), row.names = FALSE)

#-----------------------------
# Create data frame and save plots 
# G+C content 
Plot <- filtered.df %>%
        as_tibble() %>%
        ggplot(aes(x = G_plus_C*100, 
                   y = ..density..)) + 
        geom_histogram(fill = "skyblue",
                       color = "black",
                       alpha = 1,
                       breaks = seq(min(dataset$G_plus_C)*100,
                                    max(dataset$G_plus_C)*100,
                                    length.out = 51)) + 
        labs(x = paste0("G+C content (%)", "\n",
                        "Mean = ", format(round(mean(dataset$G_plus_C)*100, 2), 
                                          nsmall = 2), " ",
                        "St.Dev = ", format(round(sd(dataset$G_plus_C)*100, 2), 
                                            nsmall = 2)),
             y = "Density",
             title = title)

ggsave(width=15, height=8, 
       filename = paste0("../figures/", species, "/GC_hist.", save.as), 
       plot = Plot)
print("GC histogram plot done!", quote = FALSE)

#-----------------------------
# A+T content 
Plot <- filtered.df %>%
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

ggsave(width=15, height=8, 
       filename = paste0("../figures/", species, "/AT_hist.", save.as), 
       plot = Plot)
print("AT histogram plot done!", quote = FALSE)

#-----------------------------
# G-C vs. A-T 
x = filtered.df$A_minus_T
y = filtered.df$G_minus_C

pdf(width=6, height=6, file=paste0("../figures/", species, "/G-C_vs_A-T.pdf"))
smoothScatter(y=y, x=x,
              nrpoints=100, nbin=1000,
              bandwidth=c(diff(range(x))/500, diff(range(y))/500),
              xlab="A-T", ylab="G-C",
              colramp=colfun, main=NULL)
plot.save <- dev.off()
print("G-C vs. A-T plot done!", quote = FALSE)

#---------------------
# GC-ratio vs. AT-ratio
pdf(width=6, height=6, file=paste0("../figures/", species, "/GC-ratio_vs_AT-ratio.pdf"))
smoothScatter(y=filtered.df$GC_ratio,
              x=filtered.df$AT_ratio,
              nrpoints=100, nbin=1000,
              xlab="AT ratio", ylab="GC ratio",
              colramp=colfun, main=NULL)
plot.save <- dev.off()
print("GC vs. AT ratio plot done!", quote = FALSE)

#---------------------
# GC-skew vs. AT-skew
pdf(width=6, height=6, file=paste0("../figures/", species, "/GC-skew_vs_AT-skew.pdf"))
smoothScatter(y=filtered.df$GC_skew,
              x=filtered.df$AT_skew,
              nrpoints=100, nbin=1000,
              xlab="AT skew", ylab="GC skew",
              colramp=colfun, main=NULL)
plot.save <- dev.off()
print("GC vs. AT skew plot done!", quote = FALSE)

#---------------------
# GC-skew
Plot <- filtered.df %>%
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

ggsave(width=15, height=8, 
       filename = paste0("../figures/", species, "/GC_skew.", save.as),
       plot = Plot)
print("GC skew histogram plot done!", quote = FALSE)

# AT-skew
Plot <- filtered.df %>%
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

ggsave(width=15, height=8, 
       filename = paste0("../figures/", species, "/AT_skew.", save.as),
       plot = Plot)
print("AT skew histogram plot done!", quote = FALSE)