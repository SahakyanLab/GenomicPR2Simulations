# load dependencies
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(gridExtra))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(reshape2))

args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save_file <- as.character(args[2])
minimal <- as.character(args[3])
minimal <- as.logical(minimal)
other.species <- as.character(args[4])
setwd(my_path)

if(other.species){
  # import data
  files <- list.files(path = "../../data/Chargaff_Equilibrium/",
                      pattern = "^ChargaffEquilibriumDistribution_Other.*\\.Rdata")
  files <- c("ChargaffEquilibriumDistribution_scaling_zero.Rdata", files)

  data.sets <- lapply(files, function(x){
    readRDS(file = paste0("../../data/Chargaff_Equilibrium/", x))
  })

  # extract file names
  file.name <- lapply(files, function(x){
    str_extract(paste(str_split(string = x, 
                                pattern = "_")[[1]][3:4], 
                      collapse = "_"), pattern = "^[^\\.]+")
  })
  file.name[[1]] <- "H_sapiens"
} else {
  # import data
  files <- list.files(path = "../../data/Chargaff_Equilibrium/",
                      pattern = "^ChargaffEquilibriumDistribution_scaling.*\\.Rdata")
  files <- files[c(5,2,4,1,3)]

  data.sets <- lapply(files, function(x){
    readRDS(file = paste0("../../data/Chargaff_Equilibrium/", x))
  })

  # extract file names
  file.name <- lapply(files, function(x){
    paste(head(str_split(str_extract(string = x, pattern = "([^_]+$)"), 
                        pattern = "")[[1]], n = -6), collapse = "")
  })
}

# obtain plots for each data set
plots <- lapply(1:length(file.name), function(i){
  if(minimal){
    plot.both <- data.sets[[i]] %>%
      as_tibble() %>%
      dplyr::select(-Difference) %>%
      reshape2::melt(.) %>%
      ggplot(aes(x = value, fill = variable)) + 
      geom_histogram(color = "#e9ecef", alpha = 0.5, position = 'identity', 
                     bins = 40, show.legend = FALSE) +
      scale_fill_manual(values = c("#69b3a2", "#404080")) + 
      coord_cartesian(xlim = c(0, 10)) + 
      labs(x = "",
           y = "",
           title = "")
    
    plot.diff <- data.sets[[i]] %>%
      as_tibble() %>%
      dplyr::select(Difference) %>%
      ggplot(aes(x = Difference)) +
      geom_histogram(color = "#e9ecef", alpha = 1, position = 'identity', 
                     bins = 70, show.legend = FALSE) + 
      scale_fill_manual(values = c("#69b3a2")) + 
      coord_cartesian(xlim = c(-3, 3)) + 
      labs(x = "",
           y = "",
           title = "")
  } else {
    plot.both <- data.sets[[i]] %>%
      as_tibble() %>%
      dplyr::select(-Difference) %>%
      reshape2::melt(.) %>%
      ggplot(aes(x = value, fill = variable)) + 
      geom_histogram(color = "#e9ecef", alpha = 0.5, position = 'identity', 
                     bins = 40, show.legend = FALSE) +
      scale_fill_manual(values = c("#69b3a2", "#404080")) + 
      coord_cartesian(xlim = c(0, 10)) + 
      labs(x = "Years (Billion)",
           y = "Freq",
           title = paste0("Scaling: ", file.name[[i]]))
    
    plot.diff <- data.sets[[i]] %>%
      as_tibble() %>%
      dplyr::select(Difference) %>%
      ggplot(aes(x = Difference)) +
      geom_histogram(color = "#e9ecef", alpha = 1, position = 'identity', 
                     bins = 70, show.legend = FALSE) + 
      scale_fill_manual(values = c("#69b3a2")) + 
      coord_cartesian(xlim = c(-3, 3)) + 
      labs(x = "Years (Billion)",
           y = "",
           title = paste0("Scaling: ", file.name[[i]])) 
  }
  return(grid.arrange(plot.both, plot.diff, ncol = 2))
})

if(other.species){
  # save grid of plots
  do.call(grid.arrange, c(plots, ncol = 1)) %>%
    ggsave(width = 16, height = 14,
          file = ifelse(
            minimal,
            paste0("../../figures/Chargaff_Equilibrium/minimal-ChargaffEquilibriumDistribution_OtherSpecies.", save.as),
            paste0("../../figures/Chargaff_Equilibrium/ChargaffEquilibriumDistribution_OtherSpecies.", save.as)))
} else {
  # save grid of plots
  do.call(grid.arrange, c(plots, ncol = 1)) %>%
    ggsave(width = 16, height = 14,
          file = ifelse(
            minimal,
            paste0("../../figures/Chargaff_Equilibrium/minimal-ChargaffEquilibriumDistribution.", save.as),
            paste0("../../figures/Chargaff_Equilibrium/ChargaffEquilibriumDistribution.", save.as)))
}