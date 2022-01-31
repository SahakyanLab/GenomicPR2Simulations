args <- commandArgs(trailingOnly = TRUE)
my_path <- as.character(args[1])
save_file <- as.character(args[2])
TEST <- as.logical(as.character(args[3]))
setwd(my_path)

suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))

# ---------------------
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
# ---------------------

if(TEST){
  files <- list.files(path = "../data/Test/", 
                      pattern = "EUREQA_prediction_*") 
} else {
  files <- list.files(path = "../data/Training/", 
                      pattern = "EUREQA_prediction_*")
}

# extract file name 
file.name <- str_extract(string = files, pattern = "(?<=EUREQA_prediction_).*(?=.csv)") 

# import files
data.sets <- lapply(1:length(files), function(x){
  result <- read.table(file = ifelse(TEST,
                              paste0("../data/Test/", files[x]),
                              paste0("../data/Training/", files[x])), 
             sep = ",", header = TRUE)

  result <- result %>% 
    mutate(
      R2 = paste0("R = ", r.squared(y_true, y_pred)),
      rates = file.name[x]
    ) 
})

results <- do.call(rbind, data.sets)

# main plots
p1 <- results %>% 
  group_by(rates) %>% 
  ggplot(aes(x = y_true, y = y_pred)) +
  geom_point(size = 1, alpha = 0.6) + 
  facet_wrap(~rates, ncol = 3) + 
  geom_smooth(method = "loess", formula = y ~ x) +
  geom_abline(slope = 1, linetype = "dashed") + 
  coord_cartesian(
    xlim = c(-1, 2.5),
    ylim = c(-1, 2.5)
  ) +
  theme_bw() + 
  labs(x = "Y_True",
       y = "Y_Pred") + 
  theme(
    strip.text.x = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  filename = paste0(
    "../figures/", ifelse(
      TEST, "Test/SymbolicRegressionCorrelation_TEST_EUREQA.",
      "Train/SymbolicRegressionCorrelation_TRAINING."
    ), save_file),
  plot = p1,
  height = 9, width = 7
)