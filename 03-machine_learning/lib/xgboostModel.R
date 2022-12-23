xgboostModel <- R6::R6Class(
    classname = "xgboostModel",
    public = list(
        #' @field sim_run Data.frame object of the simulation output.
        sim_run = NULL,

        #' @field sim_run_Rates Data.frame of only the Label and rate constant columns.
        sim_run_rates = NULL,

        #' @field df_compliant Tibble of the simulations' PR-2 compliant cases.
        df_compliant = NULL,

        #' @field df_shuffled Tibble. Randomly shuffled rows of simulation results.
        df_shuffled = NULL,

        #' @field to_sample Numeric vector. Indices to sample the full simulation from.
        to_sample = NULL,

        initialize = function(species, cv, cvrep){
            if(!missing(species)) private$species <- species
            if(!missing(cv)) private$cv <- cv
            if(!missing(cvrep)) private$cvrep <- cvrep
            private$get_pr2_compliance_values()
        },

        #' @description
        #' Grid-search of xgboost hyperparameters to find the best model.
        #' @param seed Numeric vector to set the seed for the random number 
        #'  generator for extracting the non-compliant PR-2 cases.
        #' @param ml_seed Numeric vector to set the seed for the random number
        #'  generator for the machine learning process.
        #' @return None.
        find_best_model = function(seed = 2022, ml_seed = 1234){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Importing required data for xgboost training"
            l <- paste0(rep(".", 60-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            private$import_sim_run(symbolic_regression = FALSE)
            
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")

            private$filter_compliant_data()
            private$filter_noncompliant_data(seed = seed)

            # progress bar
            start.time <- Sys.time()
            cur.msg <- "Fitting xgboost model to data"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            # time taken for full processing for this experiment
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cur.msg <- paste0("\n",
                "Cross-validations:  ", private$cv, "\n",
                "Repeated CVs:       ", private$cvrep, "\n",
                "Seed set:           ", ml_seed
            )
            cat("Quick infosheet:", "\n", cur.msg, "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))

            private$train_model(
                seed = ml_seed, 
                hypertuning = TRUE
            )

            # time taken for full processing for this experiment
            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
        },

        #' @description
        #' Get average ML model behaviour
        #' @param seed Numeric vector to set the seed for the random number 
        #'  generator for extracting the non-compliant PR-2 cases.
        #' @param avg_ml_seed Numeric vector to set the seed for the random number
        #'  generator for the machine learning process.
        #' @return None.
        get_avg_behaviour = function(seed = 2021, avg_ml_seed = 123){
            # progress message
            t1 <- Sys.time()
            cur.msg <- "Importing required data for xgboost training"
            l <- paste0(rep(".", 60-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            private$import_sim_run(symbolic_regression = TRUE)
            private$filter_compliant_data()
            dir.create(
                path = "../../04-symbolic_regression/data/Test/",
                showWarnings = FALSE,
                recursive = TRUE
            )
            write.csv(
                x = self$df_compliant %>% dplyr::select(-Label),
                file = paste0("../../04-symbolic_regression/data/", 
                              "Test/simulation_batches_compliant.csv"),
                row.names = FALSE
            )
            
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")

            start.time <- Sys.time()
            cur.msg <- "Fitting best xgboost model to data for avg behaviour"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(cur.msg, l, "\n", sep = "")

            private$import_sim_run(symbolic_regression = FALSE)
            private$filter_compliant_data()         
            dir.create(
                path = "../../04-symbolic_regression/data/Train/",
                showWarnings = FALSE,
                recursive = TRUE
            )
            write.csv(
                x = self$df_compliant %>% dplyr::select(-Label),
                file = paste0("../../04-symbolic_regression/data/", 
                              "Train/simulation_batches_compliant.csv"),
                row.names = FALSE
            )

            set.seed(seed = seed)
            rng <- sample(x = 1:1000, size = 1000, replace = FALSE)
            df.samples <- lapply(1:length(rng), function(x){
                # progress message
                t1 <- Sys.time()
                cur.msg <- paste0("Running model nr: ", x, "/", length(rng))
                l <- paste0(rep(".", 60-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                # random sampling of non compliant cases
                private$filter_noncompliant_data(seed = rng[x])
                
                # train model with optimised hyper-parameters
                private$train_model(
                    seed = avg_ml_seed, 
                    hypertuning = FALSE
                )

                # return feature importance as data frame
                varimp.dat <- varImp(private$xgb_model)$importance %>%
                    tibble::rownames_to_column() %>%
                    dplyr::rename(
                        rates = rowname,
                        importance = Overall) %>%
                    dplyr::arrange(rates) %>%
                    dplyr::select(importance)

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
                return(varimp.dat)
            })
            featImp.all <- do.call(cbind, df.samples)

            saveRDS(
                featImp.all,
                file = "../../data/XGBoost/featimp_1k-iters_raw.Rdata"
            )

            # rename columns
            colnames(featImp.all) <- paste0(
                "output_", seq(from = 1, to = dim(featImp.all)[2], by = 1)
            )

            # sort rate constants alphabetically
            rates.sorted <- sort(colnames(df[[1]])[2:ncol(df[[1]])])

            # transpose tibble
            featImp.all <- featImp.all %>%
                as_tibble() %>%
                tibble::rownames_to_column() %>%
                tidyr::pivot_longer(-rowname) %>%
                tidyr::pivot_wider(
                    names_from = rowname,
                    values_from = value) %>%
                dplyr::select(-name) %>%
                dplyr::rename_with(~rates.sorted)

            saveRDS(
                featImp.all,
                file = "../../data/XGBoost/featimp_1k-iters_sorted.Rdata"
            )

            # statistical calculations
            data <- rates.sorted %>%
                as_tibble() %>%
                dplyr::mutate(
                    mean = unname(colMeans(featImp.all, na.rm = TRUE)),
                    sd   = unname(apply(featImp.all, 2, sd, na.rm = TRUE)),
                    max  = mean+sd,
                    min  = mean-sd) %>%
                dplyr::rename(rate = value)

            saveRDS(
                featImp.all,
                file = "../data/XGBoost/featimp_1k-iters_summary.Rdata"
            )

            # plot average feature importance
            data.plot <- data %>%
                dplyr::arrange(mean) %>%
                dplyr::mutate(rate = forcats::fct_inorder(rate)) %>%
                ggplot() + 
                geom_bar(
                    aes(x = rate, y = mean),
                    stat = "identity", 
                    fill = "skyblue", 
                    alpha = 0.7) + 
                geom_errorbar(
                    aes(x = rate, ymin = min, ymax = max),
                    width = 0.4, 
                    colour = "orange", 
                    alpha = 1, 
                    size = 1.5) +
                coord_flip()
                labs(
                    x = "Importance",
                    y = "Features",
                    title = "Feature Importance averaged over 1000 iterations"
                )

            ggsave(
                filename = "../figures/XGBoost/xgbtree_featimport_average_plots.pdf",
                plot = data.plot
            )

            ordering <- names(sort(apply(featImp.all, 2, median, na.rm = TRUE)))
            p1 <- featImp.all %>% 
                tidyr::gather(key, value) %>% 
                dplyr::mutate(key = factor(key, levels = ordering)) %>% 
                ggplot(aes(x = value, y = key, fill = key)) + 
                geom_boxplot() +
                geom_jitter(
                    color = "black", 
                    size = 0.4, 
                    alpha = 0.9) +
                labs(
                    x = "Importance",
                    y = "Features",
                    title = "Feature Importance averaged over 1000 iterations"
                )
            ggsave(
                filename = "../figures/XGBoost/xgbtree_featimport_boxplots.pdf",
                plot = p1
            )

            final.t <- Sys.time() - start.time
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))
            cat("Final time taken:", signif(final.t[[1]], digits = 3), 
                attr(final.t, "units"), "\n")
            cat(paste(c(rep("-", 70), "\n"), collapse = ""))           
        }
    ),
    private = list(
        #' @field species Character vector of c("eukaryotes", "prokaryotes", "viruses").
        species = "eukaryotes",

        #' @field cv Numeric vector. Number of cross-validation processes.
        cv = 6,

        #' @field cvrep Numeric vector. Number of times to repeat the cv process.
        cvrep = 1,

        #' @field species_tolerance Data.frame of the PR-2 compliance values for each kingdom.
        species_tolerance = NULL,

        #' @field best_model xgboost model with the highest AUCROC.
        best_model = NULL,

        #' @field xgb_model xgboost model with the highest AUCROC plus additional information.
        xgb_model = NULL,

        #' @description
        #' Get PR-2 compliance values for each kingdom.
        #' @return None.
        get_pr2_compliance_values = function(){
            private$species_tolerance$prokaryotes <- read.csv(
                file = paste0("../../01-genome_composition/data/", 
                              "01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv"),
                header = TRUE
            )
            private$species_tolerance$eukaryotes <- read.csv(
                file = paste0("../../01-genome_composition/data/", 
                              "02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv"),
                header = TRUE
            )
            private$species_tolerance$viruses <- read.csv(
                file = paste0("../../01-genome_composition/data/", 
                              "03-Viruses/PR2_compliance/PR2_fluctuations.csv"),
                header = TRUE
            )
        },

        #' @description
        #' Import simulation output.
        #' @return None.
        import_sim_run = function(symbolic_regression = FALSE){
            if(symbolic_regression){
                self$sim_run <- readRDS(
                    file = paste0("../../04-symbolic_regression/data/", 
                                  "Simulation/simulation_batches.Rdata")
                )
            } else {
                self$sim_run <- readRDS(
                    file = paste0("../../02-simulations/data/Main_Simulation/", 
                                  "Non_Symmetric-uniform/Non_Symmetric-uniform-", 
                                  "scaling-1.Rdata")
                )
            }
        },

        #' @description
        #' Filters the simulation dataset for PR-2 compliant cases for use in ML model.
        #' @return None.
        filter_compliant_data = function(){
            t1 <- Sys.time()
            cur.msg <- "Filtering simulation data for PR-2 compliant cases"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # get chargaff compliancy region
            fluc.species <- private$species_tolerance[[private$species]]

            # categorise data as compliant and non-compliant cases
            y <- as_tibble(self$sim_run) %>%
                dplyr::mutate(Label = ifelse((
                    (ATskew >= fluc.species[fluc.species$metadata == "AT_skew","mean"]-
                        fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                    (ATskew <= fluc.species[fluc.species$metadata == "AT_skew","mean"]+
                        fluc.species[fluc.species$metadata == "AT_skew","st.dev"]) & 
                    (GCskew >= fluc.species[fluc.species$metadata == "GC_skew","mean"]-
                        fluc.species[fluc.species$metadata == "GC_skew","st.dev"]) & 
                    (GCskew <= fluc.species[fluc.species$metadata == "GC_skew","mean"]+
                        fluc.species[fluc.species$metadata == "GC_skew","st.dev"])), 
                    1, 0)) %>%
                dplyr::relocate(Label, .before = kag) %>%
                dplyr::mutate(
                    Label = as.factor(as.numeric(Label)),
                    Label = forcats::fct_recode(Label, "NO" = "0", "YES" = "1"))

            # select only columns of Label and all rate constants
            self$sim_run_rates <- y %>%
                dplyr::select(
                    contains("Label"),
                    starts_with("k")) %>% 
                tidyr::drop_na()
            
            # filter out all compliant cases
            self$df_compliant <- self$sim_run_rates %>%
                dplyr::filter(Label == "YES")

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n") 
        },

        #' @description
        #' Filters the simulation dataset for PR-2 non-compliant cases for use in ML model.
        #' @param seed Numeric vector. Random number generator for reproducibility.
        #' @return None.
        filter_noncompliant_data = function(seed){
            # random sampling
            df.non.compliant <- self$sim_run_rates %>% 
                dplyr::filter(Label == "NO")

            set.seed(seed = seed)
            self$to_sample <- sample(
                x = nrow(df.non.compliant),
                size = nrow(self$df_compliant),
                replace = FALSE)

            df.non.compliant <- self$sim_run_rates %>%
                dplyr::filter(Label == "NO") %>%
                dplyr::slice(self$to_sample)

            # combine into df
            self$sim_run_rates <- rbind(self$df_compliant, df.non.compliant)

            # Randomly shuffle rows of tbl
            self$df_shuffled <- self$sim_run_rates[sample(nrow(self$sim_run_rates)),]
        },

        #' @description
        #' Train xgboost model.
        #' @param seed Numeric vector. Random number generator for reproducibility.
        #' @param hyerptuning Boolean. If TRUE, performs grid-search 
        #'  approach for hyperparameter tunings.
        #' @return None.
        train_model = function(seed = 123, hypertuning = TRUE){
            # check if classes are correct
            if(cv%%1!=0 | cv<=0) 
                stop("Cross-validation must be a positive integer.")
            if(cvrep%%1!=0 | cvrep<=0) 
                stop("cvrep must be a positive integer.")
            if(seed%%1!=0 | seed<=0) 
                stop("Seed must be a positive integer.")

            df.norm <- apply(self$df_shuffled[2:ncol(self$df_shuffled)], 
                             2, scale, center = TRUE, scale = TRUE)
            df.norm <- as_tibble(df.norm) %>% 
                dplyr::mutate(Label = self$df_shuffled$Label, .before = 1)
            
            set.seed(seed)
            if(hypertuning){
                # hyper-parameter search
                xgb.grid <- expand.grid(
                    nrounds = c(200, 500, 1000, 2000, 3000, 
                                4000, 5000, 7000, 9000, 11000, 
                                13000, 15000, 17000, 20000), # boosting iterations
                    max_depth = c(5, 6, 8, 10), # max tree depth
                    colsample_bytree = 1, # subsample ratio of columns
                    eta = c(0.005, 0.01, 0.02, 0.1), # shrinkage
                    gamma = c(0, 0.1, 0.2), # min loss reduction
                    min_child_weight = c(1, 5, 10), # min sum of instance weight
                    subsample = c(0.3, 0.4, 0.6, 0.8) # subsample percentage
                )
            } else {
                # hyper-parameter search
                xgb.grid <- expand.grid(
                    nrounds = private$best_model$nrounds, # boosting iterations
                    max_depth = private$best_model$max_depth, # max tree depth
                    colsample_bytree = private$best_model$colsample_bytree, # subsample ratio of columns
                    eta = private$best_model$eta, # shrinkage
                    gamma = private$best_model$gamma, # min loss reduction
                    min_child_weight = private$best_model$min_child_weight, # min sum of instance weight
                    subsample = private$best_model$subsample # subsample percentage
                )
            }
            
            # define seeds for each CV process
            seeds <- lapply(1:((private$cv*private$cvrep)+1), function(x){
                sample.int(n = 50000, dim(xgb.grid)[1])
            })
            
            # create train control
            xgb.trcontrol <- trainControl(
                method = "repeatedcv",
                repeats = private$cvrep,
                number = private$cv,
                seeds = seeds,
                classProbs = TRUE,
                summaryFunction = twoClassSummary,
                returnData = TRUE,
                verboseIter = ifelse(hypertuning, TRUE, FALSE),
                allowParallel = TRUE,
                savePredictions = TRUE
            )
            
            # train model
            xgb.model <- train(
                Label ~ ., 
                data = df.norm,
                method = "xgbTree",
                trControl = xgb.trcontrol,
                tuneGrid = xgb.grid,
                metric = "ROC",
                verbose = TRUE
            )

            # best model
            private$best_model <- xgb.model$results %>% 
                dplyr::filter(ROC == max(ROC))

            # for average ML model performance
            if(!hypertuning) private$xgb_model <- xgb.model
            
            if(hypertuning){
                # save AUC ROC curve
                res <- MLeval::evalm(xgb.model)
                res.plot <- res$roc
                dir.create(
                    path = "../figures/XGBoost/",
                    showWarnings = FALSE,
                    recursive = TRUE
                )
                ggsave(
                    filename = "../figures/XGBoost/xgbtree_aucroc_plots_new.pdf",
                    plot = res.plot
                )
                
                # feature importance
                feat.imp.plot <- varImp(xgb.model)$importance %>%
                    tibble::rownames_to_column() %>%
                    dplyr::rename("variable" = rowname) %>%
                    dplyr::arrange(Overall) %>%
                    dplyr::filter(dplyr::row_number() > max(dplyr::row_number())-10) %>%
                    dplyr::mutate(variable = forcats::fct_inorder(variable)) %>%
                    ggplot() +
                    geom_col(aes(x = variable, y = Overall)) +
                    coord_flip() +
                    scale_fill_grey() +
                    labs(
                        x = "Importance",
                        y = "Features"
                    )
                ggsave(
                    filename = "../figures/XGBoost/xgbtree_featimport_plots.pdf",
                    plot = feat.imp.plot
                )

                # save best model
                dir.create(
                    path = "../data/XGBoost/",
                    showWarnings = FALSE,
                    recursive = TRUE
                )
                saveRDS(
                    xgb.model, 
                    file ="../data/XGBoost/xgbtree.model"
                )

                # save model results plot
                trellis.device(
                    width = 14,
                    height = 10,
                    device = "pdf",
                    file = "../figures/XGBoost/xgbtree_roc-hyperparams_plots.pdf"
                )
                plot(xgb.model)
                plot.save <- dev.off()
            }
        }
    )
)