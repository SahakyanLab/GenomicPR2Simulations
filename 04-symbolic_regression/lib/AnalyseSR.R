AnalyseSR <- R6::R6Class(
    classname = "AnalyseSR",
    public = list(
        #' @field df_train Data.frame. Training data from the PR-2 compliant 
        #'  results from simulation.
        df_train = NULL,

        #' @field df_test Data.frame. Test data from the PR-2 compliant 
        #'  results from simulation.
        df_test = NULL,

        initialize = function(trek_scale){
            if(!missing(trek_scale)) private$trek_scale <- trek_scale
        },

        #' @description
        #' Analyse results of symbolic regression equations.
        #' @return None.
        run_process = function(){
            private$import_data()
            private$save_rate_const_files()
            private$get_plots(train = TRUE)
            private$get_plots(train = FALSE)
            private$get_mut_rates()
            private$plot_against_real_rate_const()
        }
    ), 
    private = list(
        #' @field trek_scale Boolean. If TRUE, will convert into trek scale.
        trek_scale = TRUE,

        #' @field lynch_rates Data.frame of mutation rate constants of real life species.
        lynch_rates = NULL,

        #' @description
        #' Import training and testing data from simulations.
        #' @return None.
        import_data = function(){
            t1 <- Sys.time()
            cur.msg <- "Loading in data sets"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # training data
            self$df_train$data <- read.csv("../data/Train/simulation_batches_compliant.csv")
            self$df_train$rate_const <- colnames(self$df_train$data)

            # testing data
            self$df_test$data <- read.csv("../data/Test/simulation_batches_compliant.csv")
            self$df_test$data <- self$df_test$data[2:length(self$df_test$data)]
            self$df_test$rate_const <- colnames(self$df_test$data)

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Obtain mutation rates from Table 1 in the following paper
        #'  https://www.pnas.org/content/107/3/961
        #' @return None.
        get_mut_rates = function(){
            private$lynch_rates <- read.csv(
                file = paste0("../../02-simulations/data/Raw/Other_species/", 
                              ifelse(private$trek_scale, "/Trek_scale_", "/"),
                              "Lynch-2010-mutation-rates_FREQUENCY.csv"),
                header = TRUE
            )
        },

        #' @description
        #' Get mutation rate constant files for symbolic regression.
        #' @return None.
        save_rate_const_files = function(){
            t1 <- Sys.time()
            cur.msg <- "Getting rate constants for symbolic regression"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            for(dataset in c("df_train", "df_test")){
                for(i in 1:length(self[[dataset]]$rate_const)){
                    write.csv(
                        x = private$equations(
                            data = self[[dataset]]$data,
                            var = self[[dataset]]$rate_const[i]
                        ),
                        row.names = FALSE,
                        file = ifelse(dataset == "df_train", 
                            paste0("../data/Train/EUREQA_prediction_", 
                                   self[[dataset]]$rate_const[i],".csv"),
                            paste0("../data/Test/EUREQA_prediction_", 
                                   self[[dataset]]$rate_const[i],".csv")
                        )
                    )
                }
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Equations obtained through Eureqa symbolic regression.
        #' @param data Data.frame containing mutation rate constants.
        #' @param var Character vector. Mutation rtae constant as a variable.
        #' @return Data.frame.
        equations = function(data, var){
            if(!is.character(var)) stop("var needs to be a character vector.")
            y_pred <- switch(var,
            'kag' = {
                0.348387923987881 + 0.816734963882394*data["ktc"] + 
                0.651414343680873*data["kta"] + 0.609291658916403*data["kgc"] + 
                0.605797402160463*data["kga"] - 0.618906078160872*data["kat"] - 
                0.67329756572978*data["kcg"] - 0.721891543764904*data["kct"]
            },
            'kat' = {
                0.241665477033022 + 0.897500217963755*data["kta"] + 
                0.514297637755254*data["ktg"] + 0.46775568703087*data["ktc"] + 
                0.463842653474543*data["kca"] + 0.393996535254421*data["kga"] - 
                0.479882925548398*data["kag"] - 0.484672348490665*data["kct"] - 
                0.487621790665817*data["kac"] - 0.516076484396038*data["kgt"]
            },
            'kac' = {
                0.22558936599309 + 0.852687563679462*data["kcg"] + 
                0.801456741790903*data["ktg"] + 0.794496499861356*data["kca"] + 
                0.787171624286753*data["kta"] - 0.753938685112867*data["kat"] - 
                0.837621207499993*data["kgc"] - 0.892761983670905*data["kgt"]
            },
            'kga' = {
                0.105020824396609 + 0.736399072631546*data["kag"] + 
                0.736399072631546*data["kat"] + 0.730854831389399*data["kct"] + 
                0.668840893845205*data["kcg"] - 0.643445226498024*data["kgc"] - 
                0.668840893845205*data["ktc"] - 0.738109376474998*data["kta"]
            },
            'kgt' = {
                0.195602186906957 + 0.795944042083753*data["kcg"] + 
                0.762762200816868*data["kta"] + 0.757989214023514*data["ktg"] + 
                0.741388062485789*data["kca"] - 0.709523759605975*data["kat"] - 
                0.758850686673801*data["kgc"] - 0.801757287841475*data["kac"]
            },
            'kgc' = {
                0.0266522601542158 + 0.826978095447326*data["kcg"] + 
                0.489286002530125*data["kag"] + 0.453348892899995*data["kct"] + 
                0.449467306920196*data["kca"] + 0.414417817271229*data["ktg"] - 
                0.385407414591448*data["ktc"] - 0.418343418385765*data["kac"] - 
                0.430961689763262*data["kgt"] - 0.432476148161079*data["kga"]
            },
            'kta' = {
                0.828368985921817*data["kat"] + 0.530245079564161*data["kct"] + 
                0.509302963224534*data["kgt"] + 0.499135884778579*data["kag"] + 
                0.467813865068305*data["kac"] - 0.0460955436591221 - 
                0.42990907018686*data["kca"] - 0.446701251437334*data["ktc"] - 
                0.451430012359511*data["ktg"] - 0.46984114233575*data["kga"]
            },
            'ktg' = {
                0.112574827222794 + 0.88412124636153*data["kgt"] + 
                0.841207273827614*data["kac"] + 0.807204591082419*data["kgc"] + 
                0.794765122731266*data["kat"] - 0.792020815272827*data["kca"] - 
                0.811438663741675*data["kcg"] - 0.835577248635635*data["kta"]
            },
            'ktc' = {
                0.302524932124489 + 0.81778486498004*data["kag"] + 
                0.765641304372685*data["kct"] + 0.620216273339829*data["kat"] + 
                0.605049864248977*data["kcg"] - 0.641586087709269*data["kga"] - 
                0.739907793867365*data["kta"] - 0.757175993792224*data["kgc"]
            },
            'kca' = {
                0.135000854348154 + 0.8899122260307*data["kgt"] + 
                0.826421028836978*data["kgc"] + 0.794875984525169*data["kac"] +
                0.756327474198048*data["kat"] - 0.786605858750487*data["kta"] -
                0.787397685580583*data["ktg"] - 0.840294389850611*data["kcg"]
            },
            'kcg' = {
                0.393304935908064 + 0.764253923040039*data["kgc"] + 
                0.426508013212467*data["kac"] + 0.422715889530294*data["kgt"] + 
                0.415201669034913*data["ktc"] + 0.393304935908064*data["kga"] - 
                0.393304935908064*data["kct"] - 0.417852963826674*data["ktg"] - 
                0.509603310796301*data["kag"] - 0.509603310796301*data["kca"]
            },
            'kct' = {
                0.23960009195787 + 0.787158743182794*data["kgc"] + 
                0.76912313021918*data["kga"] + 0.727489918806558*data["ktc"] + 
                0.695649662821069*data["kta"] - 0.0369776458407348*data["kca"] - 
                0.643942083249785*data["kat"] - 0.71631868285512*data["kcg"] - 
                0.797168947756005*data["kag"]
            })
            
            # simplified 
            # y_pred <- switch(var,
            #     "kag" = {
            #         -0.6*(data["kcg"]-data["kgc"]+data["kat"]-
            #         data["kta"]+data["kct"]-1.3*data["ktc"]-
            #         data["kga"]-0.6)
            #     },
            #     "kga" = {
            #         0.7*(data["kcg"]-data["kgc"]+data["kat"]-
            #         data["kta"]+data["kct"]-data["ktc"]+
            #         data["kag"]+0.1)
            #     },
            #     "ktc" = {
            #         0.8*(0.7*data["kcg"]-data["kgc"]+0.8*data["kat"]-
            #         data["kta"]+data["kag"]-0.8*data["kga"]+
            #         data["kct"]+0.4)
            #     },
            #     "kct" = {
            #         -0.8*(data["kcg"]-data["kgc"]+0.8*data["kat"]-
            #         data["kta"]+data["kag"]-data["kga"]+0.1*data["kca"]-
            #         data["ktc"]-0.3)
            #     },
            #     "kac" = {
            #         -0.8*(-data["kcg"]+data["kgc"]+data["kat"]-
            #         data["kta"]+data["kgt"]-data["ktg"]-
            #         data["kca"]-0.3)
            #     },
            #     "kca" = {
            #         0.8*(data["kgc"]-data["kcg"]+data["kat"]-
            #         data["kta"]+data["kgt"]-data["ktg"]+
            #         data["kac"]+0.2)
            #     },
            #     "kgt" = {
            #         -0.8*(data["kgc"]-data["kcg"]+data["kat"]-
            #         data["kta"]-data["ktg"]+data["kac"]-
            #         data["kca"]-0.2)
            #     },
            #     "ktg" = {
            #         0.8*(data["kgc"]-data["kcg"]+data["kat"]-
            #         data["kta"]+data["kgt"]+data["kac"]-
            #         data["kca"]+0.1)
            #     },
            #     "kat" = {
            #         -0.5*(-2*data["kta"]+data["kct"]-data["ktc"]+
            #         data["kag"]-data["kga"]+data["kgt"]-data["ktg"]+
            #         data["kac"]-data["kca"]-0.5)
            #     },
            #     "kta" = {
            #         0.5*(+2*data["kat"]+data["kct"]-data["ktc"]+
            #         data["kag"]-data["kga"]+data["kgt"]-data["ktg"]+
            #         data["kac"]-data["kca"]-0.1)
            #     },
            #     "kgc" = {
            #         -0.4*(-2*data["kcg"]+data["ktc"]-data["kct"]+
            #         data["kga"]-data["kag"]+data["kgt"]-data["ktg"]+
            #         data["kac"]-data["kca"]-0.06)
            #     },
            #     "kcg" = {
            #         0.4*(2*data["kgc"]+data["ktc"]-data["kct"]+
            #         data["kga"]-data["kag"]+data["kgt"]-
            #         data["ktg"]+data["kac"]-data["kca"]+0.9)
            #     }
            # )
            
            # true values
            results <- data.frame(
                y_true = as.data.frame(data[, var]),
                y_pred = y_pred
            )
            colnames(results) <- c("y_true", "y_pred")
            return(results)
        },

        #' @description
        #' Get plots.
        #' @return None.
        get_plots = function(train = TRUE){
            t1 <- Sys.time()
            cur.msg <- paste0("Plotting results for ", ifelse(train, 
                              "training", "testing"), " dataset")
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            files <- list.files(
                path = paste0("../data/", ifelse(train, "Train", "Test")), 
                pattern = "EUREQA_prediction_*",
                full.names = TRUE
            )
            file.name <- stringr::str_extract(
                string = files, 
                pattern = "(?<=EUREQA_prediction_).*(?=.csv)"
            )

            # import files
            data.sets <- lapply(1:length(files), function(x){
                result <- read.csv(files[x], header = TRUE)
                result <- result %>% 
                    dplyr::mutate(
                        R2 = paste0("R = ", private$calc_r_squared(
                            y_true, y_pred
                        )),
                        rates = file.name[x]
                    )
                return(result)
            })
            results <- do.call(rbind, data.sets)

            results <- data.frame(
                results,
                xpos = -Inf,
                ypos =  Inf,
                hjustvar = 0,
                vjustvar = 1
            )

            # main plots
            p1 <- results %>% 
                dplyr::group_by(rates) %>% 
                ggplot(aes(x = y_true, y = y_pred)) +
                geom_point(
                    size = 1, 
                    alpha = 0.6) + 
                facet_wrap(
                    ~rates, 
                    ncol = 3) + 
                geom_smooth(
                    method = "lm", 
                    formula = y ~ x) +
                geom_abline(
                    slope = 1, 
                    linetype = "dashed") + 
                geom_text(aes(x = xpos, 
                        y = ypos,
                        hjust = hjustvar,
                        vjust = vjustvar,
                        label = R2
                )) +
                coord_cartesian(
                    xlim = c(-1, 2.5),
                    ylim = c(-1, 2.5)
                ) +
                theme_bw() + 
                labs(
                    x = "Y_True",
                    y = "Y_Pred") + 
                theme(
                    strip.text.x = element_text(size = 12),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )

            dir.create(
                path = paste0("../figures/", ifelse(train, 
                    "Train/", "Test/"
                )),
                showWarnings = FALSE,
                recursive = TRUE
            )
            # save_full_name = paste0(
            #     "../figures/", ifelse(train, 
            #     "Train/SymbolicRegressionCorrelation_TRAINING_SIMPLIFIED.", 
            #     "Test/SymbolicRegressionCorrelation_TEST_EUREQA_SIMPLIFIED."),
            #     "pdf"
            # )
            save_full_name = paste0(
                "../figures/", ifelse(train, 
                "Train/SymbolicRegressionCorrelation_TRAINING.", 
                "Test/SymbolicRegressionCorrelation_TEST_EUREQA."),
                "pdf"
            )
            ggsave(
                filename = save_full_name,
                plot = p1,
                height = 9, width = 7
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Calculates r-squared between two variables.
        #' @x x Numeric vector of mutation rate constant values.
        #' @y y Numeric vector of mutation rate constant values.
        #' @return Numeric vector of R-squared value.
        calc_r_squared = function(x, y){
            return(signif(summary(lm(x~y))$r.squared, digits=3))
        },

        #' @description
        #' Plot the predicted mutation rate constants vs. true values
        #' from the 17 species.
        #' @return None.
        plot_against_real_rate_const = function(){
            t1 <- Sys.time()
            cur.msg <- "Plotting true vs. predicted mutation rate constants"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            private$lynch_rates <- as_tibble(private$lynch_rates) %>%
                dplyr::mutate(
                    MUT = plyr::mapvalues(
                        x = private$lynch_rates$MUT,
                        from = c("AG|TC", "CT|GA", "AT|TA", 
                                 "CA|GT", "AC|TG", "CG|GC"),
                        to = c("kag|ktc", "kct|kga", "kat|kta", 
                               "kca|kgt", "kac|ktg", "kcg|kgc")
                    )
                ) %>% 
                tidyr::separate(MUT, c("MUT_1", "MUT_2"))

            # generate plots
            rates <- c("kag", "ktc", "kct", "kga", 
                       "kat", "kta", "kca", "kgt", 
                       "kac", "ktg", "kcg", "kgc")
            lynch.rates.columns <-  colnames(
                private$lynch_rates[3:length(private$lynch_rates)]
            )

            df <- private$lynch_rates %>% 
                tidyr::pivot_longer(1:2) %>% 
                dplyr::select(-name) %>% 
                dplyr::rename(rates = value) %>% 
                tidyr::pivot_longer(-rates) %>% 
                tidyr::pivot_wider(names_from = rates, values_from = value)

            results <- lapply(1:length(lynch.rates.columns), function(cols){
                out <- sapply(rates, function(i){
                    return(private$equations(data = df[cols, ], var = i))
                })

                filtered.df <- as.data.frame(t(out)) %>% 
                    dplyr::mutate(dplyr::across(where(is.list), unlist)) %>% 
                    dplyr::rename_with(~c("y_true", "y_pred")) %>% 
                    dplyr::mutate(
                        species = lynch.rates.columns[cols],
                        rates = colnames(out)
                    )

                return(filtered.df)
            })
            results <- do.call(rbind, results)
            rownames(results) <- NULL

            euk <- c(
                "H.sapiens", "D.melanogaster", "C.elegans", "A.thaliana",	
                "S.cerevisiae", "M.m.domesticus", "D.pulex", "P.pacificus",	
                "D.magna", "P.troglodytes", "A.nancymaae"
            )
            prok <- c(
                "E.coli", "P.luminescens ATCC29999", "T.turnerae", 
                "M.smegmatis", "P.fluorescens ATCC948", "R. toruloides"
            )

            results <- results %>% 
                dplyr::mutate(
                    kingdom = ifelse(species %in% euk, 
                        "eukaryotes", 
                        "prokaryotes"
                    ),
                    kingdom.col = ifelse(kingdom == "eukaryotes", 
                        "purple", 
                        "darkgreen"
                    )
                )
            
            # one plot with 12 rate constants for each species
            p <- results %>% 
                dplyr::group_by(species) %>% 
                dplyr::mutate(R2 = paste0("R = ", 
                    private$calc_r_squared(y_true, y_pred))) %>% 
                ggplot(aes(x = y_true, y = y_pred, col = rates)) +
                geom_abline(
                    intercept = 0, 
                    slope = 1, 
                    linetype = "dashed") + 
                geom_point() + 
                facet_wrap(
                    ~species, 
                    ncol = 6) + 
                coord_cartesian(
                    xlim = c(0, 1.7), 
                    ylim = c(0, 1.7)) +
                labs(
                    x = "Y_True", 
                    y = "Y_Pred") + 
                theme_bw() + 
                theme(
                    strip.text.x = element_text(size = 9),
                    axis.text.x = element_text(angle = 90, 
                                               vjust = 0.5, 
                                               hjust = 1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )

            dir.create(
                path = "../figures/",
                showWarnings = FALSE,
                recursive = TRUE
            )
            ggsave(
                filename = paste0("../figures", 
                    ifelse(private$trek_scale, "/Trek_scale_", "/"), 
                    "ExperimentalRateConstants_by_species.pdf"), 
                plot = p,
                height = 6, width = 11
            )

            # one plot with all species for each rate constant
            p <- results %>% 
                dplyr::group_by(rates) %>% 
                dplyr::mutate(R2 = paste0("R = ", 
                    private$calc_r_squared(y_true, y_pred))) %>% 
                ggplot(aes(x = y_true, y = y_pred, col = species)) +
                geom_abline(
                    intercept = 0, 
                    slope = 1, 
                    linetype = "dashed") + 
                geom_point() + 
                facet_wrap(
                    ~rates, 
                    ncol = 6) + 
                # coord_cartesian(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
                coord_cartesian(
                    xlim = c(0, 1.5), 
                    ylim = c(0, 1.5)) +
                labs(
                    x = "Y_True", 
                    y = "Y_Pred") + 
                theme_bw() + 
                theme(
                    strip.text.x = element_text(size = 12),
                    axis.text.x = element_text(angle = 90, 
                                               vjust = 0.5, 
                                               hjust = 1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )
            
            ggsave(
                filename = paste0("../figures", 
                    ifelse(private$trek_scale, "/Trek_scale_", "/"), 
                    "ExperimentalRateConstants_by_rates.pdf"), 
                plot = p,
                height = 5, width = 13
            )

            # col by kingdoms
            p <- results %>% 
                dplyr::group_by(rates) %>% 
                dplyr::mutate(R2 = paste0("R = ", 
                    private$calc_r_squared(y_true, y_pred))) %>% 
                ggplot(aes(x = y_true, y = y_pred)) +
                geom_abline(
                    intercept = 0, 
                    slope = 1, 
                    linetype = "dashed") + 
                geom_point(
                    colour = results$kingdom.col) + 
                facet_wrap(
                    ~rates, 
                    ncol = 6) + 
                geom_text(
                    aes(x = 1.15, y = 0.05, label = R2),
                    show.legend = FALSE) + 
                coord_cartesian(
                    xlim = c(0, 1.5), 
                    ylim = c(0, 1.5)) +
                labs(
                    x = "Y_True", 
                    y = "Y_Pred") + 
                theme_bw() + 
                theme(
                    strip.text.x = element_text(size = 12),
                    axis.text.x = element_text(angle = 90, 
                                               vjust = 0.5, 
                                               hjust = 1),
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                )
            
            ggsave(
                filename = paste0("../figures", 
                    ifelse(private$trek_scale, "/Trek_scale_", "/"), 
                    "ExperimentalRateConstants_by_rates_bykingdom.pdf"), 
                plot = p,
                height = 4, width = 10
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        }
    )
)