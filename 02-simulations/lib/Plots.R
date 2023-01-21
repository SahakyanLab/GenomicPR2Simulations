Plots <- R6::R6Class(
    classname = "Plots",
    public = list(
        #' @field sim_run Data.frame object of the simulation output.
        sim_run = NULL,

        initialize = function(muttype, distribution, scale_fac, equilibration, 
                              random_init_bases, edge_skew_cases){
            if(!missing(muttype)) private$muttype <- muttype
            if(!missing(distribution)) private$distribution <- distribution
            if(!missing(scale_fac)) private$scale_fac <- scale_fac
            if(!missing(equilibration)) private$equilibration <- equilibration
            if(!missing(random_init_bases)) private$random_init_bases <- random_init_bases
            if(!missing(edge_skew_cases)) private$edge_skew_cases <- edge_skew_cases
            if(private$edge_skew_cases) private$random_init_bases <- FALSE
        },

        #' @description
        #' Generates all plots for the main simulation.
        #' @param save_as Character vector of c("png", "pdf") to save the plots.
        #' @return None.
        generate_sim_plots = function(save_as = "png"){
            private$get_sim_output() 
            private$save_as <- save_as
            private$plot_GC_hist()
            private$plot_AT_hist()
            private$plot_difference()
            if(private$muttype != "Strand_Symmetric") private$plot_ratios()
            private$plot_skews()
            private$plot_gc_skew_vs_rates()
            private$plot_gc_content_vs_rates(tolerance = FALSE)
            if(private$muttype == "Strand_Symmetric"){
                private$plot_ratio_constants_for_strand_symmetric_cases()
            } else {
                private$generate_equil_const_plots()
                private$plot_ratio_constants()
            }
            if(private$muttype == "Strand_Symmetric"){
                private$plot_skew_evolution()
            }
        },

        #' @description
        #' Generates all plots for the equilibration runs
        #' @param save_as Character vector of c("png", "pdf") to save the plots.
        #' @param other_species Boolean. If TRUE, using mutation rate constants
        #'  from specific types of species.
        #' @return None.
        generate_equil_plots = function(save_as = "png", other_species = FALSE){
            private$get_sim_output(other_species = other_species) 
            private$save_as <- save_as
            private$plot_equilibrations(other_species = other_species)
        }
    ),
    private = list(
        #' @field muttype Character vector of c("Strand_Symmetric", "Non_Symmetric") for
        #'  with and without symmetry-constrained mutation rate constants, respectively. 
        muttype = "Non_Symmetric",

        #' @field distribution Character vector of c("uniform", "normal"). 
        #'  The type of distribution from which the mutation rate constants to randomly draw from. 
        distribution = "uniform",

        #' @field Numeric vector of c(1,2,5,10), indicating the standard deviation of the
        #'  random drawing of mutation rate constants.
        scale_fac = 1,

        #' @field species_tolerance Data.frame of the PR-2 compliance values for each kingdom.
        species_tolerance = NULL,

        #' @field equilibration Boolean. If TRUE, will import simulation results 
        #'  from Chargaff equilibrium runs.
        equilibration = FALSE,
        
        #' @field EQtolerance Numeric vector of the tolerance value for reaching genome equilibration
        EQtolerance = NULL,

        #' @field Character of c("png", "pdf") to save the plots.
        save_as = "png",

        #' @field random_init_bases Boolean. If TRUE, will use random initial nucleotide contents 
        #'  for the generation of a given system.
        random_init_bases = FALSE,

        #' @field edge_skew_cases Boolean. If TRUE, will run a demonstration simulation where the
        #'  initial base contents for each generation is at edge values of (1,1), (-1,1), (1,-1), (-1,-1).
        edge_skew_cases = FALSE,

        #' @description
        #' Load simulation results.
        #' @return None.
        get_sim_output = function(other_species = FALSE){
            t1 <- Sys.time()
            cur.msg <- "Loading in data sets"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))
            
            if(private$equilibration){
                if(other_species){
                    files <- list.files(
                        path = "../data/Chargaff_Equilibrium",
                        pattern = "^ChargaffEquilibriumDistribution_Other.*\\.Rdata",
                        full.names = TRUE
                    )
                    file.name <- stringr::str_remove(
                        string = files,
                        pattern = paste0("../data/Chargaff_Equilibrium/", 
                                         "ChargaffEquilibriumDistribution_Other_")
                    )
                    file.name <- stringr::str_remove(
                        string = file.name, 
                        pattern = ".Rdata"
                    )
                    file.name <- c("H_sapiens", file.name)
                    files <- c(
                        paste0("../data/Chargaff_Equilibrium/", 
                               "ChargaffEquilibriumDistribution_scaling_one.Rdata"), 
                        files
                    )
                } else {
                    files <- list.files(
                        path = "../data/Chargaff_Equilibrium",
                        pattern = "^ChargaffEquilibriumDistribution_scaling.*\\.Rdata",
                        full.names = TRUE
                    )
                    files <- files[c(2,4,1,3)]
                    scaling <- c("one" = 1, "two" = 2, "five" = 5, "ten" = 10)
                    file.name <- as.character(unname(scaling))
                }

                data.sets <- lapply(1:length(files), function(x){
                    res <- readRDS(file = files[x])
                    if(other_species){
                        scaling.col <- file.name[x]
                    } else {
                        scaling.col <- paste0("scaling_", file.name[x])
                    }
                    res <- dplyr::mutate(res, scalings = scaling.col)
                    return(res)
                })
                data.sets <- do.call(rbind, data.sets)
                if(other_species){
                    level.labels <- file.name
                } else {
                    level.labels <- paste0("scaling_", file.name)
                }
                self$sim_run <- as_tibble(data.sets) %>% 
                    dplyr::mutate(scalings = factor(
                            scalings, 
                            levels = level.labels)) %>% 
                    tidyr::gather(key, value, -scalings)
            } else {
                self$sim_run <- readRDS(
                    file = paste0("../data/Main_Simulation/", private$muttype, 
                                  "-", private$distribution, ifelse(
                                    private$random_init_bases, 
                                    "/random_init_bases_", ifelse(
                                        private$edge_skew_cases,
                                        "/edge_skew_cases_",
                                    "/")),
                                  private$muttype,"-", private$distribution, 
                                  "-scaling-",private$scale_fac,".Rdata")
                )
                self$sim_run$nC <- (self$sim_run$GC/100)/(1+self$sim_run$GCratio)
                self$sim_run$nG <- (self$sim_run$GC/100)-self$sim_run$nC
                self$sim_run$nT <- (1-(self$sim_run$GC/100)) / (1+self$sim_run$ATratio)
                self$sim_run$nA <- (1-(self$sim_run$GC/100)) - self$sim_run$nT            

                # chargaff equilibrium values from simulations
                Fluc.tol <- readRDS(file = "../data/Chargaff_Equilibrium/ChargaffEquilibrium.Rdata")
                private$EQtolerance <- mean(Fluc.tol$Mean, na.rm = TRUE)*((1/100)*25)

                # get PR-2 compliance values for each kingdom
                private$species_tolerance$prokaryotes <- read.csv(
                    file = paste0("../../01-genome_composition/data/01-Prokaryotes/", 
                                  "PR2_compliance/PR2_fluctuations.csv"),
                    header = TRUE
                )
                private$species_tolerance$eukaryotes <- read.csv(
                    file = paste0("../../01-genome_composition/data/02-Eukaryotes/", 
                                  "PR2_compliance/PR2_fluctuations.csv"),
                    header = TRUE
                )
                private$species_tolerance$viruses <- read.csv(
                    file = paste0("../../01-genome_composition/data/03-Viruses/", 
                                  "PR2_compliance/PR2_fluctuations.csv"),
                    header = TRUE
                )
            }
            
            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")                          
        },

        #' @description
        #' Generates all plots for the equilibration runs.
        #' @return None.
        plot_equilibrations = function(other_species){
            t1 <- Sys.time()
            cur.msg <- "Generating equilibration plots"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            p1 <- self$sim_run %>% 
                dplyr::filter(key != "Difference") %>% 
                ggplot(aes(x = value, fill = key)) + 
                geom_histogram(
                    color = "#e9ecef", 
                    alpha = 0.5, 
                    position = 'identity', 
                    bins = 40, 
                    show.legend = FALSE) +
                geom_vline(
                    xintercept = 4.28, 
                    linetype = "dashed") +
                facet_wrap(vars(scalings), ncol = 1, scales = "free_y") + 
                scale_fill_manual(values = c("#69b3a2", "#404080")) + 
                coord_cartesian(xlim = c(0, 10)) + 
                labs(
                    x = "Years (Billion)",
                    y = ""
                )

            p2 <- self$sim_run %>% 
                dplyr::filter(key == "Difference") %>% 
                ggplot(aes(x = value)) + 
                geom_histogram(
                    color = "#e9ecef", 
                    alpha = 1, 
                    position = 'identity', 
                    bins = 70, 
                    show.legend = FALSE) +
                facet_wrap(vars(scalings), ncol = 1, scales = "free_y") + 
                scale_fill_manual(values = c("#69b3a2")) + 
                coord_cartesian(xlim = c(-3, 3)) + 
                labs(
                    x = "Difference in years (Billion)",
                    y = ""
                )

            dir.create(
                path = "../figures/Chargaff_Equilibrium/",
                showWarnings = FALSE,
                recursive = TRUE
            )

            save_plot_name = paste0("../figures/Chargaff_Equilibrium/ChargaffEquilibrium", 
                ifelse(other_species, "_OtherSpecies.", "Distribution."), private$save_as)
            if(private$save_as == "pdf"){
                pdf(width = 15, height = 16, file = save_plot_name)
            } else if(private$save_as == "png"){
                png(width = 1000, height = 700, file = save_plot_name)
            }
            gridExtra::grid.arrange(p1, p2, ncol = 2)
            pic.saved <- dev.off()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Generate GC content as histograms.
        #' @return None.
        plot_GC_hist = function(){
            t1 <- Sys.time()
            cur.msg <- "Generating GC content histograms"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            p1 <- self$sim_run %>%
                as_tibble() %>%
                ggplot(aes(
                    x = GC,
                    y = after_stat(density))) +
                geom_histogram(
                    fill = "skyblue",
                    color = "black",
                    alpha = 1,
                    breaks = seq(min(self$sim_run$GC, na.rm = TRUE),
                                 max(self$sim_run$GC, na.rm = TRUE),
                                 length.out = 51)) +
                coord_cartesian(xlim = c(0, 100)) + 
                labs(x = paste0(
                    "GC hist (%)", "\n",
                    "Mean = ", signif(mean(self$sim_run$GC, na.rm = TRUE), 2), " ",
                    "St.Dev = ", signif(sd(self$sim_run$GC, na.rm = TRUE), 2)),
                    y = "Density",
                    title = paste0("Scaling: ", private$scale_fac)
                )

            dir.create(
                path = paste0("../figures/Main_Simulation/",
                              private$muttype, "-", private$distribution,
                              "/Scaling_", private$scale_fac),
                showWarnings = FALSE,
                recursive = TRUE
            )
            ggsave(
                filename = paste0("../figures/Main_Simulation/",
                                  private$muttype, "-", private$distribution,
                                  "/Scaling_", private$scale_fac,
                                  ifelse(private$random_init_bases,
                                    "/random_init_bases_", ifelse(
                                        private$edge_skew_cases,
                                        "/edge_skew_cases_",
                                    "/")),
                                  "GC_hist.", private$save_as),
                plot = p1,
                width=15, height=8
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")            
        },

        #' @description
        #' Generate AT content as histograms.
        #' @return None.
        plot_AT_hist = function(){
            t1 <- Sys.time()
            cur.msg <- "Generating AT content histograms"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            p1 <- self$sim_run %>%
                as_tibble() %>%
                ggplot(aes(
                    x = AT,
                    y = after_stat(density))) +
                geom_histogram(
                    fill = "skyblue",
                    color = "black",
                    alpha = 1,
                    breaks = seq(min(self$sim_run$AT, na.rm = TRUE),
                                 max(self$sim_run$AT, na.rm = TRUE),
                                 length.out = 51)) +
                coord_cartesian(xlim = c(0, 100)) + 
                labs(x = paste0(
                    "AT hist (%)", "\n",
                    "Mean = ", signif(mean(self$sim_run$AT, na.rm = TRUE), 2), " ",
                    "St.Dev = ", signif(sd(self$sim_run$AT, na.rm = TRUE), 2)),
                    y = "Density",
                    title = paste0("Scaling: ", private$scale_fac)
                )

            ggsave(
                filename = paste0("../figures/Main_Simulation/",
                                  private$muttype, "-", private$distribution,
                                  "/Scaling_", private$scale_fac,
                                  ifelse(private$random_init_bases,
                                    "/random_init_bases_", ifelse(
                                        private$edge_skew_cases,
                                        "/edge_skew_cases_",
                                    "/")),
                                  "AT_hist.", private$save_as),
                plot = p1,
                width=15, height=8
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")            
        },

        #' @description
        #' Plot the G-C vs. A-T differences.
        #' @return None.
        plot_difference = function(){
            t1 <- Sys.time()
            cur.msg <- "Generating density plot for G-C vs. A-T"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            save_plot_name <- paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                ifelse(private$random_init_bases,
                "/random_init_bases_", ifelse(
                    private$edge_skew_cases,
                    "/edge_skew_cases_",
                "/")),
                "G-C_vs_A-T.", private$save_as
            )
            x = self$sim_run$nA-self$sim_run$nT
            y = self$sim_run$nG-self$sim_run$nC

            if(private$save_as == "pdf"){
                pdf(width = 6, height = 6, file = save_plot_name)
            } else if(private$save_as == "png"){
                png(width = 550, height = 550, file = save_plot_name)
            }         
            smoothScatter(
                y=y, x=x,
                nrpoints=100, nbin=1000,
                bandwidth=c(
                    diff(range(x, na.rm = TRUE))/500, 
                    diff(range(y, na.rm = TRUE))/500
                ),
                xlim=c(-1,1),ylim=c(-1,1),
                xlab="A-T", ylab="G-C", 
                main=paste0("Scaling: ", private$scale_fac),
                colramp=colorRampPalette(c(
                    "white","blue","skyblue",
                    "chartreuse3","green","yellow",
                    "orange","red","darkred"
                )),
                cex.axis = 1.2, cex.lab = 1.1
            )
            abline(v = 0, lty = 2, lwd = 2)
            abline(h = 0, lty = 2, lwd = 2)  
            pic.saved <- dev.off()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")              
        },

        #' @description
        #' Plot the GC-ratio vs. AT-ratio
        #' @return None.
        plot_ratios = function(){
            t1 <- Sys.time()
            cur.msg <- "Generating density plot for GC ratio vs. AT ratio"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            plot_func = function(xlim=c(0,100), ylim=c(0,100)){                
                smoothScatter(
                    y=self$sim_run$GCratio, 
                    x=self$sim_run$ATratio,
                    nrpoints=100, nbin=1000,
                    xlim=xlim,ylim=ylim,
                    xlab="AT-ratio", ylab="GC-ratio", 
                    main=paste0("Scaling: ", private$scale_fac),
                    colramp=colorRampPalette(c(
                        "white","blue","skyblue",
                        "chartreuse3","green","yellow",
                        "orange","red","darkred"
                    )),
                    cex.axis = 1.2, cex.lab = 1.1
                )
                abline(v = 1, lty = 2, lwd = 2)
                abline(h = 1, lty = 2, lwd = 2)              
            }

            save_plot_name_full <- paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                ifelse(private$random_init_bases,
                "/random_init_bases_", ifelse(
                    private$edge_skew_cases,
                    "/edge_skew_cases_",
                "/")),
                "GC-ratio_vs_AT-ratio.", private$save_as
            )
            save_plot_name_zoomed <- paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                ifelse(private$random_init_bases,
                "/random_init_bases_", ifelse(
                    private$edge_skew_cases,
                    "/edge_skew_cases_",
                "/")),
                "GC-ratio_vs_AT-ratio_0-10-range.", private$save_as
            )

            if(private$save_as == "pdf"){
                pdf(width = 15, height = 8, file = save_plot_name_full)
                plot_func()
                pic.saved <- dev.off()

                pdf(width = 15, height = 8, file = save_plot_name_zoomed)
                plot_func(xlim=c(0,10), ylim=c(0,10))
                pic.saved <- dev.off()                
            } else if(private$save_as == "png"){
                png(width = 850, height = 850, file = save_plot_name_full)
                plot_func()
                pic.saved <- dev.off()

                png(width = 550, height = 550, file = save_plot_name_zoomed)
                plot_func(xlim=c(0,10), ylim=c(0,10))
                pic.saved <- dev.off()  
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")  
        },

        #' @description
        #' Plot the GC-skew vs. AT-skew
        #' @return None.
        plot_skews = function(){
            t1 <- Sys.time()
            cur.msg <- "Generating density plot for GC skew vs. AT skew"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            save_plot_name <- paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                ifelse(private$random_init_bases,
                "/random_init_bases_", ifelse(
                    private$edge_skew_cases,
                    "/edge_skew_cases_",
                "/")),
                "GC-skew_vs_AT-skew.", private$save_as
            )    

            if(private$save_as == "pdf"){
                pdf(width = 6, height = 6, file = save_plot_name)
            } else if(private$save_as == "png"){
                png(width = 550, height = 550, file = save_plot_name)
            }
            smoothScatter(
                y=self$sim_run$GCskew, 
                x=self$sim_run$ATskew,
                nrpoints=100, nbin=1000,
                xlab="AT skew", ylab="GC skew", 
                main=paste0("Scaling: ", private$scale_fac),
                colramp=colorRampPalette(c(
                    "white","blue","skyblue",
                    "chartreuse3","green","yellow",
                    "orange","red","darkred"
                )),
                cex.axis = 1.2, cex.lab = 1.1
            )
            abline(v = 0, lty = 2, lwd = 2)
            abline(h = 0, lty = 2, lwd = 2)
            pic.saved <- dev.off()  

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")                        
        },

        #' @description
        #' Plot the G+C content vs. mutation rate constants
        #' @param tolerance Boolean. If TRUE, will apply PR-2 tolerance.
        #' @return None.
        plot_gc_content_vs_rates = function(tolerance = FALSE){
            t1 <- Sys.time()
            cur.msg <- "Generating density plot for GC content vs. rates"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            #' @description
            #' Helper function for generating plots.
            #' @param species Character vector of c("Prokaryotes", "Eukaryotes", "Viruses")
            #' @param tolerance Boolean. If TRUE, will apply PR-2 tolerance.
            #' @return List of plots.
            plot.gc.cont.func <- function(species = NULL, tolerance = FALSE){
                rates <- c("kag", "kat", "kac", "kga", "kgt", "kcg")
                xlab <- switch(private$muttype,
                    "Strand_Symmetric" = paste(rates, ", byr-1", sep = ""),
                    "Non_Symmetric" = rates
                )

                if(tolerance){
                    species_tol_values <- private$species_tolerance[[species]]
                    at.tol <- species_tol_values[species_tol_values$metadata == "AT_skew", "st.dev"]
                    gc.tol <- species_tol_values[species_tol_values$metadata == "GC_skew", "st.dev"]
                    ind <- which(
                        abs(self$sim_run$GCskew-0) <= gc.tol & abs(self$sim_run$ATskew-0) <= at.tol
                    )

                    plots <- lapply(1:length(rates), function(x){
                        x.values <- self$sim_run[ind, rates[x]]
                        y.values <- self$sim_run$GC[ind]
                        return(smoothScatter(
                            x=x.values, y=y.values, 
                            nrpoints=100, 
                            nbin=400,
                            bandwidth=c(diff(range(x.values, na.rm = TRUE))/200, 
                                        diff(range(self$sim_run$GC, na.rm = TRUE))/200),
                            xlab=xlab[x],
                            ylab="G+C content, %",
                            main=paste0("Scaling: ", private$scale_fac),
                            colramp=colorRampPalette(c(
                                "white","blue","skyblue",
                                "chartreuse3","green","yellow",
                                "orange","red","darkred"
                            )), 
                            ylim = c(0,100),
                            cex.axis = 1.5, cex.lab = 1.5
                        ))
                    })
                } else {
                    plots <- lapply(1:length(rates), function(x){
                        x.values <- self$sim_run[,rates[x]] 
                        return(smoothScatter(
                            x=x.values, y=self$sim_run$GC, 
                            nrpoints=100, 
                            nbin=400,
                            bandwidth=c(diff(range(x.values, na.rm = TRUE))/200, 
                                        diff(range(self$sim_run$GC, na.rm = TRUE))/200),
                            xlab=xlab[x],
                            ylab="G+C content, %",
                            main=paste0("Scaling: ", private$scale_fac),
                            colramp=colorRampPalette(c(
                                "white","blue","skyblue",
                                "chartreuse3","green","yellow",
                                "orange","red","darkred"
                            )),
                            ylim = c(0,100),
                            cex.axis = 1.5, cex.lab = 1.5
                        ))
                    })
                }
                return(plots)
            }

            save_plot_name <- paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                ifelse(tolerance, 
                    ifelse(private$random_init_bases,
                    "/random_init_bases_tolerance-GC-content_vs_rates-all.", 
                    ifelse(private$edge_skew_cases, 
                    "/edge_skew_cases_tolerance-GC-content_vs_rates-all.",
                    "/tolerance-GC-content_vs_rates-all.")),
                ifelse(private$random_init_bases, 
                    "/random_init_bases_GC-content_vs_rates-all.",
                    ifelse(private$edge_skew_cases, 
                    "/edge_skew_cases_GC-content_vs_rates-all.",
                "/GC-content_vs_rates-all."))),
                private$save_as
            )

            if(private$save_as == "pdf"){
                pdf(width = 6, height = 15, file = save_plot_name)
            } else if(private$save_as == "png"){
                png.width <- ifelse(tolerance, 1200, 500)
                png(width = png.width, height = 1100, file = save_plot_name)
            }
            if(tolerance){
                par(mfcol=c(6,3))
                plot.gc.cont.func(species = "eukaryotes", tolerance = tolerance)
                plot.gc.cont.func(species = "prokaryotes", tolerance = tolerance)
                plot.gc.cont.func(species = "viruses", tolerance = tolerance)
            } else {
                par(mfcol=c(6,1))
                plot.gc.cont.func()
            }
            pic.saved <- dev.off()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")                    
        },

        #' @description
        #' Plot the G+C skew vs. mutation rate constants
        #' @param tolerance Boolean. If TRUE, will apply PR-2 tolerance.
        #' @return None.
        plot_gc_skew_vs_rates = function(tolerance = FALSE){
            t1 <- Sys.time()
            cur.msg <- "Generating density plot for GC skew vs. rates"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            #' @description
            #' Helper function for generating plots.
            #' @param species Character vector of c("Prokaryotes", "Eukaryotes", "Viruses")
            #' @param tolerance Boolean. If TRUE, will apply PR-2 tolerance.
            #' @return List of plots.
            plot.gc.skew.func <- function(species = NULL, tolerance = FALSE){
                rates <- c("kag", "kat", "kac", "ktc", "kgt", "kgc",  
                           "kga", "kta", "kca", "kct", "ktg", "kcg")
                rates.lims <- c(
                    ceiling(max(self$sim_run[, "kag"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kat"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kac"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kga"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kgt"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kgc"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kag"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kat"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kac"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kga"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kgt"], na.rm = TRUE)),
                    ceiling(max(self$sim_run[, "kgc"], na.rm = TRUE))                  
                )                           
                xlab <- switch(private$muttype,
                    "Strand_Symmetric" = paste(rates, ", byr-1", sep = ""),
                    "Non_Symmetric" = rates
                )

                if(tolerance){
                    species_tol_values <- private$species_tolerance[[species]]
                    at.tol <- species_tol_values[species_tol_values$metadata == "AT_skew", "st.dev"]
                    gc.tol <- species_tol_values[species_tol_values$metadata == "GC_skew", "st.dev"]
                    ind <- which(
                        abs(self$sim_run$GCskew-0) <= gc.tol & abs(self$sim_run$ATskew-0) <= at.tol
                    )

                    plots <- lapply(1:length(rates), function(x){
                        x.values <- self$sim_run[ind, rates[x]]
                        y.values <- self$sim_run$GCskew[ind]
                        smoothScatter(
                            x=x.values, y=y.values, 
                            nrpoints=100, 
                            nbin=400,
                            bandwidth=c(diff(range(x.values, na.rm = TRUE))/200, 
                                        diff(range(self$sim_run$GC, na.rm = TRUE))/200),
                            xlim=c(0, ceiling(max(x.values, na.rm=TRUE))),
                            xlab=xlab[x],
                            ylab="GC skew",
                            main=paste0("Scaling: ", private$scale_fac),
                            colramp=colorRampPalette(c(
                                "white","blue","skyblue",
                                "chartreuse3","green","yellow",
                                "orange","red","darkred"
                            )), 
                            cex.axis = 1.5, cex.lab = 1.5
                        )
                        abline(h = 0, lty = 2, lwd = 2)
                    })
                } else {
                    plots <- lapply(1:length(rates), function(x){
                        x.values <- self$sim_run[,rates[x]] 
                        smoothScatter(
                            x=x.values, y=self$sim_run$GCskew, 
                            nrpoints=100, 
                            nbin=400,
                            bandwidth=c(diff(range(x.values, na.rm = TRUE))/200, 
                                        diff(range(self$sim_run$GCskew, na.rm = TRUE))/200),
                            xlim=c(0, rates.lims[x]),                                       
                            xlab=xlab[x],
                            ylab="GC skew",
                            main=paste0("Scaling: ", private$scale_fac),
                            colramp=colorRampPalette(c(
                                "white","blue","skyblue",
                                "chartreuse3","green","yellow",
                                "orange","red","darkred"
                            )),
                            cex.axis = 1.5, cex.lab = 1.5
                        )
                        abline(h = 0, lty = 2, lwd = 2)
                    })
                }
                return(plots)
            }

            save_plot_name <- paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                ifelse(tolerance, 
                    ifelse(private$random_init_bases,
                    "/random_init_bases_tolerance-GC-skew_vs_rates-all.", 
                    ifelse(private$edge_skew_cases, 
                    "/edge_skew_cases_tolerance-GC-skew_vs_rates-all.",
                    "/tolerance-GC-skew_vs_rates-all.")),
                ifelse(private$random_init_bases, 
                    "/random_init_bases_GC-skew_vs_rates-all.",
                    ifelse(private$edge_skew_cases, 
                    "/edge_skew_cases_GC-skew_vs_rates-all.", 
                "/GC-skew_vs_rates-all."))),
                private$save_as
            )

            if(private$save_as == "pdf"){
                pdf(width = 6, height = 15, file = save_plot_name)
            } else if(private$save_as == "png"){
                png.width <- ifelse(tolerance, 1200, 700)
                png(width = png.width, height = 1700, file = save_plot_name)
            }
            if(tolerance){
                par(mfcol=c(12,3))
                plot.gc.skew.func(species = "eukaryotes", tolerance = tolerance)
                plot.gc.skew.func(species = "prokaryotes", tolerance = tolerance)
                plot.gc.skew.func(species = "viruses", tolerance = tolerance)
            } else {
                par(mfcol=c(6,2))
                plot.gc.skew.func()
            }
            pic.saved <- dev.off()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")                    
        },        

        #' @description
        #' Plot the mutation rate constants as ratio distribution plots
        #' @return None.
        plot_ratio_constants = function(){
            #' @description
            #' Helper function for generating plots.
            #' @param species Character vector of c("Prokaryotes", "Eukaryotes", "Viruses")
            #' @param tolerance Boolean. If TRUE, will apply PR-2 tolerance.
            #' @return List of plots.
            plot.ratio.func <- function(species = NULL, tolerance = FALSE){
                length.out <- 10*5+1
                df.filtered <- as_tibble(self$sim_run)

                # chargaff compliance
                if(tolerance){
                    species_tol_values <- private$species_tolerance[[species]]
                    at.tol <- species_tol_values[species_tol_values$metadata == "AT_skew", "st.dev"]
                    gc.tol <- species_tol_values[species_tol_values$metadata == "GC_skew", "st.dev"]
                    ind <- which(
                        abs(self$sim_run$GCskew-0) <= gc.tol & abs(self$sim_run$ATskew-0) <= at.tol
                    )

                    df.filtered <- df.filtered %>% 
                        dplyr::slice(ind)
                }

                df.filtered <- df.filtered %>% 
                    dplyr::summarise(
                        `kag/ktc` = kag/ktc,
                        `kat/kta` = kat/kta,
                        `kac/ktg` = kac/ktg,
                        `kct/kga` = kct/kga,
                        `kca/kgt` = kca/kgt,
                        `kcg/kgc` = kcg/kgc) %>% 
                    tidyr::gather(key, frac) %>% 
                    dplyr::filter(frac <= 5)

                if(nrow(df.filtered) > 1){
                    p1 <- df.filtered %>% 
                        ggplot(aes(x = frac, after_stat(density))) + 
                        geom_histogram(
                            fill = "skyblue", 
                            color = "black",
                            alpha = 1,
                            breaks = seq(min(df.filtered$frac, na.rm = TRUE), 
                                        max(df.filtered$frac, na.rm = TRUE), 
                                        length.out = length.out)) +
                        geom_vline(
                            xintercept = 1,
                            linetype = "dashed",
                            linewidth = 1) + 
                        facet_wrap(~key, ncol = 1) + 
                        coord_cartesian(xlim = c(0, 5)) + 
                        scale_fill_manual(values = c("#69b3a2")) + 
                        labs(
                            x = "Rate constant ratio",
                            y = "Density"
                        )
                } else {
                    p1 <- NULL
                }
                return(p1)
            }

            for(tolerance in c(FALSE, TRUE)){
                t1 <- Sys.time()
                cur.msg <- paste0("Generating ratio distribution plots with tolerance: ", tolerance)
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                save_plot_name = paste0(
                    "../figures/Main_Simulation/",
                    private$muttype, "-", private$distribution,
                    "/Scaling_", private$scale_fac,
                    ifelse(tolerance, 
                        ifelse(private$random_init_bases,
                        "/random_init_bases_tolerance-rate_ratios.", 
                        ifelse(private$edge_skew_cases, 
                        "/edge_skew_cases_tolerance-rate_ratios.",
                        "/tolerance-rate_ratios.")),
                    ifelse(private$random_init_bases, 
                        "/random_init_bases_allind-rate_ratios.",
                        ifelse(private$edge_skew_cases, 
                        "/edge_skew_cases_allind-rate_ratios.", 
                    "/allind-rate_ratios."))), 
                    private$save_as
                )

                if(tolerance){
                    p1 <- plot.ratio.func(species = "eukaryotes", tolerance = TRUE)
                    p2 <- plot.ratio.func(species = "prokaryotes", tolerance = TRUE)
                    p3 <- plot.ratio.func(species = "viruses", tolerance = TRUE)

                    if(private$save_as == "pdf"){
                        pdf(width = 18, height = 14, file = save_plot_name)
                    } else if(private$save_as == "png"){
                        png(width = 1000, height = 700, file = save_plot_name)
                    }
                    gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
                    pic.saved <- dev.off()
                } else {
                    if(private$save_as == "pdf"){
                        pdf(width = 9, height = 14, file = save_plot_name)
                    } else if(private$save_as == "png"){
                        png(width = 700, height = 1000, file = save_plot_name)
                    }
                    p1 <- plot.ratio.func(tolerance = FALSE)
                    gridExtra::grid.arrange(p1, ncol = 1)
                    pic.saved <- dev.off()
                }

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }
        },

        #' @description
        #' Plot the mutation rate constants as ratio distribution plots
        #' @return None.
        plot_ratio_constants_for_strand_symmetric_cases = function(){
            #' @description
            #' Helper function for generating plots.
            #' @return List of plots.
            plot.ratio.func <- function(){
                length.out <- 10*5+1
                df.filtered <- as_tibble(self$sim_run)
                df.filtered <- df.filtered %>% 
                    dplyr::summarise(
                        m_n = kct/kag,
                        i_j = kca/kac) %>% 
                    tidyr::gather(key, frac) %>% 
                    dplyr::filter(frac <= 5)

                if(nrow(df.filtered) > 1){
                    p1 <- df.filtered %>% 
                        ggplot(aes(x = frac, after_stat(density))) + 
                        geom_histogram(
                            fill = "skyblue", 
                            color = "black",
                            alpha = 1,
                            breaks = seq(min(df.filtered$frac, na.rm = TRUE), 
                                         max(df.filtered$frac, na.rm = TRUE), 
                                         length.out = length.out)) +
                        geom_vline(
                            xintercept = 1,
                            linetype = "dashed",
                            linewidth = 1) + 
                        facet_wrap(~key, ncol = 1, scales = "free_y") + 
                        coord_cartesian(xlim = c(0, 5)) + 
                        scale_fill_manual(values = c("#69b3a2")) + 
                        labs(
                            x = "Rate constant ratio",
                            y = "Density"
                        )
                } else {
                    p1 <- NULL
                }
                return(p1)
            }

            t1 <- Sys.time()
            cur.msg <- "Generating ratio distribution plots"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            save_plot_name = paste0(
                "../figures/Main_Simulation/",
                private$muttype, "-", private$distribution,
                "/Scaling_", private$scale_fac,
                "/allind-rate_ratios.", 
                private$save_as
            )

            if(private$save_as == "pdf"){
                pdf(width = 9, height = 14, file = save_plot_name)
            } else if(private$save_as == "png"){
                png(width = 700, height = 1000, file = save_plot_name)
            }
            p1 <- plot.ratio.func()
            gridExtra::grid.arrange(p1, ncol = 1)
            pic.saved <- dev.off()

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")
        },

        #' @description
        #' Plots the evolution of the GC-AT skew in intervals of 1 billion years
        #' for the strand-symmetric mutation rate simulation outcomes.
        #' @return None.
        plot_skew_evolution = function(){
            # edit simulation data set
            to.keep <- which(!is.na(stringr::str_extract(
                string = names(self$sim_run), 
                pattern = "skew"
            )))
            self$sim_run <- self$sim_run[,to.keep]
            self$sim_run <- self$sim_run[, c(3:10, 1:2)]

            # extract year timelines
            bn.years <- stringr::str_extract(
                string = names(self$sim_run), 
                pattern = ".+(?=.bn.)"
            )
            bn.years <- unique(bn.years)
            bn.years <- c(bn.years[1:(length(bn.years)-1)], "four_28")

            first.col <- seq(from = 1, to = length(self$sim_run), by = 2)
            second.col <- seq(from = 2, to = length(self$sim_run), by = 2)

            #' @description
            #' Plot the GC vs. AT skews.
            #' @return Plot.
            plot_func <- function(){
                smoothScatter(
                    x=self$sim_run[,second.col[i]],
                    y=self$sim_run[,first.col[i]], 
                    nrpoints=100, nbin=1000,
                    xlab="AT skew", ylab="GC skew", 
                    main=paste0("Scaling: ", private$scale_fac),
                    colramp=colorRampPalette(c(
                        "white","blue","skyblue",
                        "chartreuse3","green","yellow",
                        "orange","red","darkred"
                    )),
                    xlim=c(-1e-15,1e-15),
                    ylim=c(-1e-15,1e-15),
                    cex.axis = 1.2, cex.lab = 1.1
                )
                abline(v = 0, lty = 2, lwd = 2)
                abline(h = 0, lty = 2, lwd = 2)
            }

            for(i in 1:length(bn.years)){            
                png(width=600, height=600, 
                    paste0("../figures/Main_Simulation/Strand_Symmetric-normal/scaling_",
                            private$scale_fac,"/evolution_GC-skew_vs_AT-skew_", 
                            bn.years[i], "bn.png"))
                plot_func()
                pic.saved <- dev.off()
            }

            # make animation movie out of strand symmetry evolution images
            for(package in c("magick", "animation")){
                if(!package %in% rownames(installed.packages())){
                    install.packages(package)
                    suppressPackageStartupMessages(suppressWarnings(
                        library(package, character.only = TRUE)
                    ))
                }
            }

            # load photos for video clip
            SLF.photos <- list.files(
                path = paste0("../../figures/Main_Simulation/", ,
                              "Strand_Symmetric-normal/scaling_", 
                              private$scale_fac),
                pattern = "evolution", 
                all.files = TRUE, 
                ignore.case = TRUE,
                full.names = TRUE
            )
            SLF.photos <- SLF.photos[c(3,5,4,2,1)]

            #This extracts the underlying height, width, and type of image.
            img.height <- magick::image_info(image_read(SLF.photos[1]))$height
            img.width <- magick::image_info(image_read(SLF.photos[1]))$width
            img.type <- magick::image_info(image_read(SLF.photos[1]))$format

            # display each picture for 0.25 seconds
            animation::ani.options(
                interval = 0.25,
                ani.height = img.height,
                ani.width = img.width,
                ani.dev = tolower(img.type),
                ani.type = tolower(img.type)
            )

            # increasing video dimensions for better image quality
            opts <- paste("-s ", img.height * 1.5, "x", img.width * 1.5, sep = "")

            animation::saveVideo(
                for(i in 1:length(SLF.photos)){
                    SLF.image <- magick::image_read(SLF.photos[i])
                    plot(SLF.image)
                },
                video.name = paste0("../figures/Main_Simulation/", 
                                    "Strand_Symmetric-normal/scaling_",
                                    private$scale_fac, "/SkewEvolutionAnimation.mp4"),
                other.opts = paste0("-pix_fmt yuv420p -b 300k ", opts)
            )
        },

        #' @description
        #' Plot the mutation rate constants as ratio distribution plots.
        #' @return None.
        generate_equil_const_plots = function(){
            #' @description
            #' Helper function for generating plots.
            #' @param species Character vector of c("Prokaryotes", "Eukaryotes", "Viruses")
            #' @param tolerance Boolean. If TRUE, will apply PR-2 tolerance.
            #' @return List of plots.
            plot.equil.func <- function(species = NULL, tolerance = FALSE){
                length.out <- 10*5+1
                df.filtered <- as_tibble(self$sim_run)

                # chargaff compliance
                if(tolerance){
                    species_tol_values <- private$species_tolerance[[species]]
                    at.tol <- species_tol_values[species_tol_values$metadata == "AT_skew", "st.dev"]
                    gc.tol <- species_tol_values[species_tol_values$metadata == "GC_skew", "st.dev"]
                    ind <- which(
                        abs(self$sim_run$GCskew-0) <= gc.tol & abs(self$sim_run$ATskew-0) <= at.tol
                    )

                    df.filtered <- df.filtered %>% 
                        dplyr::slice(ind)
                }

                df.filtered <- df.filtered %>% 
                    dplyr::summarise(
                        `kct/ktc (m/n)` = kct/ktc,
                        `kag/kga (m/n)` = kag/kga,
                        `kca/kac (i/j)` = kca/kac,
                        `kgt/ktg (i/j)` = kgt/ktg,
                        `kcg/kgc (k/k)` = kcg/kgc,
                        `kat/kta (l/l)` = kat/kta) %>% 
                    tidyr::gather(key, frac) %>% 
                    dplyr::filter(frac <= 5)

                if(nrow(df.filtered) > 1){
                    p1 <- df.filtered %>% 
                        ggplot(aes(x = frac, after_stat(density))) + 
                        geom_histogram(
                            fill = "skyblue", 
                            color = "black",
                            alpha = 1,
                            breaks = seq(min(df.filtered$frac, na.rm = TRUE), 
                                        max(df.filtered$frac, na.rm = TRUE), 
                                        length.out = length.out)) +
                        geom_vline(
                            xintercept = 1,
                            linetype = "dashed",
                            linewidth = 1) + 
                        facet_wrap(~key, ncol = 1) + 
                        coord_cartesian(xlim = c(0, 5)) + 
                        scale_fill_manual(values = c("#69b3a2")) + 
                        labs(
                            x = "Rate constant ratio",
                            y = "Density"
                        )
                } else {
                    p1 <- NULL
                }
                return(p1)
            }

            for(tolerance in c(FALSE, TRUE)){
                t1 <- Sys.time()
                cur.msg <- paste0("Generating equil constant plots with tolerance: ", tolerance)
                l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
                cat(paste0(cur.msg, l))

                save_plot_name = paste0(
                    "../figures/Main_Simulation/",
                    private$muttype, "-", private$distribution,
                    "/Scaling_", private$scale_fac,
                    ifelse(tolerance, 
                        ifelse(private$random_init_bases,
                        "/random_init_bases_tolerance-rate_ratios_fwd-div-rev.", 
                        ifelse(private$edge_skew_cases, 
                        "/edge_skew_cases_tolerance-rate_ratios_fwd-div-rev.",
                        "/tolerance-rate_ratios_fwd-div-rev.")),
                    ifelse(private$random_init_bases, 
                        "/random_init_bases_allind-rate_ratios_fwd-div-rev.",
                        ifelse(private$edge_skew_cases, 
                        "/edge_skew_cases_allind-rate_ratios_fwd-div-rev.", 
                    "/allind-rate_ratios_fwd-div-rev."))), 
                    private$save_as
                )

                if(tolerance){
                    p1 <- plot.equil.func(species = "eukaryotes", tolerance = TRUE)
                    p2 <- plot.equil.func(species = "prokaryotes", tolerance = TRUE)
                    p3 <- plot.equil.func(species = "viruses", tolerance = TRUE)

                    if(private$save_as == "pdf"){
                        pdf(width = 18, height = 14, file = save_plot_name)
                    } else if(private$save_as == "png"){
                        png(width = 1000, height = 700, file = save_plot_name)
                    }
                    gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
                    pic.saved <- dev.off()
                } else {
                    if(private$save_as == "pdf"){
                        pdf(width = 9, height = 14, file = save_plot_name)
                    } else if(private$save_as == "png"){
                        png(width = 700, height = 1000, file = save_plot_name)
                    }
                    p1 <- plot.equil.func(tolerance = FALSE)
                    gridExtra::grid.arrange(p1, ncol = 1)
                    pic.saved <- dev.off()
                }

                total.time <- Sys.time() - t1
                cat("DONE! --", signif(total.time[[1]], 2), 
                    attr(total.time, "units"), "\n")
            }
        }
    )
)
