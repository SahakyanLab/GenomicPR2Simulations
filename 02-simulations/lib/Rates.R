Rates <- R6::R6Class(
    classname = "Rates",
    public = list(
        initialize = function(trek_scale){
            if(!missing(trek_scale)) private$trek_scale <- trek_scale
        },

        #' @description
        #' Generates all plots for the main simulation.
        #' @param save_as Character vector of c("png", "pdf") to save the plots.
        #' @return None
        process_rates = function(save_as = "png"){
            private$save_as <- save_as
            private$get_mut_rates()
            private$get_pr2_compliance_values()
            private$plot_lynch_rates()
            private$get_data()
            private$get_gc_content_curves()
        }
    ),
    private = list(
        #' @field trek_scale Boolean. If TRUE, will convert into trek scale.
        trek_scale = NULL,

        #' @field notes Data.frame of mutation rate constants from trek paper.
        notes = NULL,

        #' @field lynch_rates Data.frame of mutation rate constants of real life species.
        lynch_rates = NULL,

        #' @field species_tolerance Data.frame of the PR-2 compliance values for each kingdom.
        species_tolerance = NULL,

        #' @field species_gc Data.frame of the GC content for each species.
        species_gc = NULL, 

        #' @field save_as Character vector of c("png", "pdf") to save the plots.
        save_as = "png",

        #' @description
        #' Checks if a given URL exists or not
        #' @param url_in Character vector of the URL to download.
        #' @param t Numeric vector of max. time until timeout reached.
        #' @return Boolean
        valid_url = function(url_in, t = 1000){
            con <- url(url_in)
            check <- suppressWarnings(try(open.connection(
                con, open = "rt", timeout = t), silent = T)[1])
            suppressWarnings(try(close.connection(con), silent = T))
            return(ifelse(is.null(check), TRUE, FALSE))
        },

        #' @description
        #' Import the pre-calculated mutation rate constants from trek paper.
        #' @return None.
        get_mut_rates = function(){
            t1 <- Sys.time()
            cur.msg <- "Loading in trek-based mutation rate constants"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            private$notes$note_two <- read.csv(
                file = "../data/Raw/Trek-paper-Note-2-mutation-rates.csv",
                header = TRUE
            )

            # obtain mutation rates from Table 1 in the following paper
            # https://www.pnas.org/content/107/3/961
            private$lynch_rates <- read.csv(
                file = "../data/Raw/Michael_Lynch/Lynch-2010-mutation-rates_FREQUENCY.csv",
                header = TRUE
            )
            private$lynch_rates <- private$lynch_rates[match(
                private$notes$note_two$MUT, private$lynch_rates$MUT
            ),]
            # match rate constants from lynch with trek
            private$lynch_rates <- private$lynch_rates %>%
                dplyr::mutate(
                    MEAN = private$notes$note_two$MEAN, 
                    .after = 1
                )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n")  
        },

        #' @description
        #' Get PR-2 compliance values for each kingdom.
        #' @return None.
        get_pr2_compliance_values = function(){
            cur.msg <- "Loading in PR-2 compliance values"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            private$species_tolerance$prokaryotes <- read.csv(
                file = "../../01-genome_composition/data/01-Prokaryotes/PR2_compliance/PR2_fluctuations.csv",
                header = TRUE
            )
            private$species_tolerance$eukaryotes <- read.csv(
                file = "../../01-genome_composition/data/02-Eukaryotes/PR2_compliance/PR2_fluctuations.csv",
                header = TRUE
            )
            private$species_tolerance$viruses <- read.csv(
                file = "../../01-genome_composition/data/03-Viruses/PR2_compliance/PR2_fluctuations.csv",
                header = TRUE
            )
            private$species_gc <- read.csv(
                file = "../data/Raw/Michael_Lynch/GC_vs_Rates.csv",
                header = TRUE
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n") 
        },

        #' @description
        #' Plot mutation rate constants from lynchs' paper with trek 
        #' @return None.
        plot_lynch_rates = function(){
            t1 <- Sys.time()
            cur.msg <- "Plot mutation rate constants PNAS paper"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            lynch.plot <- private$lynch_rates %>%
                ggplot(aes(x = MEAN)) +
                geom_point(aes(y = H.sapiens, color = "H.sapiens")) +
                geom_line(aes(y = H.sapiens, color = "H.sapiens")) +
                geom_point(aes(y = D.melanogaster, color = "D.melanogaster")) +
                geom_point(aes(y = C.elegans, color = "C.elegans")) +
                geom_point(aes(y = A.thaliana, color = "A.thaliana")) +
                geom_point(aes(y = S.cerevisiae, color = "S.cerevisiae")) +
                geom_point(aes(y = E.coli, color = "E.coli")) +
                xlim(0, 1.2) + 
                ylim(0, 0.6) + 
                labs(
                    x = "Trek (Mean), byr",
                    y = "Lynch, 2010"
                )

            ggsave(
                filename = "../figures/Chargaff_Equilibrium/LynchRateFreq.pdf", 
                plot = lynch.plot
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n") 
        },

        #' @description
        #' Get fasta files of species and process them.
        #' @return None.
        get_data = function(){
            cur.msg <- "Analysing GC content and rate constants from species"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # conversion of rate frequencies to time bound rate constants
            Coef <- private$lynch_rates %>%
                lm(formula = MEAN ~ H.sapiens-1) %>%
                coefficients

            conversion.formula <- function(x){
                return(Coef*x)
            }

            final.rates <- private$lynch_rates %>%
                dplyr::select(1:3) %>%
                cbind(apply(
                    private$lynch_rates[4:length(private$lynch_rates)],
                    2, conversion.formula)) %>%
                dplyr::mutate(
                    MEDIAN = private$notes$note_two$MEDIAN,
                    SD = private$notes$note_two$SD,
                    .after = MEAN
                )

            # save converted rate constants as csv
            write.table(
                x = final.rates, 
                file = "../data/Raw/Michael_Lynch/Lynch-2010-converted-mutation-rates.csv", 
                sep = ",", 
                row.names = FALSE
            )

            # fasta files to download 
            species.names <- c(
                "Homo_sapiens", 
                "Arabidopsis_thaliana", 
                "Caenorhabditis_elegans", 
                "Drosophila_melanogaster",
                "Escherichia_coli", 
                "Saccharomyces_cerevisiae",
                "Photorhabdus_luminescens",
                "Teredinibacter_turnerae",
                "Mycobacterium_smegmatis",
                "Pseudomonas_fluorescens",
                "Rhodosporidium_toruloides",
                "Mus_musculus",
                "Daphnia_pulex",
                "Pristionchus_pacificus",
                "Daphnia_magna",
                "Pan_troglodytes",
                "Aotus_nancymaae"
            )

            download.files.url <- c(
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/", 
                       "405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/", 
                       "735/GCF_000001735.4_TAIR10.1/GCF_000001735.4_TAIR10.1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/", 
                       "985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/", 
                       "215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_", 
                       "plus_ISO1_MT_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/157/", 
                       "115/GCF_000157115.1_ASM15711v1/GCF_000157115.1_ASM15711v1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/756/", 
                       "235/GCA_000756235.1_ASM75623v1/GCA_000756235.1_ASM75623v1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/900/102/", 
                       "985/GCF_900102985.1_IMG-taxon_2597490348_annotated_assembly/GCF_900102985.1", 
                       "_IMG-taxon_2597490348_annotated_assembly_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/964/", 
                       "255/GCF_000964255.1_ASM96425v1/GCF_000964255.1_ASM96425v1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/013/349/", 
                       "145/GCF_013349145.1_ASM1334914v1/GCF_013349145.1_ASM1334914v1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/414/", 
                       "285/GCF_001414285.1_ATCC948-1/GCF_001414285.1_ATCC948-1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/320/", 
                       "785/GCF_000320785.1_RHOziaDV1.0/GCF_000320785.1_RHOziaDV1.0_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/", 
                       "635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/911/175/", 
                       "335/GCA_911175335.1_PA42_4.2/GCA_911175335.1_PA42_4.2_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/180/", 
                       "635/GCA_000180635.4_El_Paco_v._4/GCA_000180635.4_El_Paco_v._4_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/020/631/", 
                       "705/GCF_020631705.1_ASM2063170v1.1/GCF_020631705.1_ASM2063170v1.1_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/880/", 
                       "755/GCF_002880755.1_Clint_PTRv2/GCF_002880755.1_Clint_PTRv2_genomic.fna.gz"),
                paste0("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/952/", 
                       "055/GCF_000952055.2_Anan_2.0/GCF_000952055.2_Anan_2.0_genomic.fna.gz")
            )

            # big files need more time to download
            if(getOption('timeout') < 1000){
                options(timeout = 1000)
            }

            # download files
            for(i in 1:length(species.names)){
                # allow several attempts to donwload fasta files
                # as internet disruptions may cause failure
                attempt <- 1
                while(attempt <= 3){
                    url.test <- private$valid_url(url_in = download.files.url[i])
                    
                    if(isTRUE(url.test)){
                        try(download.file(
                            download.files.url[i], 
                            paste0("../data/Raw/Michael_Lynch/", i, "-", 
                                  species.names[i], ".dna.fna.gz")
                        ))
                        
                        # If file doesn't exist, it'll likely be a server problem so 
                        # allow it to retry up to 3 times
                        # If it requires >3 times, likely problems are:
                        #   (1) no internet connection
                        #   (2) another pattern that's unaccounted for
                        if(!file.exists(paste0("../data/Raw/Michael_Lynch/", i, "-", 
                                               species.names[i], ".dna.fna.gz"))){
                            Sys.sleep(3)
                            attempt <- attempt + 1
                    } else {
                        break
                    }
                    } 
                    if(isFALSE(url.test)){
                        try(download.file(download.files.url[i], 
                                        paste0("../data/Raw/Michael_Lynch/", i, "-", 
                                               species.names[i], ".dna.fna.gz")))
                    
                    if(!file.exists(paste0("../data/Raw/Michael_Lynch/", i, "-", 
                                           species.names[i], ".dna.fna.gz"))){
                        Sys.sleep(3)
                        attempt <- attempt + 1
                    } else { 
                        break
                    }
                    }
                }
            }

            # unzip all files if needed
            system("gunzip ../data/Raw/Michael_Lynch/*")

            # read fasta files
            files <- list.files(
                path = "../../data/Raw/Michael_Lynch", 
                pattern = ".fna$",
                full.names = TRUE
            )
            files <- stringr::str_sort(files, numeric = TRUE)

            # obtain GC content from species
            base.values <- lapply(files, function(x){
                fai <- readDNAStringSet(x)
                
                # Initialise data frame for base calculations
                genome_length <- width(fai)
                
                # extract base contents
                all.letters <- letterFrequency(fai, letters="ACGT", OR=0)
                
                # obtain GC average
                all.letters <- cbind(all.letters, genome_length)
                all.letters <- colSums(all.letters)
                G_plus_C  <- (all.letters["G"]+all.letters["C"])/all.letters["genome_length"]
                G_minus_C <- (all.letters["G"]-all.letters["C"])/all.letters["genome_length"]
                GC_skew <- G_minus_C/G_plus_C

                A_plus_T  <- (all.letters["A"]+all.letters["T"])/all.letters["genome_length"]
                A_minus_T <- (all.letters["A"]-all.letters["T"])/all.letters["genome_length"]
                AT_skew <- A_minus_T/A_plus_T

                return(list(G_plus_C, GC_skew, AT_skew))
            })
            GC.content <- sapply(base.values, `[[`, 1)
            GC.skew <- sapply(base.values, `[[`, 2)
            AT.skew <- sapply(base.values, `[[`, 3)

            other.private$species_gc.avg <- data.frame(
                species = files,
                GC.content = unlist(GC.content),
                GC.skew = unlist(GC.skew),
                AT.skew = unlist(AT.skew)
            )

            other.private$species_gc.avg %>% 
                dplyr::select(-AT.skew) %>% 
                write.csv(
                    file = "../data/Raw/Michael_Lynch/GC_values.csv",
                    row.names = TRUE
                )

            if(private$trek_scale){
                private$lynch_rates <- private$lynch_rates %>% 
                    dplyr::select(1:3) %>% 
                    cbind(apply(
                        private$lynch_rates[4:length(private$lynch_rates)], 
                        2, conversion.formula)) %>% 
                    as_tibble() %>% 
                    dplyr::mutate(
                        RATES = c("j", "n", "l", "i", "k", "m"),
                        .after = MUT
                    )
            } else {
                private$lynch_rates <- private$lynch_rates %>% 
                    as_tibble() %>% 
                    dplyr::mutate(
                        RATES = c("j", "n", "l", "i", "k", "m"),
                        .after = MUT
                    )
            }

            # calculate rate constant ratios
            ind.n <- which(private$lynch_rates$RATES == "n")
            ind.j <- which(private$lynch_rates$RATES == "j")
            ind.m <- which(private$lynch_rates$RATES == "m")
            ind.i <- which(private$lynch_rates$RATES == "i")

            # assign new GC contents
            rate.constant.ratios <- sapply(4:length(private$lynch_rates), function(x){
                return(
                    sum(private$lynch_rates[ind.n, x],private$lynch_rates[ind.j, x])/
                    sum(private$lynch_rates[ind.i, x],private$lynch_rates[ind.m, x])
                )
            })

            # assign ratios to new column 
            other.private$species_gc.avg <- data.frame(
                Species = colnames(private$lynch_rates)[4:length(private$lynch_rates)],
                Rates = rate.constant.ratios,
                GC.average = other.private$species_gc.avg$GC.content*100,
                GC.skew = other.private$species_gc.avg$GC.skew,
                AT.skew = other.private$species_gc.avg$AT.skew
            )

            # save new data frame as csv
            other.private$species_gc.avg %>%
            dplyr::select(-c(GC.skew, AT.skew)) %>% 
            write.csv(
                file = paste0("../data/Raw/Michael_Lynch", 
                ifelse(private$trek_scale, "/Trek_scale_", "/"), 
                "GC_vs_Rates.csv"),
                row.names = FALSE
            )

            other.private$species_gc.avg %>%
            dplyr::select(-GC.average) %>% 
            write.csv(
                file = paste0("../data/Raw/Michael_Lynch", 
                ifelse(private$trek_scale, "/Trek_scale_", "/"), 
                "GC_AT_skew_vs_Rates.csv"),
                row.names = FALSE
            )
            
            # save new data frame as csv
            if(private$trek_scale){
                private$lynch_rates %>%
                    select(-c(2:3)) %>%
                    write.csv(
                        file = paste0("../data/Raw/Michael_Lynch/Trek_scale_Lynch",
                                      "-2010-mutation-rates_FREQUENCY.csv"), 
                        row.names = FALSE
                    )
            }

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n") 
        },

        #' @description
        #' Plotting the mutation rate constants from the species in the PNAS paper
        #' along the theoretical GC curve.
        #' @return None.
        get_gc_content_curves = function(){
            cur.msg <- "Plotting species mutation rates along theoretical GC curve"
            l <- paste0(rep(".", 70-nchar(cur.msg)), collapse = "")
            cat(paste0(cur.msg, l))

            # helper functions
            theoretical <- function(x){
                return(x/(1-x))
            }

            hist.count <- function(data){
                p <- hist(data[, "GC"], plot = FALSE)
                return(p$counts)
            }

            euk <- c(
                "H.sapiens", "D.melanogaster", "C.elegans",	"A.thaliana",	
                "M.m.domesticus", "D.pulex", "P.pacificus",	
                "D.magna", "P.troglodytes", "A.nancymaae"
            )

            prok <- c(
                "E.coli", "P.luminescens ATCC29999", "T.turnerae", "S.cerevisiae", 
                "M.smegmatis", "P.fluorescens ATCC948", "R. toruloides"
            )

            private$species_gc <- private$species_gc %>% 
                dplyr::mutate(
                    kingdom = ifelse(Species %in% euk, "eukaryotes", "prokaryotes"),
                    kingdom.col = ifelse(kingdom == "eukaryotes", "purple", "darkgreen")
                ) %>% 
                dplyr::arrange(dplyr::desc(kingdom))
            prok.hist <- data.frame(GC = private$species_tolerance$prokaryotes$G_plus_C)
            euk.hist <- data.frame(GC = private$species_tolerance$eukaryotes$G_plus_C)
            vir.hist <- data.frame(GC = private$species_tolerance$viruses$G_plus_C)

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
            df <- rbind(theo.gc.line, private$species_gc)

            # # euclidean distance from point to tangent
            # spl <- smooth.spline(theo.gc.line$GC.average, theo.gc.line$Rates, spar=0.3)
            # newx <- seq(min(theo.gc.line$GC.average), max(theo.gc.line$Rates), 0.1)
            # pred <- predict(spl, x=newx, deriv=0)

            # eucl.dist <- sapply(1:nrow(private$species_gc), function(x){
            #   newx  <- private$species_gc[x, "GC.average"]
            #   pred0 <- predict(spl, x=newx, deriv=0)
            #   pred1 <- predict(spl, x=newx, deriv=1)
            #   yint  <- pred0$y - (pred1$y*newx)
            #   xint  <- -yint/pred1$y

            #   tangent <- data.frame(
            #     x = df$GC.average,
            #     y = pred1$y*df$GC.average + yint) %>% 
            #     dplyr::arrange(x)

            #   distance <- tangent[match(newx, tangent$x), "y"]-private$species_gc[x, "Rates"]

            #   return(distance)
            # })
            # eucl.df <- cbind(private$species_gc, eucl.dist) 
            # p.den <- eucl.df %>% 
            #   ggplot(aes(x = eucl.dist, fill = kingdom)) + 
            #   geom_density(alpha = 0.5) +  
            #   scale_fill_manual(
            #     values = unique(eucl.df$kingdom.col),
            #     limits = unique(eucl.df$kingdom)
            #   ) + 
            #   theme_bw() + 
            #   theme(
            #     legend.position = "none",
            #     panel.grid.major = element_blank(),
            #     panel.grid.minor = element_blank()
            #   ) + 
            #   labs(
            #     x = "Euclidean Distance",
            #     y = "Density"
            #   )

            # ggsave(
            #   filename = paste0("../figures/Chargaff_Equilibrium/
            #                     "Theoretical_GC/Euclidean_Dist.", private$save_as),
            #   plot = p.den,
            #   height = 7, width = 10
            # )

            # theoretical vs. actual GC content for each species
            eucl.dist <- sapply(1:nrow(private$species_gc), function(x){
                actual.rates  <- private$species_gc[x, "Rates"]
                nearest.ind <- which.min(abs(theo.gc.line$Rates-actual.rates))
                abs.diff <- abs(
                    theo.gc.line$GC.average[nearest.ind]-private$species_gc[x, "GC.average"]
                )
                return(abs.diff)
            })
            eucl.df <- cbind(private$species_gc, eucl.dist)

            as_tibble(eucl.df) %>% 
                dplyr::arrange(eucl.dist) %>% 
                dplyr::group_by(kingdom) %>% 
                dplyr::summarise(avg.dist = mean(eucl.dist))

            p.den <- as_tibble(eucl.df) %>% 
                dplyr::arrange(eucl.dist) %>% 
                dplyr::mutate(Species = forcats::fct_inorder(Species)) %>% 
                ggplot(aes(x = eucl.dist, y = Species, fill = kingdom)) + 
                geom_bar(stat = "identity") + 
                scale_fill_manual(
                    values = unique(eucl.df$kingdom.col),
                    limits = unique(eucl.df$kingdom)
                ) + 
                theme_bw() + 
                theme(
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank()
                ) + 
                labs(
                    x = "GC content difference (theoretical vs. actual)",
                    y = "Species"
                )
            
            dir.create(
                path = "../figures/Chargaff_Equilibrium/Theoretical_GC/",
                showWarnings = FALSE,
                recursive = TRUE
            )
            ggsave(
                filename = paste0("../figures/Chargaff_Equilibrium/", 
                                  "Theoretical_GC/GCcontent_Dist.", 
                                  private$save_as),
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
                size = 3) + 
                coord_cartesian(
                    xlim = c(0, 100), 
                    ylim = c(-1.5, 12)
                )

            ggsave(
                filename = paste0("../figures/Chargaff_Equilibrium/", 
                                  "Theoretical_GC/GCcontenthist_other", 
                                  "species_withlabels.", private$save_as),
                plot = p1,
                height = 8, width = 10
            )

            p2 <- p + 
                coord_cartesian(
                    xlim = c(0, 100), 
                    ylim = c(0, 12)
                )

            ggsave(
                filename = paste0("../figures/Chargaff_Equilibrium/", 
                                  "Theoretical_GC/GCcontenthist_other", 
                                  "species.", private$save_as),
                plot = p2,
                height = 7, width = 10
            )

            total.time <- Sys.time() - t1
            cat("DONE! --", signif(total.time[[1]], 2), 
                attr(total.time, "units"), "\n") 
        }
    )
)