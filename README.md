# Genomic PR-2 Simulations

## R&D pipeline for: Generalised interrelations among mutation rates drive the genomic compliance of Chargaff’s second parity rule.

Chargaff’s second parity rule (PR-2), where the complementary base (and k-mer) contents are matching within the same strand of a double stranded DNA (dsDNA), is a phenomenon that invited many elaborations for its explanation. The strict compliance of nearly all nuclear dsDNA genomes to PR-2 implies that the explanation should also be similarly strict and constraining. In this work, we revisited the possibility of mutation rates being the drivers behind PR-2. Starting from the assumption-free approach, we constructed a set of kinetic equations for unconstrained simulations. The results were analysed in terms of their PR-2 compliance by employing symbolic regression and machine learning techniques.

## Installation

Download the source code from the [GitHub](https://github.com/SahakyanLab/GenomicPR2Simulations) repository. You can also do that *via* a Linux/Unix/OSX command line, given that git is installed, by typing the following:

```
git clone https://github.com/SahakyanLab/GenomicPR2Simulations
```

## Requirements

* [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) >= 1.3.2 (for dplyr, tidyr, ggplot2, stringr)
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) >= 2.3
* [deSolve](https://cran.r-project.org/web/packages/deSolve/index.html) >= 1.30
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html) >= 1.5.1
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)  >= 1.0.16
* [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) >= 1.8.2
* [truncnorm](https://cran.r-project.org/web/packages/truncnorm/index.html) >= 1.0-8
* [caret](https://cran.r-project.org/web/packages/caret/index.html) >= 6.0-90
* [MLeval](https://cran.r-project.org/web/packages/MLeval/index.html) >= 0.3
* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html) >= 2.58.0 (required for processing raw data, but processed versions already deposited here.)

## Workflow to get similar plots in the publication
Execute the `submit.sh` file for sequential execution of the complete workflow.

### [01-genome_composition](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/01-genome_composition)

Downloads, filters and plots various meta-data of the species for the prokaryotes, eukaryotes and DNA virus kingdoms. Please note the following. First, the DNA Virus sequences are automatically downloaded from the latest updated list on NCBI, whereas the publication worked with a list of viruses released until Dec 2020. Second, if you execute the main `submit.sh` file, it will only plot the already processed data sets. If you wish to perform the full pipeline for the 3 kingdoms, please uncomment the lines within the `script.sh` file in this folder.

### [02-simulations](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/02-simulations)
Numerical simulations to produce the time evolution of genomic base content within 4.28-byr period. Please note, that each simulation run solves an ODE-based system and calculates the associated meta-data. In the main paper, the core simulations were run 25 million times each and took approximately 3-4 weeks each. If you wish to entirely replicate this, please amend the `Process.R` script in this folder so as to run this in parallel. By default, each simulation is run sequentially. Alternatively, for a simple demonstration purpose, please amend the `run` argument within the `script.sh` file in this folder from `25,000,000` to a lower number like `100,000`. Please also note that if the default of `25,000,000` is not used, all the subsequent results below will be different from the main paper.

### [03-machine_learning](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/03-machine_learning)
Performs the machine learning strategy for classification of compliance and non-compliance with the second parity rule solutions. 

### [04-symbolic_regression](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/04-symbolic_regression)

Symbolic regression of mutation rate constants under the tolerance region of the second parity rule compliance. Equations were generated with the Eureqa modelling engine.

## Interactive web application

The corresponding interactive web application for demonstration purposes is accessible on [GenomicPR2SimsWebApp GitHub](https://github.com/SahakyanLab/GenomicPR2SimsWebApp) repository.
