# Genomic PR-2 Simulations

## R&D pipeline for: Generalised interrelations among mutation rates drive the genomic compliance of Chargaff’s second parity rule.

Chargaff’s second parity rule (PR-2), where the complementary base (and k-mer) contents are matching within the same strand of a double stranded DNA (dsDNA), is a phenomenon that invited many elaborations for its explanation. The strict compliance of nearly all nuclear dsDNA genomes to PR-2 implies that the explanation should also be similarly strict and constraining. In this work, we revisited the possibility of mutation rates being the drivers behind PR-2. Starting from the assumption-free approach, we constructed a set of kinetic equations for unconstrained simulations. The results were analysed in terms of their PR-2 compliance by employing symbolic regression and machine learning techniques.

## Installation

Download the source code from the [GitHub](https://github.com/SahakyanLab/GenomicPR2Simulations) repository. You can also do that *via* a Linux/Unix/OSX command line, given that git is installed, by typing the following:

```
git clone https://github.com/SahakyanLab/GenomicPR2Simulations
```

## Requirements

* [deSolve](https://cran.r-project.org/web/packages/deSolve/index.html) >= 1.30
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) >= 1.0.7
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) >= 1.3.2
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) >= 3.3.5
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) >= 2.3
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html) >= 1.5.1
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)  >= 1.0.16
* [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) >= 1.8.2
* [truncnorm](https://cran.r-project.org/web/packages/truncnorm/index.html) >= 1.0-8
* [caret](https://cran.r-project.org/web/packages/caret/index.html) >= 6.0-90
* [MLeval](https://cran.r-project.org/web/packages/MLeval/index.html) >= 0.3

## Workflow to get similar plots in the publication
Execute the `submit.sh` file for sequential execution of the complete workflow.

### [01-genome_composition](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/01-genome_composition)

Downloads, filters and plots various meta-data of the species for the prokaryotes, eukaryotes and DNA virus kingdoms. Please note that the DNA Virus sequences are automatically downloaded from the latest updated list on NCBI, whereas the publication worked with a list of viruses released until Dec 2020.

### [02-simulations](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/02-simulations)
Numerical simulations to produce the time evolution of genomic base content within 4.28-byr period. 

### [03-machine_learning](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/03-machine_learning)
Performs the machine learning strategy for classification of compliance and non-compliance with the second parity rule solutions. 

### [04-symbolic_regression](https://github.com/SahakyanLab/GenomicPR2Simulations/tree/master/04-symbolic_regression)

Symbolic regression of mutation rate constants under the tolerance region of the second parity rule compliance. Equations were generated with the Eureqa modelling engine.

## Interactive web application

The corresponding interactive web application for demonstration purposes is accessible on [GenomicPR2SimsWebApp GitHub](https://github.com/SahakyanLab/GenomicPR2SimsWebApp) repository.