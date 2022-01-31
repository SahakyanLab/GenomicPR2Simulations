# ATGC Network Research

## R&D pipeline for the simulation-based explanation for the Chargaff's second parity rule phenomenon.

Simulation-based models are constructed based on the presence of natural equalities of mutation rates in the double stranded DNA. The models express the individual base contents through the underlying mutation rate constants.

## Installation

Download the source code from the [GitHub](https://github.com/SahakyanLab/ATGCNetworkResearch) repository. You can also do that *via* a Linux/Unix/OSX command line, given that git is installed, by typing the following:

```
git clone https://github.com/SahakyanLab/ATGCNetworkResearch
```

## Requirements

* [deSolve](https://cran.r-project.org/web/packages/deSolve/index.html) >= 1.30
* [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) >= 1.0.7
* [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html) >= 3.3.5
* [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html) >= 2.3
* [foreach](https://cran.r-project.org/web/packages/foreach/index.html) >= 1.5.1
* [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html)  >= 1.0.16
* [doRNG](https://cran.r-project.org/web/packages/doRNG/index.html) >= 1.8.2
* [truncnorm](https://cran.r-project.org/web/packages/truncnorm/index.html) >= 1.0-8
* [caret](https://cran.r-project.org/web/packages/caret/index.html) >= 6.0-90
* [MLeval](https://cran.r-project.org/web/packages/MLeval/index.html) >= 0.3

## Workflow to get similar plots in the publication

Execute the `submit.sh` file for sequential execution of the complete workflow. Change the parameters as follows:

* parallel = FALSE (default) or TRUE for parallel execution of simulations.
* pca = FALSE (default) or TRUE for PCA performance.
* ncpu = 4 (default). Input any positive integer for processes requiring multiple central processing units.

### [01-genome_composition](https://github.com/SahakyanLab/ATGCNetworkResearch/tree/master/01-genome_composition)

Downloads, filters and plots various meta-data of the species for the prokaryotes, eukaryotes and DNA virus kingdoms. Please note that the DNA Virus sequences are automatically downloaded from the latest updated list on NCBI, whereas the publication worked with a list of viruses released until Dec 2020.

### [02-simulations](https://github.com/SahakyanLab/ATGCNetworkResearch/tree/master/02-simulations)

Numerical simulations to produce the time evolution of genomic base content within 4.28-byr period. Performs the machine learning strategy for classification of compliance and non-compliance with the second parity rule solutions. Optional PCA performance.

### [03-symbolic_regression](https://github.com/SahakyanLab/ATGCNetworkResearch/tree/master/03-symbolic_regression)

Symbolic regression of mutation rate constants under the tolerance region of the second parity rule compliance. Equations were generated with the Eureqa modelling engine.

## Interactive web application

The interactive ATGC dynamics solver web application is accessible on [ATGCsolver GitHub](https://github.com/SahakyanLab/ATGCsolver) repository.