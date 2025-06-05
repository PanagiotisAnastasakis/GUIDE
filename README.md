# Decomposition of Genome-Phenome Associations

This repository contains the analysis of methods for decomposing genetic associations between variants and phenotypes using GWAS summary statistics data. The following methods are included:

1. [Genetic Umixing by Independent Decomposition (GUIDE)](https://www.biorxiv.org/content/10.1101/2024.05.03.592285v2)

2. [Decomposition of Genetic Associations (DeGAs)](https://www.nature.com/articles/s41467-019-11953-9)

3. [Bayesian Non-Negative Matrix Factorization (bNMF)](https://pubmed.ncbi.nlm.nih.gov/30240442/)

The first two methods are extensively compared, with a particular focus a GUIDE. The implementation of GUIDE has been extended through the use of ICASSO framework. 


Details on this extension can be found in the file `Thesis.pdf`. 


### Files

All analyses in this repository, located in the `scripts` folder, were used in my MSc thesis, which can be found in the file `Thesis.pdf`.

* `GUIDE.R`: Contains functions for running GUIDE, extended with ICASSO.
* `DeGAs.R`: Contains the function `get_tsvd`, which is used to run DeGAs.
* `bNMF.R`: Contains functions for running bNMF, originating from [this](https://github.com/gwas-partitioning/bnmf-clustering) repository. The implementation of bNMF is not part of this work, but may be included in future analyses.
* `Simulation_functions.R`: Contains functions related to the simulation of genetic architectures, used to evaluate GUIDE and compare it with DeGAs across diverse scenarios.
* `Additional_helper_functions.R`: Here, additional functions used throughout are included.
* `Example_usage.R`: This file has a simple example demonstrating how one can use this implementation of GUIDE and DeGAs.
* `Block_simulation.Rmd`, `Polygenic_simulation_1`, `Polygenic_simulation_2`, `Polygenic_simulation_3`: The 4 R Markdown files contain all the simulation analyses, from generating the data, storing it and visualizing it.
* `T2D_data_analysis`: The script of this file was used to analyze T2D data from the paper ['Multi-ancestry polygenic mechanisms of
type 2 diabetes'](https://www.nature.com/articles/s41591-024-02865-3) by K. Smith et al. (2024).


