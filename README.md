# Introduction

`locater` is an R package for testing tree structures (represented via discrete and/or via relatedness matrices) for association with traits of interest. Please see our pre-print describing the underlying LOCATER methodology here: https://www.biorxiv.org/content/10.1101/2024.09.30.615852 .

While `locater` natively interacts with the ancestry inference software package `kalis`, available here https://github.com/louisaslett/kalis/, the package exposes all internal testing functions so that LOCATER may be used in conjunction with any ancestry inference software without requiring the user to re-implement the core testing routines. The exported function `TestLoci` is a wrapper function used for testing target variants along a chromosome; it provides a template for users hoping to adapt `locater` to use a different ancestry inference engine or otherwise customize the testing routine.



# Installation Instructions


## Install R 
If not already installed, please visit https://posit.co/download/rstudio-desktop/ for guidance. Once R is installed, we suggest installing them in the order below by copying and pasting the following R commands into R.

## Install `kalis`
Our default ancestry inference engine `kalis` may be installed from our [public Github repository](https://github.com/louisaslett/kalis/) as follows.

```
if(!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
BiocManager::install("rhdf5") # a suggested dependency for kalis that will be required for running LOCATER

install.packages("remotes")
remotes::install_github("louisaslett/kalis",
  configure.vars = c(kalis = "PKG_CFLAGS='-march=native -mtune=native -O3'"))
if(!require(kalis)){stop("kalis was not properly installed")}
```

## Install Dependencies and LOCATER from out public r-universe respository

```{r install_from_r_universe,eval=FALSE}
install.packages("QForm", repos = "https://ryanchrist.r-universe.dev")
if(!require(QForm)){stop("Qform was not properly installed")}

install.packages("renyi", repos = "https://ryanchrist.r-universe.dev")
if(!require(renyi)){stop("renyi was not properly installed")}

install.packages("locater", repos = "https://ryanchrist.r-universe.dev")
if(!require(locater)){stop("locater was not properly installed")}
```

## Post Installation
Please open `vignettes/simple_locater_example.html` for a vignette showing simple example of `locater` in action and a demonstration of how the `locater` package API can test inferred clades and relatedness matrices provided by any ancestry inference method.

## Replicating Paper Figures

To replicate our power results, we require a couple of dependencies to enable haplotype simulation.

### Installing dependencies for haplotype simulation

First, [install msprime](
https://tskit.dev/msprime/docs/stable/installation.html).

Then install our custom python package, champs, for making simulation calls to msprime by first installing the pandas package for python and then installing champs from the source tar.gz provided. From command line, this can be easily done with pip as follows.

```{bash install, eval=FALSE}
pip3 install pandas
pip3 install champs.tar.gz
```

### Replicating Paper Power Figures

All scripts used to run our power simulations are then available within the `shark` directory. All simulations were run using the IBM LSF job management system on an academic university cluster using a job submission template found in `shark/simulation/association_sims/ryan_run_locater_sim.sh`. That script can be used to call to any of the four simulation R scripts in `shark/simulation/ryan_sims_parameters` described below.

- locater12.R was used to run our null simulations
- locater14.R was used to run our power simulations assuming all causal variants were observed and fell within a 10kb window
- locater15.R was used to run our power simulations assuming all causal variants were hidden and fell within a 10kb window
- locater16.R was used to run our power simulations assuming all causal variants were observed and fell within a 100kb window
- locater17.R was used to run our power simulations assuming all causal variants were hidden and fell within a 100kb window

The results of these five scripts were then written to intermediary files (.rds files) which were then read into R and used to produce our paper figures using the corresponding visualization scripts found in `shark/simulation/association_sims/ryan_viz/`.
