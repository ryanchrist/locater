# Introduction

`locater` is an R package for testing tree structures (represented via a set of discrete clades and/or a relatedness matrix) for association with traits of interest. Please see our pre-print describing the underlying LOCATER methodology here: https://www.biorxiv.org/content/10.1101/2024.09.30.615852 . We recently applied `locater` to a richly phenotyped whole genome sequencing described in our pre-print here https://www.medrxiv.org/content/10.1101/2024.11.04.24316696.

While `locater` natively interacts with the local ancestry inference software package `kalis`, available here https://github.com/louisaslett/kalis/, the package exports functions exposing all internal testing routines so that LOCATER may be used in conjunction with other local genealogy inference engines/software. The exported function `TestLoci` is a wrapper function used for testing target variants along a chromosome; it provides a template for users hoping to adapt `locater` to use a different ancestry inference engine or otherwise customize the testing routine.

# Installation Instructions

## Docker Image
To try out locater immediately, check out our docker image available at https://hub.docker.com/repositories/rchrist7 .  On a system with Docker installed, R can be launched in a interactive session under this image by running:  
```{bash docker, eval=FALSE}
docker run -it rchrist7/mini-shark /bin/bash R
```
Call  `require(locater)` in the R session to load locater.  To get started, a simple vignette is available under the Articles tab at https://ryanchrist.github.io/locater/ alongside the documentation. 

## Install R 
If not already installed, please visit https://posit.co/download/rstudio-desktop/ for guidance. Once R is installed, we suggest installing them in the order below by copying and pasting the following R commands into R.

## Install `kalis`
Our default ancestry inference engine `kalis` may be installed from our [public Github repository](https://github.com/louisaslett/kalis/) as follows.

```
if(!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
# install rdhf5, this requires zlib (if not already installed, try sudo apt install zlib1g-dev on Linux; brew install zlib on MacOS)
BiocManager::install("rhdf5") # a suggested dependency for kalis that will be required for running locater

install.packages("remotes")
remotes::install_github("louisaslett/kalis",
  configure.vars = c(kalis = "PKG_CFLAGS='-march=native -mtune=native -O3'"))
if(!require(kalis)){stop("kalis was not properly installed")}
```

## Install Dependencies from public r-universe respository

```{r install_from_r_universe,eval=FALSE}
install.packages("RcppRoll")
install.packages("QForm", repos = "https://ryanchrist.r-universe.dev", dependencies=TRUE)
if(!require(QForm)){stop("Qform was not properly installed")}

install.packages("renyi", repos = "https://ryanchrist.r-universe.dev")
if(!require(renyi)){stop("renyi was not properly installed")}
```

## Install LOCATER from GitHub
```{r install_locater,eval=FALSE}
remotes::install_github("ryanchrist/locater")
if(!require(locater)){stop("locater was not properly installed")}
```

# Vignettes
Please see the introductory package vignette in `vignettes/simple_locater_example.Rmd` for a simple example deploying `locater` and a demonstration of how the `locater` package API can be used to test inferred clades and relatedness matrices generated by other ancestry inference methods. This vignette can be built and viewed in R by calling the following.

```{r build_vignette,eval=FALSE}
remotes::install_github("ryanchrist/locater", build_vignettes=TRUE)
vignette("simple_locater_example",package="locater")
```
