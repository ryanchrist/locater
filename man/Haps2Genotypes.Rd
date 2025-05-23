% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{Haps2Genotypes}
\alias{Haps2Genotypes}
\title{Collapse a matrix of haplotypes to genotypes}
\usage{
Haps2Genotypes(haps, ploidy = 2L, method = "additive")
}
\arguments{
\item{haps}{a vector or \code{p} x \code{N} matrix of 1s and 0s for \code{N} haplotypes observed at \code{p} loci}

\item{ploidy}{the ploidy of the underlying organisms, haplotypes corresponding to the same organism are assumed to be stored consecutively in \code{haps}. Default = 2 (diploid)}

\item{method}{a character in "additive", "recessive", or "dominant" giving the model of inheritance that should be tested, default = "additive"}
}
\value{
a matrix of genotypes with dimensions \code{p} x \code{n} = \code{N} / \code{ploidy}
}
\description{
Collapse a matrix of haplotypes to genotypes
}
