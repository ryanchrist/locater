#' Fit Null Ordinary Least Squares Model
#' @param y a n x m matrix of m quantitative phenotypes (one phenotype / column)
#' @param A a n x q matrix of q background covariates
#'@export
FitNull <- function(y, A = NULL){

  if(!is.matrix(y)){y <- as.matrix(y)}
  if(is.null(A)){A = matrix(1,nrow(y),1)}

  if(nrow(y) != nrow(A)){
    stop("nrow(y) must equal nrow(A)")
  }

  if(!any(apply(A,2,function(x){all(x == x[1])},simplify = T))){
    A <- cbind(1,A) # add intercept to A if no intercept in input
  }

  Q <- qr.Q(qr(A))
  y.resids <- y - Q %*% crossprod(Q,y)
  sumsq <- colSums(y.resids^2)

  list("y" = y,
       "Q" = Q,
       "sumsq" = sumsq)
}

#' Perform Single Marker Testing (SMT) on a particular variant
#' @param x a null model object as returned by \link{FitNull}
#' @param g a genotype vector
#' @param add.noise if 'raw' return raw residuals, if FALSE simply rank normalize residuals, if TRUE add structured noise along the subspace spanned by the background covs and genotype so that the residuals are 'full rank': truly from N(0,I), not a curved Gaussian.
#' @return a list with three elements
#' \itemize{
#' \item `p-value`: single marker test p-value
#' \item `y`: the residual phenotype vector after testing the genotype
#' \item `Q`: an orthogonal matrix with columns spanning the column space of the tested genotype and the background covariates in the null model, useful for further testing
#' }
#'@export
TestMarker <- function(x, g, add.noise = FALSE){
  if(all(g==g[1])){
    return(list("p.value" = rep(NA_real_,ncol(x$y)),
                "y" = x$y, # iid Gaussians under null
                "Q" = x$Q))}

  if(is.character(add.noise)){
    add.noise <- switch(add.noise,
                        "TRUE" = TRUE,
                        "FALSE" = FALSE,
                        "raw" = "raw",
                        stop("invalid add.noise parameter.  Must be TRUE, FALSE or 'raw' "))
  }

  n <- nrow(x$Q)
  g <- g - x$Q %*% crossprod(x$Q,g)
  g <- g / sqrt(sum(g^2))
  Q.local <- cbind(g,x$Q)
  temp <- crossprod(Q.local,x$y)
  Z2 <- c(temp[1,])^2


  y.resids.local <- x$y - Q.local %*% temp # just return this if add.noise == "raw"

  if(is.logical(add.noise)){
    if(add.noise){
      y.resids.local <- scale(y.resids.local, center = FALSE) + c(Q.local %*% crossprod(Q.local,rnorm(n)))
    } else {
      y.resids.local <- apply(y.resids.local,2,function(x){rank2gauss(rank(x,ties.method = "random"))})
    }
  }

  # use statmod::qresiduals (randomized quantile residuals to extend this to GLMs)

  nu <- n-ncol(x$Q)-1
  list("p.value" = exp(pf(q = Z2*nu / (x$sumsq - Z2) , df1 = 1, df2 = nu, lower.tail = FALSE, log.p = TRUE)),
       "y" = y.resids.local, # iid Gaussians under null
       "Q" = Q.local)
}

#' Test All Cached Markers
#' Fast method for performing single marker testing across all variants stored in the \code{kalis::CacheHaplotypes}
#' @param y a \code{n} x \code{m} matrix of \code{m} quantitative phenotypes (one phenotype / column)
#' @param A a \code{n} x \code{q} matrix of \code{q} background covariates
#' @param from a positive integer giving the index of the first variant in the cache that should be tested, default = 1L.
#' @param to a positive integer giving the index of the last variant in the cache that should be tested, default = \link{\code{kalis::L()}}, the last variant.
#' @param ploidy a positive integer giving the ploidy of the organisms whose haplotypes were cached, default = 2L.
#' @param model a character in "additive", "recessive", or "dominant" giving the model of inheritance that should be tested, default = "additive"
#' @return a data.frame of -log10 p-values. Each row corresponds to a variant tested (in order starting with the variant with index \code{from}) and each of the \code{m} columns corresponds to a phenotype provided in \code{y} (in the same order as the columns of \code{y}).
#' @export
TestCachedMarkers <- function(y, A = NULL, from = 1L, to = kalis::L(), ploidy = 2L, model = "additive"){

  if(is.null(kalis::N())){stop("haplotypes must be cached with kalis::CacheHaplotypes before running FindTargetVars")}
  if(length(from)!=1 || from > kalis::L() || from < 1){
    stop("from must be an integer in [1,kalis::L()]")}
  if(length(to)!=1 || to > kalis::L() || to < 1){
    stop("to must be an integer in [1,kalis:L()]")}
  if(from>to){
    stop("from must be <= to")}

  if(!is.matrix(y)){y <- as.matrix(y)}
  if(is.null(A)){A = matrix(1,nrow(y),1)}

  if(nrow(y) != nrow(A)){
    stop("nrow(y) must equal nrow(A)")
  }

  X <- Haps2Genotypes(QueryCache(loci.idx = from:to),
                      ploidy = ploidy,
                      method = model)
  res <- fast_smt(y,X,qr.Q(qr(A)))
  if(!is.null(colnames(y))){colnames(res) <- colnames(y)}
  rownames(res) <- from:to

  as.data.frame(res)
}


# Test showing that TestMarker, TestCachedMarkers and R's native lm all give same result
# require(locater)
# require(kalis)
# target.l <- 10
# CacheHaplotypes(SmallHaps)
#
# p <- 1000
# A <- matrix(1,nrow=N()/2,1)
# y <- matrix(rnorm(p*N()/2),ncol = p)
# temp <- TestCachedMarkers(y,A)
#
# g <- t(locater:::Haps2Genotypes(QueryCache(target.l), ploidy = 2L))
# x <- locater:::FitNull(y,A)
#
# # TestMarker and TestCachedMarkers disagree.
#
# temp.res <- data.frame("lm" = numeric(p),
#                        "tm" = -log10(locater:::TestMarker(x,g)$p.value),
#                        "tcm" = as.numeric(unlist(temp[target.l,])))
# for(i in 1:p){
#   ms <- summary(lm(I(y[,i]) ~ g))
#   temp.res$lm[i] <- -pf(ms$fstatistic[1],ms$fstatistic[2],ms$fstatistic[3],lower.tail = FALSE,log.p = TRUE)/log(10)
#   print(i)
# }
#
# plot(temp.res$lm,temp.res$tm); abline(0,1)
# plot(temp.res$lm,temp.res$tcm); abline(0,1)

