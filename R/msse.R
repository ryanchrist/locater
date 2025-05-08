
internal_mpset <- function(cdf.hat,k,a.r,b.r,x,lower.tail,log.p){
  # Confirm algebraic solution to root matches what we get from machine:
  intersection.x <- (a.r + log(k) ) / (1-b.r)
  #uniroot(function(x){max(x-log(k),0) + cdf.hat(x,density = FALSE,lower.tail = FALSE, log.p = TRUE)},interval = c(5,50))

  y <- cdf.hat(x, density = FALSE, lower.tail = lower.tail, log.p = log.p)

  if(!is.finite(intersection.x) || intersection.x <= log(k)){ # the extrapolated CDF never crosses over the Bonferroni line
    # so we can just return the estimate
    return(y)
  }

  if(length(x) <=0){
    warning("x provided has length 0, returning NA")
    res <- NA_real_
  }

  if(length(x) == 1){
    res <- if(!log.p){
      if(lower.tail){
        if(x > intersection.x){1-k*exp(-x)}else{y}
      } else {
        if(x > intersection.x){k*exp(-x)}else{y}
      }
    } else {
      if(lower.tail){
        if(x > intersection.x){log1p(-k*exp(-x))}else{y}
      } else {
        if(x > intersection.x){-x+log(k)}else{y}
      }
    }
  }

  if(length(x) > 1){
    res <- if(!log.p){
      if(lower.tail){
        ifelse(x > intersection.x,1-k*exp(-x),y)
      } else {
        ifelse(x > intersection.x,k*exp(-x),y)
      }
    } else {
      if(lower.tail){
        ifelse(x > intersection.x,log1p(-k*exp(-x)),y)
      } else {
        ifelse(x > intersection.x,-x+log(k),y)
      }
    }
  }

  res
}


mpset_cdf_include1 <- function(x, lower.tail = TRUE, log.p = FALSE){
  cdf.hat <- QForm:::wrap.QFcdf(MSSETestData[["include1"]])
  k <- MSSETestData[["include1"]]$k
  a.r <- MSSETestData[["include1"]]$a.r
  b.r <- MSSETestData[["include1"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}

mpset_cdf_exclude1 <- function(x, lower.tail = TRUE, log.p = FALSE){
  # cdf.hat <- QForm:::wrap.QFcdf(MSSETestData[["exclude1"]])
  # k <- MSSETestData[["exclude1"]]$k
  # a.r <- MSSETestData[["exclude1"]]$a.r
  # b.r <- MSSETestData[["exclude1"]]$b.r
  cdf.hat <- QForm:::wrap.QFcdf(MSSETestData[[2]])
  k <- MSSETestData[[2]]$k
  a.r <- MSSETestData[[2]]$a.r
  b.r <- MSSETestData[[2]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}


#' MSSE Test
#' A p-value combination test that combines three independent p-values (element wise if three vectors are given) while prioritizing the first \code{x}, as proposed in our pre-print \link{www.biorxiv.org/content/10.1101/2024.09.30.615852}. See Details
#'
#' Let \code{f} be a Fisher combination function that maps a set of -log10 p-values to a single -log10 p-value using the classic Fisher combination approach.
#' In a setting where we wish to combine three -log10 p-values -- \code{x}, \code{y}, and \code{z} -- using \code{f(x,y,z)} treats all three -log10 p-values exchangeably.
#' However, in some settings, we may be using \code{y} and \code{z} to boost the signal in \code{x}. In other words, we may be able to safely assume that \code{f(x,y,z)} has no chance of being significant is \code{x} is not at least somewhat significant.
#' We run into this setting with LOCATER where \code{x} corresponds to the single marker test (SMT) -log10 p-value while \code{y} and \code{z} correspond to tests performed via Stable Distillation and QForm respectively.
#' In this context, MSSE gains power over Fisher combination by ignoring cases where \code{y} and/or \code{z} would be significant without \code{x} being somewhat significant.
#'
#' @param x vector of -log10 p-values
#' @param y vector of -log10 p-values
#' @param z vector of -log10 p-values
#' @param test.1.solo logical, if FALSE ignore situations where \code{x} might be significant by itself and focus statistical power on cases where \code{x} is only significant in combination with \code{y} and/or \code{z}. By default = TRUE.
#' @return a vector of -log10 p-values
#' @export
msse.test <- function(x,y,z, test.1.solo = TRUE){

  if(!is.vector(x) || !is.vector(y) || !is.vector(z) || length(x)!=length(y) || length(x)!=length(z)){
    stop("x,y,z must all be vectors of equal length")}

  x <- x * log(10)
  y <- y * log(10)
  z <- z * log(10)

  to.exclude <- which(!is.finite(x) | !is.finite(y) | !is.finite(z))

  if(length(to.exclude) == length(x)){return(rep(NA_real_,length(x)))}

  if(length(to.exclude)){
    res <- rep(NA_real_,length(x))
    x <- x[-to.exclude]
    y <- y[-to.exclude]
    z <- z[-to.exclude]
  }

  tm <- pmax(-pgamma(x+y,shape = 2,lower.tail = FALSE,log.p = TRUE),
             -pgamma(x+z,shape = 2,lower.tail = FALSE,log.p = TRUE),
             -pgamma(x+y+z,shape = 3,lower.tail = FALSE,log.p = TRUE))

  if(test.1.solo){
    tm <- pmax(x,tm)
    if(length(to.exclude)){
      res[-to.exclude] <- -mpset_cdf_include1(tm, lower.tail = FALSE, log.p = TRUE)/log(10)
    } else {
      res <- -mpset_cdf_include1(tm, lower.tail = FALSE, log.p = TRUE)/log(10)
    }
  } else {
    if(length(to.exclude)){
      res[-to.exclude] <- -mpset_cdf_exclude1(tm, lower.tail = FALSE, log.p = TRUE)/log(10)
    } else {
      res <- -mpset_cdf_exclude1(tm, lower.tail = FALSE, log.p = TRUE)/log(10)
    }
  }
  res
}

# n <- 1e3
# x <- rexp(n,log(10))
# y <- rexp(n,log(10))
# z <- rexp(n,log(10))
# t1 <- locater::msse.test(x,y,z,test.1.solo = TRUE)
# shapiro.test(qnorm(10^-t1))
# options(error=recover)
# t2 <- locater::msse.test(x,y,z,test.1.solo = FALSE)
# shapiro.test(qnorm(10^-t2))
