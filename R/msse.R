
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
  cdf.hat <- QForm:::wrap.QFcdf(MSSETestData[["exclude1"]])
  k <- MSSETestData[["exclude1"]]$k
  a.r <- MSSETestData[["exclude1"]]$a.r
  b.r <- MSSETestData[["exclude1"]]$b.r
  internal_mpset(cdf.hat,k,a.r,b.r,x,lower.tail,log.p)
}


#' MSSE Test
#' Test MSSE
#' @param x vector of p-values
#' @param y vector of p-values
#' @param z vector of p-values
#' @export
msse.test <- function(x,y,z, test.1.solo = TRUE){

  if(!is.vector(x) || !is.vector(y) || !is.vector(z) || length(x)!=length(y) || length(x)!=length(z)){
    stop("x,y,z must all be vectors of equal length")}

  if(any(x>1 | x<0 | y>1 | y<0 | z>1 | z<0)){
    stop("all entries of x, y, z must be in [0,1]")}

  x <- -log(x)
  y <- -log(y)
  z <- -log(z)

  tm <- pmax(-pgamma(x+y,shape = 2,lower.tail = FALSE,log.p = TRUE),
                -pgamma(x+z,shape = 2,lower.tail = FALSE,log.p = TRUE),
                -pgamma(x+y+z,shape = 3,lower.tail = FALSE,log.p = TRUE))

  if(test.1.solo){
    tm <- pmax(x,tm)
    mpset_cdf_include1(tm, lower.tail = FALSE)
  } else {
    mpset_cdf_exclude1(tm, lower.tail = FALSE)
  }
}
# ans <- replicate(5e3,mpse.test(rexp(32)))
# shapiro.test(qnorm(ans))
#
# ans <- rep(0,1e4)
# for(i in 1:length(ans)){
#   ans[i] <- mpse.test(rexp(64))
#   print(i)
# }
# F.hat <- ecdf(-log10(ans))
# tail.est <- function(x){-log10(1-F.hat(x))}
# xx <- seq(0,5,len=1e3)
# plot(xx,tail.est(xx),type="l")
# abline(0,1,lty=2)


