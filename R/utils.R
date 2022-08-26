Haps2Genotypes <- function(haps, # vector or p x N matrix of 1s and 0s
                           ploidy = 2L,
                           method = "additive"){
  ploidy.in <- ploidy
  ploidy <- as.integer(ploidy)
  if(ploidy!=ploidy.in || ploidy <= 0){stop("ploidy must be a positive integer.")}

  if(is.vector(haps)){haps <- matrix(haps,1)}

  N <- ncol(haps)

  if(ploidy > 1){
    genotypes <- haps[,seq.int(1,N,by = ploidy)]
    for(j in 2:ploidy){
      if(method == "additive"){
        genotypes <- genotypes + haps[,seq.int(j,N,by = ploidy)]
      } else if(method == "dominant"){
        genotypes <- genotypes | haps[,seq.int(j,N,by = ploidy)]
      } else {
        stop("method mis-specified: must be additive or dominant")
      }
    }
    storage.mode(genotypes) <- "integer"
  } else {
    genotypes <- haps
  }

  if(! "matrix" %in% class(genotypes)){
    genotypes <- matrix(genotypes,ncol=N/ploidy)
  }

  genotypes # a p x n matrix
}

fish <- function(..., na.rm = FALSE){
  # input vectors of -log10 p-values and output single vector of -log10 p-values
  # if all entries in a row of tests are NA, 0 will be returned for those rows
  tests <- do.call(cbind,list(...))
  if(na.rm){
    na.tests <- is.na(tests)
    tests[na.tests] <- 0
  }
  -pgamma(log(10) * rowSums(tests), shape = if(na.rm){ncol(tests) - c(rowSums(na.tests))}else{ncol(tests)},
          lower.tail = FALSE, log.p = T)/log(10)
}

fishv <- function(v){
  # single vector of negative log10 p-values (automatically ignores any non-finite entries)
  to.combine <- is.finite(v)
  -pgamma(log(10) * sum(v[to.combine]), shape = sum(to.combine),lower.tail = FALSE, log.p = T)/log(10)
}


fast_smt <- function(y,x,Q){
  # LRT for \beta_i = 0 for many OLS
  # models y_j ~ A\alpha + x_i \beta_i + \epsilon
  # testing many outcomes y_j against many x_i separately
  # y is a n by m matrix for n samples with m outcomes
  # x is a p by n matrix for n samples with p predictors
  # Q is a n by q matrix for n samples with q background predictors

  # returns p by m matrix with the -log10 p-value for every pair of predictors and outcomes

  if(nrow(y) != nrow(Q) | nrow(y) != ncol(x) ){
    stop("nrow(y) must equal nrow(Q) and ncol(x)")
  }

  y.resids <- y - Q %*% crossprod(Q,y)
  sumsq0 <- colSums(y.resids^2)
  nu <- nrow(Q)-ncol(Q)-1

  x <- x - (x %*% Q) %*% t(Q)

  Z2 <- (x%*%y.resids)^2 / rowSums(x^2)

  -pf(q = - Z2*nu / sweep(Z2,2,sumsq0,"-"), df1 = 1, df2 = nu, lower.tail = FALSE, log.p = TRUE)/log(10)
}
#
# # Example:
# n <- 10000
# p <- 2000
# #p <- 1
# m <- 200
# x <- matrix(sample(c(0L,1L),n * p,TRUE),p,2*n)
# x <- x[,seq(1,2*n,2)] + x[,seq(2,2*n,2)]
# y <- matrix(rnorm(n*m),n,m)
# A <- cbind(1,matrix(rnorm(n * 5),n,5))
# Q <- qr.Q(qr(A))
# system.time(neglog10pval <- fast_smt(y,x,Q))
# hist(10^-neglog10pval) # uniform p-values under the null

rank2gauss <- function(x){
  # x is a vector of ranks with no ties
  # See ?quantile Type 9 to help justify below transformation
  qnorm((seq_len(length(x)) - 3/8)/(length(x) + 1/4))[x]
}






# Below is general version for univariate and multivariate testing
CalcBounds <- function(f = function(k,args){eigs_sym(args[[1]],k)},
                       args,
                       obs,
                       traces,
                       k=c(20,200),
                       neg.log10.cutoff=6,
                       tau = 1, # may be a vector as long as obs
                       delta = 0, # may be a vector as long as obs
                       only.point.est = FALSE, # if TRUE, then the first element of k is used for calculated trunc part
                       only.bounds.est = FALSE,
                       eval.for.popstruc = FALSE,
                       unfinished = NULL,
                       lower.tail = TRUE,
                       other.test.res = NULL, # -log10 pvalues of other tests w/ list element per observation
                       parallel.sapply = base::sapply){

  res <- matrix(NA_real_,nrow=5,ncol=length(obs))
  res[1,] <- 1

  one.inflation.setting <- length(unique(tau)) == 1 & length(unique(delta)) == 1

  if(only.point.est){
    e <- f(k[1], args)

    if((length(traces$diag)/length(e$values))*sum(e$values^2) < traces$hsnorm2){
      # the leading e$values^2 can't sum up to the estimated hsnorm2 -- eigendecomposition is too
      # unstable and so the QForm test must be dropped for all phenotypes: we exit with all NAs
      # for the bounds and the point estimate
      return(rbind(1,k[1],matrix(NA_real_,nrow=3,ncol=length(obs))))
    }


    res[2,] <- k[1]

    if(one.inflation.setting){  # calculate function only once
      calc.func <- CalcBounds2(traces, e, tau[1], delta[1], only.point.est, only.bounds.est, parallel.sapply = parallel.sapply)
    }

    for(p in 1:length(obs)){
      if(!one.inflation.setting){  # calculate function separately for each phenotype
        calc.func <- CalcBounds2(traces, e, tau[p], delta[p], only.point.est, only.bounds.est, parallel.sapply = parallel.sapply)
      }
      res[5,p] <- calc.func(obs[p],lower.tail)
    }

    return(res)
  }


  if(eval.for.popstruc){ k <- max(k)} # eval.for.popstruc means just go directly to evaluating the largest k

  if(is.null(unfinished)){unfinished <- rep(TRUE,length(obs))}
  j <- 0

  while(any(unfinished)){

    j <- j+1 # advance to next k
    e <- f(k[j], args)

    if((length(traces$diag)/length(e$values))*sum(e$values^2) < traces$hsnorm2){
      # the leading e$values^2 can't sum up to the estimated hsnorm2 -- eigendecomposition is too
      # unstable and so the QForm test must be dropped for all phenotypes: we exit with all NAs
      # for the bounds and the point estimate
      return(rbind(1,k[j],matrix(NA_real_,nrow=3,ncol=length(obs))))
    }


    if(one.inflation.setting){ # calculate function only once
      if(exists("calc.func.backup")){rm(calc.func.backup)}
      calc.func <- CalcBounds2(traces, e, tau[1], delta[1], only.bounds.est = only.bounds.est, parallel.sapply = parallel.sapply)
    }

    for(p in which(unfinished)){
      res[2,p] <- k[j]

      if(!one.inflation.setting){ # calculate function separately for each phenotype
        if(exists("calc.func.backup")){rm(calc.func.backup)}
        calc.func <- CalcBounds2(traces, e, tau[p], delta[p], only.bounds.est = only.bounds.est, parallel.sapply = parallel.sapply)
      }

      res[3:5,p] <- calc.func(obs[p],lower.tail)

      # THIS SECTION ALLOWS POINT ESTIMATES TO BE USED WHENEVER THE BOUNDS / CALC.FUNC
      # FAILS -- THE CALC.FUNC.BACKUP FUNCTION BELOW
      if(
        if(only.bounds.est){any(is.na(res[3:4,p]))}else{any(is.na(res[3:5,p]))}
        # if only.bounds.est, then res[5,p] will always be NA!
      ){

        if(!exists("calc.func.backup")){
          calc.func.backup <- CalcBounds2(traces, e,
                                          tau[if(one.inflation.setting){1}else{p}],
                                          delta[if(one.inflation.setting){1}else{p}],
                                          only.point.est = TRUE, parallel.sapply = parallel.sapply)
        }
        res[3:5,p] <- calc.func.backup(obs[p],lower.tail)

      }


      if(!is.null(other.test.res)){
        bound.est <- fishv(c(other.test.res[[p]],res[ if(lower.tail){3}else{4}, p ]))
      } else {
        bound.est <- res[ if(lower.tail){3}else{4}, p ]
      }

      if(is.na(bound.est)){
        unfinished[p] <- FALSE # if we're still getting a NA for the bound, just exit here for now.
      } else if(bound.est  < neg.log10.cutoff ){
        res[1,p] <- 0
        unfinished[p] <- FALSE
      }
    }

    if(j==length(k)){unfinished <- FALSE} # we've run out of k to evaluate for now
  }

  # first row: still interesting indicator
  # second row: final k used to calculate this particular bound and point estimate
  # third row: -log10 lower bound
  # fourth row: -log10 upper bound
  # fifth row: -log10 pvalue point estimate

  res
}


CalcBounds2 <- function(traces, e = NULL, tau = 1, delta = 0, only.point.est = FALSE, only.bounds.est = FALSE, parallel.sapply = base::sapply){

  # this returns a function that is NOT vectorized (it can only take observed value at a time)
  # the returned function returns a vector of three numbers:
  # first entry: -log10 lower bound
  # second entry: -log10 upper bound (NULL if e is NULL)
  # third entry: -log10 pvalue point estimate (NULL if e is NULL)

  if(only.point.est & only.bounds.est){
    stop("only.point.est and only.bounds.est cannot both be TRUE")
  }

  if(tau < 0 | delta < 0 | length(tau) != 1 | length(delta) != 1){
    stop("tau and delta must be greater than or equal to 0 and have length 1.")
  }

  delta2 <- delta^2


  if(is.null(e)){ # Use fully Gaussian approximation
    mu <- tau * (1 + delta2) * traces$trace
    sigma <- tau * sqrt((2 + 4*delta2) * traces$hsnorm2)

    return(function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -pnorm(obs, mean = mu, sd = sigma,
             lower.tail = lower.tail, log.p = TRUE)/log(10)
    })
  }


  # IF THE AMOUNT OF VARIANCE LEFT IN THE REMAINDER IS ZERO OR NEGLIGIBLE IN RELATIVE OR ABSOLUTE TERMS,
  # JUST DROP THE REMAINDER AND USE A POINT ESTIMATE BASED ON THE TOP K EIGENVALUES
  if(length(e$values)==length(traces$diag) ||
     (1 - sum(e$values^2)/traces$hsnorm2) <= .Machine$double.eps^.25 ||
     traces$hsnorm2 - sum(e$values^2) <= sqrt(.Machine$double.eps)){
    # bounds are not needed because either we have all of the eigenvalues or the eigenvalues that
    # are in the remainder are non-zero due to numerical imprecision.
    gauss.tcdf <- QForm::QFGauss(tau * e$values, parallel.sapply = parallel.sapply)

    return(function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      rep(-gauss.tcdf(obs, lower.tail = lower.tail, log.p = TRUE)/log(10),if(only.point.est){1}else{3})
    })
  }


  # Calc Required Traces
  R.max.abs.eta <- tau * min(abs(e$values))
  R.sum.eta <-  tau * (traces$trace - sum(e$values))
  R.sum.etasq <- tau^2 * (traces$hsnorm2 - sum(e$values^2))
  if(R.sum.etasq < 0){stop("Squared HS norm of matrix is less than sum of eigenvalues squared -- check that the provided traces are from the same matrix as the eigenvalues")}
  R.sum.eta.deltasq <-  delta2 * R.sum.eta
  R.sum.etasq.deltasq <- delta2 * R.sum.etasq


  if(only.point.est){

    # Point Estimate Function
    gauss.tcdf <- suppressWarnings(QForm::QFGauss(tau * e$values, sigma = sqrt((2 + 4*delta2)*R.sum.etasq),parallel.sapply = parallel.sapply))
    E.R <- (1 + delta2) * R.sum.eta
    gauss.approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) gauss.tcdf(x-E.R, density = density, lower.tail = lower.tail, log.p = log.p)
    attr(gauss.approxfullcdf,"mu") <- attr(gauss.tcdf,"mu") + E.R
    attr(gauss.approxfullcdf,"Q.sd") <- attr(gauss.tcdf,"Q.sd")


    function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -gauss.approxfullcdf(obs, lower.tail = lower.tail, log.p=T)/log(10)
    }

  } else if(only.bounds.est){

    # Bound Function
    tcdf <- suppressWarnings(QForm::QFGauss(e$values,parallel.sapply = parallel.sapply))
    bound.func <- QForm::QFGaussBounds(tcdf,"identity", R.max.abs.eta, R.sum.eta, R.sum.etasq, R.sum.eta.deltasq, R.sum.etasq.deltasq)

    function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -c(log(bound.func(obs)[1,1:2 + if(lower.tail){0}else{2}]),NA_real_)/log(10)
    }

  } else {

    # Point Estimate Function
    gauss.tcdf <- suppressWarnings(QForm::QFGauss(tau * e$values, sigma = sqrt((2 + 4*delta2)*R.sum.etasq),parallel.sapply = parallel.sapply))
    E.R <- (1 + delta2) * R.sum.eta
    gauss.approxfullcdf <- function(x, density = FALSE, lower.tail = TRUE, log.p = FALSE) gauss.tcdf(x-E.R, density = density, lower.tail = lower.tail, log.p = log.p)
    attr(gauss.approxfullcdf,"mu") <- attr(gauss.tcdf,"mu") + E.R
    attr(gauss.approxfullcdf,"Q.sd") <- attr(gauss.tcdf,"Q.sd")

    # Bound Function
    tcdf <- suppressWarnings(QForm::QFGauss(e$values,parallel.sapply = parallel.sapply))
    bound.func <- QForm::QFGaussBounds(tcdf,"identity", R.max.abs.eta, R.sum.eta, R.sum.etasq, R.sum.eta.deltasq, R.sum.etasq.deltasq)

    function(obs, lower.tail = FALSE){
      if(length(obs) != 1){stop("while this function could be easily parallelized,
                              to keep the code simple, for now it only takes one obs at a time.")}
      -c(log(bound.func(obs)[1,1:2 + if(lower.tail){0}else{2}]),gauss.approxfullcdf(obs,lower.tail = lower.tail,log.p=T))/log(10)
    }
  }
}


calc_traces <- function(M, Q, sym.M = FALSE,
                        from_recipient = 1, nthreads = 1){
  if(sym.M){M <- 0.5*(M + t(M))}
  J <- crossprod(Q, M)
  tX <- t((Q %*% (J%*%Q)) - (M %*% Q))

  kalis:::CalcTraces(M,tX,t(Q),J,from_recipient,nthreads)
}



SimpleCalcBounds <- function(y,
                             matmul,
                             traces,
                             min.prop.var = 0.98,
                             k=c(0),
                             neg.log10.cutoff = NULL, #6
                             other.test.res = NULL, # -log10 pvalues of other tests w/ list element per observation
                             lower.tail = FALSE,
                             parallel.sapply = base::sapply){

  # this function uses the schedule k to evaluate more and more eigenvalues for the implicitly provided matrix until either
  # at least the top min.prop.var of the variance of the matrix is explained or the last k is reached.
  # If a neg.log10.cutoff is given, then bounds are calculated for each k and further eigendecomposition is stopped if all of the observed statistics have bounds that exclude them from being
  # more significant than the provided neg.log10.cutoff.

  # Motivation behind this approach: we are targeting the full quadratic form as our test statistic -- hence the bounds are exact.  But we typically stop early to get a truncated
  # estimate of the statistic and just report that p-value.  We ensure that our p-values are well calibrated by computing the truncated test statistic and comparing it to the truncated distribution.
  # The bounds here are simply a means of conserving compute

  # If any of the bounds are not finite due to some sort of numerical instability, then we simply fall back on evaluating min.prop.var of the variance and returning the point estimates

  m <- ncol(y)
  res <- data.frame("prop.var" = rep(NA_real_,m), "k.qform" = rep(NA_integer_,m), "qform" = rep(NA_real_,m))

  obs <- c(colSums(y * matmul(y,1)))

  if(length(k)==1 && k[1] == 0){
    res$prop.var <- 0
    res$k.qform <- 0L
    if(traces$hsnorm2 <=0){
      res$qform <- NA_real_
    } else {
      a0 <- traces$hsnorm2 / traces$trace
      nu0 <- traces$trace^2 / traces$hsnorm2
      res$qform <- -pchisq(obs/a0,df = nu0,lower.tail = FALSE, log.p = TRUE)/log(10)
    }
    return(res)
  }

  n <- nrow(y)
  k <- sort(unique(c(pmin(k,floor(0.9*n))))) # so we stop at max k even if we can't get the min variance desired
  f <- function(k,args){RSpectra::eigs_sym(matmul,
                                           k = k, n = n, args = args,
                                           opts = list("ncv" = min(n, max( 4*((2*k+1)%/%4+1), 20)) ))}


  #res <- matrix(NA_real_,nrow=15,ncol=length(obs))

  unfinished <- rep(TRUE,length(obs))
  j <- 0

  while(any(unfinished)){

    j <- j+1 # advance to next k
    e <- f(k[j], 0)

    res$prop.var <- sum(e$values^2)/traces$hsnorm2
    res$k.qform <- k[j]

    # check if traces or eigendecomposition look stable
    if(traces$hsnorm2 <=0 |
       sum(e$values^2)/length(e$values) < traces$hsnorm2/length(traces$diag) |
       abs(traces$trace - sum(e$values)) > min(abs(e$values))*(n-length(e$values))){
      # the leading e$values^2 can't sum up to the estimated hsnorm2 -- eigendecomposition is too
      # unstable and so the QForm test must be dropped for all phenotypes: we exit with all NAs
      # for the bounds and the point estimate
      return(res)}

    # if we've obtained min.prop.var or we've hit the max k, return final estimates
    if(sum(e$values^2)/traces$hsnorm2 >= min.prop.var | j == length(k) | any(e$values <= 0)){

      # subset eigenvalues and eigenvectors so that all negative eigenvalues are removed
      # and if that's not a constraint, that only the min.prop.var is obtained (to guard against low-rank clade matrices)
      k.max1 <- match(TRUE,cumsum(e$values^2)/traces$hsnorm2 >= min.prop.var)
      if(is.na(k.max1)){k.max1 <- k[j]}

      k.max2 <- match(TRUE,e$values<=0)-1L
      if(k.max2 == 0){ # the first eigenvalue is negative, so we exit with all NAs
        return(res)}
      if(is.na(k.max2)){k.max2 <- k[j]}

      k.max <- min(k.max1,k.max2)

      e$values <- e$values[1:k.max]
      e$vectors <- e$vectors[,1:k.max]

      g <- SimpleCalcQFGauss(e$values,parallel.sapply)
      z2 <- crossprod(e$vectors, y)^2
      obs.qf <- c(colSums(e$values * z2))
      res$qform <- g(obs.qf) # note we don't need to do projection of Q here b/c already baked into e$vectors/values
#
#       a0 <- sum(e$values^2) / sum(e$values)
#       nu0 <- sum(e$values)^2 / sum(e$values^2)
#
#       res[6,] <- -pchisq(obs.qf/a0,df = nu0,lower.tail = FALSE, log.p = TRUE)/log(10)
#
#       u <- pchisq(z2,df = 1,lower.tail = FALSE)
#
#       res[7,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 4),function(x){getElement(x,"p.value")}))
#       res[8,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 16),function(x){getElement(x,"p.value")}))
#       res[9,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 64),function(x){getElement(x,"p.value")}))
#
#       w <- e$values
#       w <- floor(w/min(w))
#       res[10,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 4, w = w),function(x){getElement(x,"p.value")}))
#       res[11,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 16, w = w),function(x){getElement(x,"p.value")}))
#       res[12,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 64, w = w),function(x){getElement(x,"p.value")}))
#
#       w <- e$values^2
#       w <- floor(w/min(w))
#       res[13,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 4, w = w),function(x){getElement(x,"p.value")}))
#       res[14,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 16, w = w),function(x){getElement(x,"p.value")}))
#       res[15,] <- -log10(sapply(apply(u,2, ro::renyi.test,k = 64, w = w),function(x){getElement(x,"p.value")}))

      return(res)
    }

    if(is.null(neg.log10.cutoff)){next}

    # otherwise, calculate bounds
    calc.func <- SimpleCalcBounds2(traces, e, parallel.sapply)
    temp.bounds <- calc.func(obs)

    # if any of the bounds are unstable (are not finite, then go back to top of loop for more eigenvalues)
    if(any(!is.finite(temp.bounds))){next}

    # combine QForm bounds with other test results to obtain upper bounds estimates on significance
    bound.est <- rep(NA_real_,length(obs))
    if(!is.null(other.test.res)){
      for(p in 1:length(obs)){
        bound.est[p] <- fishv(c(other.test.res[[p]],temp.bounds[p,if(lower.tail){3}else{4}]))
      }
    } else {
      bound.est <- temp.bounds[,if(lower.tail){3}else{4}]
    }

    for(p in 1:length(obs)){
      if(!is.na(bound.est[p]) && bound.est[p] < neg.log10.cutoff){
        unfinished[p] <- FALSE}
    }

  }

  # first row: percent variance captured with final k
  # second row: final k used to calculate this particular bound and point estimate
  # third row: -log10 lower bound
  # fourth row: -log10 upper bound
  # fifth row: -log10 pvalue point estimate

  res
}


SimpleCalcQFGauss <- function(e.values, parallel.sapply = base::sapply){
  gauss.tcdf <- QForm::QFGauss(e.values, parallel.sapply = parallel.sapply)
  return(function(obs, lower.tail = FALSE){-gauss.tcdf(obs, lower.tail = lower.tail, log.p = TRUE)/log(10)})
}



SimpleCalcBounds2 <- function(traces,
                              e,
                              parallel.sapply = base::sapply){

  # this returns a vectorized function that returns a matrix with two columns:
  # first column: -log10 lower bound
  # second column: -log10 upper bound

  # IF THE AMOUNT OF VARIANCE LEFT IN THE REMAINDER IS ZERO OR NEGLIGIBLE IN RELATIVE OR ABSOLUTE TERMS,
  # JUST DROP THE REMAINDER AND USE A POINT ESTIMATE BASED ON THE TOP K EIGENVALUES
  if(length(e$values)==length(traces$diag) ||
     (1 - sum(e$values^2)/traces$hsnorm2) <= .Machine$double.eps^.25 ||
     traces$hsnorm2 - sum(e$values^2) <= sqrt(.Machine$double.eps)){
    # bounds are not needed because either we have all of the eigenvalues or the eigenvalues that
    # are in the remainder are non-zero due to numerical imprecision.
    f <- SimpleCalcQFGauss(e, parallel.sapply = base::sapply)
    return(function(obs,lower.tail = FALSE){
      a <- f(obs, lower.tail)
      a <- cbind(a,a)
      colnames(a) <- c("-log10_lower_bound","-log10_upper_bound")
      a})
  }


  # Calc Required Traces
  R.max.abs.eta <- min(abs(e$values))
  R.sum.eta <-  traces$trace - sum(e$values)
  R.sum.etasq <- traces$hsnorm2 - sum(e$values^2)

  # Bound Function
  tcdf <- suppressWarnings(QForm::QFGauss(e$values,parallel.sapply = parallel.sapply))
  bound.func <- QForm::QFGaussBounds(tcdf,"identity", R.max.abs.eta, R.sum.eta, R.sum.etasq)

  function(obs, lower.tail = FALSE){
    a <- -log10(bound.func(obs)[,1:2 + if(lower.tail){0}else{2}])
    if(is.vector(a)){a <- matrix(a,1,2)} # to cover case of only one obs
    colnames(a) <- c("-log10_lower_bound","-log10_upper_bound")
    a}
}




