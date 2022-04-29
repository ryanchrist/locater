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

fish <- function(...){
  # vectors of negative log10 p-values
  tests <- list(...)
  -pgamma(log(10) * rowSums(do.call(cbind,tests)), shape = length(tests),lower.tail = FALSE, log.p = T)/log(10)
}

fishv <- function(v){
  # single vector of negative log10 p-values
  -pgamma(log(10) * sum(v), shape = length(v),lower.tail = FALSE, log.p = T)/log(10)
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


