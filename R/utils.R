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





calc_remainder_parameters <- function(traces, sum_evalues, sum_evalues_2, abs_max_lambda){

  C1 <- (traces$trace - sum_evalues) / abs_max_lambda # mean constraint
  C2 <- (traces$hsnorm2 - sum_evalues_2) / (abs_max_lambda^2) # variance constraint

  if(C2 <= 0){
    return(list("a"=0,"b"=0,"mu"=0))
  }

  if(C1 > C2){
    a <- C2
    b <- 0
    mu <- C1 - C2
  } else if(C1 < -C2){
    a <- 0
    b <- C2
    mu <- C1 + C2
  } else {
    a <- (C2 + C1) / 2
    b <- C2 - a
    mu <- 0
  }

  return(list("a"=a,"b"=b,"mu"=mu))
}

update_res_qform <- expression({

  k_for_T <- if(is.finite(kstar)){kstar}else{length(evalues)}

  if(calc.obs.T){
    res$obs.T <- c(colSums(evalues[1:k_for_T]*(crossprod(e$vectors[,order_evalues[1:k_for_T]],y)^2)))
  }

  R_params <- calc_remainder_parameters(traces, sum_evalues, sum_evalues_2, abs_last_evalue)

  cdf <- QForm::QFGauss(f.eta = c(evalues,
                                  if(R_params$a){abs_last_evalue}else{double()},
                                  if(R_params$b){-abs_last_evalue}else{double()}),
                        delta2 = c(rep(delta2.T,k_for_T),
                                   rep(delta2.R,length(evalues) - k_for_T),
                                   if(R_params$a){R_params$a * if(is.finite(kstar)){delta2.R}else{min(delta2.R,delta2.T)}}else{double()},
                                   if(R_params$b){R_params$b * if(is.finite(kstar)){delta2.R}else{min(delta2.R,delta2.T)}}else{double()}),
                        df = c(rep(1,length(evalues)),
                               if(R_params$a){R_params$a}else{double()},
                               if(R_params$b){R_params$b}else{double()}),
                        parallel.sapply = parallel.sapply)
  temp_obs <- obs/nu - R_params$mu*abs_last_evalue*(1 + delta2.R)

  res$qform <- -cdf(temp_obs, lower.tail = FALSE, log.p = TRUE)/log(10)

  if(!is.null(attr(cdf,"tail.features")$a.l) && is.na(attr(cdf,"tail.features")$a.l)){ # left tail extrapolation failed if a.l is NA (need to check it's non-null first b/c is.na(NULL) returns logical(0))
    to_substitute <- temp_obs < attr(cdf,"tail.features")$extrapolation.point.l # find any observations that are NA because they are below the left extrapolation point
    res$qform[to_substitute] <- runif(sum(to_substitute),min = 0,max = cdf(attr(cdf,"tail.features")$extrapolation.point.l))
    res$exit.status[to_substitute] <-  3L
  }

  details$a <- R_params$a
  details$b <- R_params$b
  details$extrapolation.point.l <- attr(cdf,"tail.features")$extrapolation.point.l
})


# see TestCladeMat help
SimpleCalcBounds <- function(y,
                             matmul,
                             traces,
                             k = NULL, # vector of positive integers, if null, we just do SW approx taking k=0.
                             prop.var.goal = 0.95,
                             var.ratio.goal = 0.95,
                             stop.eval.func = NULL,
                             nu = 1, # can be a VECTOR
                             delta2.T = 0, # for now MUST be a SCALAR
                             delta2.R = 0, # for now MUST be a SCALAR
                             calc.obs.T = FALSE, # for now, just an indicator that we should return extra information
                             parallel.sapply = base::sapply){

  # if stop.eval.func is NULL, then we evaluate the quadratic form until we hit the max k or
  # we hit the var.ratio.goal or we hit the prop.var.goal

  # stop.eval.func allows us to stop evaluation of the quadratic form early.
  # it must take in a vector of p-values (or bounds on those p-values) and return
  # a bool where TRUE stops evaluation for all phenotypes


  # return.status may be
  # 0 -- "ok" -- either we've achieved prop.var.goal or var.ratio.goal or we hit the max k
  # 1 -- traces$hsnorm2 <= 0 so "no_clade_structure" and return NAs
  # 2 -- "numerically_unstable eigendecomposition" -- at some point, results from eigendecomposition were incompatible with traces
  # 3 -- estimation of left tail failed so returned p-value is based on a random uniform sampled between 0 and the left extrapolation point

  # we keep the return.status the same across all phenotypes because they're all relevant to all phenotypes.
  # even if say one phenotype can be ruled as uninteresting using k=5 but we don't run into numerical instability
  # until say k-10, it's still important for all phenotypes that numerical instability was encountered at this locus
  # so if eigendecomposition ever disagrees with the traces, we send up a flag at that locus across all phenotypes.

  # check input
  if(!is.null(stop.eval.func)){
    if(!is.function(stop.eval.func)){
      stop("stop.eval.func must be a function with either 1 or 2 arguments, see Details")
    }
    num_args <- length(formals(args(stop.eval.func)))

    if(!(num_args %in% c(1L,2L))){
      stop("stop.eval.func must be a function with either 1 or 2 arguments, see Details")
    }
  }


  m <- ncol(y)

  if(length(nu)>1 & length(nu)!=m){stop("nu must either be a single scalar or a vector of length = ncol(y)")}
  if(any(nu<=0)){stop("all nu must be >= 0")}
  if(length(delta2.T)>1){stop("for now delta2.T must be a scalar")}
  if(length(delta2.R)>1){stop("for now delta2.R must be a scalar")}

  if(any(delta2.T<0) | any(delta2.R<0)){stop("all delta2.T and delta2.R must be >= 0")}


  # Initialize Output
  res <- data.frame("qform" = rep(NA_real_,m),
                    "obs" = rep(NA_real_,m),
                    "obs.T" = rep(NA_real_,m),
                    "precise" = rep(FALSE,m),  # can these results be interpreted as precise QForm -log10 p-values?
                    "exit.status" = rep(0L,m)) # explanation of any NAs one might see in qform, eg:
                                               # num instability in trace calculation or eigendecomposition
                                               # or in running FFT / estimating right tail using QForm package.

  # Details needed to tune inflation parameters
  details <- list(
    # features of the clade matrix
    "trace" = traces$trace,
    "hsnorm2" = traces$hsnorm2,
    "evalues" = NA_real_,
    # some helper details about how the distribution was approximated and k.star defined, this information
    # could be recovered from the details above and the input to this function, but it's helpful to have them
    # precalculated
    "kstar" = NA_integer_, # will be NA unless we get a precise estimate (we have a kstar <= k that satisfies the prop.var or var.ratio goals)
    "a" = NA_real_,
    "b" = NA_real_,
    "extrapolation.point.l"= NA_real_)



  if(traces$hsnorm2 <=0){
    res$exit.status <- 1L

    attr(res,"details") <- details
    return(res)
  }

  # Calculate observed statistics
  obs <- c(colSums(y * matmul(y,0)))
  res$obs <- obs


  # Start with SW approximation
  a0 <- traces$hsnorm2 / traces$trace
  nu0 <- traces$trace^2 / traces$hsnorm2
  # here we divide by min of nu.T and nu.R in order to make this approximation as liberal as possible for initial screening
  res$qform <- -pchisq((obs/nu)/a0,df = nu0,lower.tail = FALSE, log.p = TRUE)/log(10)


  if(is.null(k) || (length(k)==1 && k==0)){
    attr(res,"details") <- details
    return(res)
  }

  if(!is.null(stop.eval.func)){
    if((num_args==1 && stop.eval.func(10^-res$qform)) || (num_args==2 && stop.eval.func(10^-res$qform,0))){

      attr(res,"details") <- details
      return(res)
    }
  }


  # try to refine approximation
  ###################################

  n <- nrow(y)

  if(!all.equal(as.integer(k),k) || any(k<0) || any(k>floor(0.9*n))){
    stop("k must be a vector of positive integers in [0,floor(0.9*n)] where n x n is the dimension of the matrix.")}
  k <- sort(unique(as.integer(k)))
  if(k[1]==0){k <- k[-1]} # drop k=0 case if included


  # define helper function and expression

  f <- function(k,args){
    RSpectra::eigs_sym(matmul,
                       k = k, n = n, args = args,
                       opts = list("ncv" = min(n, max( 4*((2*k+1)%/%4+1), 20)) ,
                                   "retvec" = calc.obs.T))
  }



  # keep checking if we can / need to improve the approximation
  for(j in 1:length(k)){

    e <- f(k[j], 0) # returned pre-sorted from largest to smallest magnitude

    order_evalues <- order(abs(e$values),decreasing = TRUE)
    evalues <- e$values[order_evalues]
    sum_evalues <- sum(evalues)
    sum_evalues_2 <- sum(evalues^2)
    abs_last_evalue <- abs(evalues[k[j]])

    # check whether returned eigenvalues look incompatible with traces (numerical instability)
    if(sum_evalues_2/k[j] < traces$hsnorm2/length(traces$diag) |
       abs(traces$trace - sum_evalues) >  abs_last_evalue*(n-k[j])){

      res$exit.status <- 2L

      attr(res,"details") <- details
      return(res)
    }


    details$evalues <- evalues

    # calculate kstar

    kstar_prop_explained <- match(TRUE,cumsum(evalues^2)>=traces$hsnorm2*prop.var.goal)
    if(is.na(kstar_prop_explained)){kstar_prop_explained <- Inf}
    kstar_var_ratio <- if(k[j] == 1){Inf} else {
      match(TRUE,(evalues[-1]/evalues[-k[j]])^2 >= var.ratio.goal)
    }
    if(is.na(kstar_var_ratio)){kstar_var_ratio <- Inf}

    kstar <- min(kstar_prop_explained,kstar_var_ratio)


    # TEST IF WE'VE ACHIEVED A PRECISE APPROX SO WE CAN STOP
    if(is.finite(kstar)){
      eval(update_res_qform)

      res$precise <- TRUE

      details$kstar <- kstar

      attr(res,"details") <- details
      return(res)
    }
    # TRY TO STOP EARLY IF NO PHENOTYPE LOOKS SIGNIFICANT

    if(!is.null(stop.eval.func)){
      eval(update_res_qform)

      if((num_args==1 && stop.eval.func(10^-res$qform)) || (num_args==2 && stop.eval.func(10^-res$qform,0))){
        attr(res,"details") <- details
        return(res)
      }
    }

  }

  # WE'VE HIT MAX K
  # so evaluate and exit with precise=FALSE.
  eval(update_res_qform)

  attr(res,"details") <- details
  return(res)


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


  # Code below is for calculating and screening based on iterative bounds and used to be a part of the
  # for loop above
  # START
  ###############
  #
  #     if(is.null(neg.log10.cutoff)){next}
  #
  #     # otherwise, calculate bounds
  #     calc.func <- SimpleCalcBounds2(traces, e, parallel.sapply)
  #     temp.bounds <- calc.func(obs)
  #
  #     # if any of the bounds are unstable (are not finite, then go back to top of loop for more eigenvalues)
  #     if(any(!is.finite(temp.bounds))){next}
  #
  #     # combine QForm bounds with other test results to obtain upper bounds estimates on significance
  #     bound.est <- rep(NA_real_,length(obs))
  #     if(!is.null(other.test.res)){
  #       for(p in 1:length(obs)){
  #         bound.est[p] <- fishv(c(other.test.res[[p]],temp.bounds[p,if(lower.tail){3}else{4}]))
  #       }
  #     } else {
  #       bound.est <- temp.bounds[,if(lower.tail){3}else{4}]
  #     }
  #
  #     for(p in 1:length(obs)){
  #       if(!is.na(bound.est[p]) && bound.est[p] < neg.log10.cutoff){
  #         unfinished[p] <- FALSE}
  #     }
  #########
  # END
  ###############

}


SimpleCalcQFGauss <- function(e.values, mu.gauss = 0, sigma.gauss = 0, parallel.sapply = base::sapply){
  gauss.tcdf <- QForm::QFGauss(e.values, sigma = sigma.gauss, parallel.sapply = parallel.sapply)
  return(function(obs, lower.tail = FALSE){-gauss.tcdf(obs-mu.gauss, lower.tail = lower.tail, log.p = TRUE)/log(10)})
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
    f <- SimpleCalcQFGauss(e$values, parallel.sapply = base::sapply)
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


extract_carriers <- function(x,i){
  if(i == length(x$order)){return(integer())} # i == n indicates i is a leaf
  i1 <- i2 <- i
  while(i1 > 0){i1 <- x$merge[i1,1]}
  while(i2 > 0){i2 <- x$merge[i2,2]}
  x$order[match(-i1,x$order):match(-i2,x$order)]
}

make_design_matrix <- function(x,cl,n){

  if(!length(cl)){
    return(Matrix::sparseMatrix(integer(),integer(),dims = c(n,0)))
  }

  i_list <- j_list <- as.list(1:length(cl))

  for(k in 1:length(cl)){
    i_list[[k]] <- extract_carriers(x,cl[k])
    j_list[[k]] <- rep(k,length(i_list[[k]]))
  }

  Matrix::sparseMatrix(unlist(i_list),unlist(j_list),dims = c(n,length(cl)))
}

dist2design <- function(d){
  n <- length(d$height) + 1L
  d_height <- d$height
  d_height[n] <- 0
  d_merge <- d$merge
  d_merge[d_merge < 0] <- n
  d_merge <- rbind(d_merge,n)
  jumps1 <- d_height - d_height[d_merge[,1]]
  jumps2 <- d_height - d_height[d_merge[,2]]

  cutoff <- 1

  cl1 <- d_merge[jumps1 >= cutoff,1]
  cl1 <- cl1[cl1!=n]
  cl2 <- d_merge[jumps2 >= cutoff,2]
  cl2 <- cl2[cl2!=n]

  X1 <- make_design_matrix(d,cl1,n)
  X2 <- make_design_matrix(d,cl2,n)
  X1 <- X1[seq(1,n,by=2),] + X1[seq(2,n,by=2),]
  X2 <- X2[seq(1,n,by=2),] + X2[seq(2,n,by=2),]

  cbind(X1,X2)
}


rdd <- function(y,X){

  if(!is.matrix(y)){
    if(!is.numeric(y)){stop("y must be a numeric vector or a matrix")}
    y <- matrix(y,ncol=1)
  }
  m <- ncol(y)
  if(!ncol(X)){return(rep(NA_real_,m))}
  X <- rdistill::std_sparse_matrix(X)
  p1 <- locater:::ortho_part_clades(X)

  p <- ncol(X)

  k_raw <- 2^floor(log2(p))

  k_filt <- min(ceiling(0.1*p),10)

  rd_log10_pval <- rep(0,m)
  for(ii in 1:m){
    rd_log10_pval[ii] <- -log10(rdistill::rdistill(y = y[,ii], x = X, l = p1, scale_col = FALSE,
                                                   filt_opts = list("method" = "thresh", "t" = qbeta(0.05,k_filt,p-k_filt+1)),
                                                   test_opts = list("k" = min(32,k_raw)))$gpval_layers)
  }
  rd_log10_pval
}

