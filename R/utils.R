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


# see TestCladeMat help
SimpleCalcBounds <- function(y,
                             matmul,
                             traces,
                             k = NULL, # vector of positive integers, if null, we just do SW approx taking k=0.
                             min.prop.var = 0.98,
                             var.ratio.goal = 0.95,
                             stop.eval.func = NULL,
                             cs.approx = FALSE,
                             parallel.sapply = base::sapply){

  # if stop.eval.func is NULL, then we evaluate the quadratic form until we hit the max k or
  # we hit the var.ratio.goal or we hit the min.prop.var

  # stop.eval.func allows us to stop evaluation of the quadratic form early.
  # it must take in a vector of p-values (or bounds on those p-values) and return
  # a bool where TRUE stops evaluation for all phenotypes


  # return.status may be
  # 2 -- "no_clade_structure" and return NAs
  # 1 -- "numerically_unstable" -- at some point, results from eigendecomposition were incompatible with traces
  # 0 -- "ok" -- either we've achieved min.prop.var or var.ratio.goal or we hit the max k
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



  # Initialize Output
  m <- ncol(y)
  res <- data.frame("prop.var" = rep(NA_real_,m),
                    "var.ratio" = rep(NA_real_,m),
                    "k.qform" = rep(NA_integer_,m),
                    "qform" = rep(NA_real_,m),
                    "exit.status" = rep(0L,m),
                    "precise" = rep(FALSE,m)) # num instability in trace calculation or eigendecomposition

  if(traces$hsnorm2 <=0){
    res$exit.status <- 2L
    return(res)
  }

  obs <- c(colSums(y * matmul(y,0)))

  # Start with SW approximation
  res$prop.var <- 0
  res$k.qform <- 0L
  a0 <- traces$hsnorm2 / traces$trace
  nu0 <- traces$trace^2 / traces$hsnorm2
  res$qform <- -pchisq(obs/a0,df = nu0,lower.tail = FALSE, log.p = TRUE)/log(10)

  if(is.null(k) || (length(k)==1 && k==0)){
    return(res)
  }

  if(!is.null(stop.eval.func)){
    if((num_args==1 && stop.eval.func(10^-res$qform)) || (num_args==2 && stop.eval.func(10^-res$qform,0))){
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
    e_values <- RSpectra::eigs_sym(matmul,
                                   k = k, n = n, args = args,
                                   opts = list("ncv" = min(n, max( 4*((2*k+1)%/%4+1), 20)) ,
                                               "retvec" = FALSE))$values
    e_values[order(abs(e_values),decreasing = TRUE)]
  }

  update_res_qform <- expression({

    mu.R <- traces$trace - sum_evalues
    var.R <- 2*(traces$hsnorm2 - sum_evalues_2)

    if(var.R <= 0){
      g <- SimpleCalcQFGauss(e_values,
                             parallel.sapply)

    } else {

      if(cs.approx){

        new.cs.approx <- FALSE

        if(new.cs.approx){

          # here we assume that k is at least 3 here.

          x <- abs(tail(e_values,3))

          if(x[3] >= x[1]){
            # assume eigenvalues have reached plateau based on x[1]
            lambda <- double()
            abs_last_lambda <- x[1]
            last_df <- floor(var.R / (2*abs_last_lambda^2))

          } else if(x[3] >= x[2]){
            # assume eigenvalues have reached plateau based on x[2]
            lambda <- double()
            abs_last_lambda <- x[2]
            last_df <- floor(var.R / (2*abs_last_lambda^2))

          } else {

            if(x[2] >= 0.5*(x[3]+x[1]) ){
              # assume eigenvalues are exponentially decaying based on x[3] and x[1] with floor of zero
              # extrapolate out 30 eigenvalues
              lambda <- (x[3])*sqrt(x[3]/x[1])^(1:30)

            } else {
              # estimate exponential decay and floor
              alpha = max(0,(x[2]^2 - x[1]*x[3])/(2-x[1]-x[3])) # alpha must be >= 0
              beta = min(1-sqrt(.Machine$double.eps),max(sqrt(.Machine$double.eps),(x[3]-alpha)/(x[2]-alpha)))
              # beta must be in (0,1) and we back it off zero and 1 just a bit to have some numerical precision and robustness

              # extrapolate out 30 eigenvalues
              lambda <- (x[3]-alpha)*beta^(1:30) + alpha
            }

            if(var.R < 2*lambda[1]^2){

              lambda <- double() # lambda will not quite capture all of the variance
              abs_last_lambda <- 0
              last_df <- 0L

            } else {

              lambda <- lambda[1:(match(TRUE,var.R < 2*cumsum(lambda^2),nomatch = length(lambda)+1L)-1L)]
              abs_last_lambda <- lambda[length(lambda)]
              last_df <- max(0L,floor((var.R - 2*sum(lambda^2))/(2*abs_last_lambda^2)))

            }

          }


          # ASSIGN ANY REMAINING VARIANCE TO GAUSSIAN
          # ATTENTION: sigma.gauss here IS allowed to be zero. mu.gauss directly shifts the distribution, so the
          # the Gaussian component can be degenerate and the mean shift still occurs
          sigma.gauss <- max(0,var.R - 2*sum(lambda^2) - 2*last_df*abs_last_lambda^2)



          # ASSIGN SIGNS TO UNSIGNED EIGENVALUES
          if(length(lambda) | last_df){

            imputed_e_values <- c(rep(abs_last_lambda,ceiling(last_df/2)), lambda[seq(1,length(lambda),by=2)],
                                  rep(abs_last_lambda,last_df-ceiling(last_df/2)), lambda[seq(2,length(lambda),by=2)])

            cs_imputed_e_values <- c(0,cumsum(imputed_e_values))

            set_to_make_positive <- seq_len(findInterval(mu.R, cs_imputed_e_values - rev(cs_imputed_e_values)))

            imputed_e_values <- -imputed_e_values
            imputed_e_values[set_to_make_positive] <- -imputed_e_values[set_to_make_positive]

            sum_imputed_e_values <- sum(imputed_e_values)

            extended_e_values <- c(e_values,imputed_e_values)
          } else {
            extended_e_values <- e_values
          }

          browser()
          #print(extended_e_values)

          # MAKE ANY FINAL CORRECTIONS TO THE MEAN
          mu.gauss <- traces$trace - sum(extended_e_values)

          #
          #
          #   # Calculate size of each "resource" to draw on for negative eigenvalues
          #
          #   odd_lambda <- lambda[seq(1,length(lambda),by=2)]
          #   even_lambda <- lambda[seq(2,length(lambda),by=2)]
          #
          #   A <- ceiling(last_df/2)*abs_last_lambda
          #   B <- if(length(lambda)){sum(odd_lambda)}else{0}
          #   C <- (last_df - A)*abs_last_lambda
          #   D <- if(length(lambda)>=2){sum(even_lambda)}else{0}
          #
          #   sign_strategy <- findInterval(mu.R, c(-A-B-C-D, -A-B-C+D, -A-B+C+D, -A+B+C+D, A+B+C+D))
          #
          #   if(sign_strategy == 0){
          #     # flip all signs
          #     lambda <- -lambda
          #     a <- 0
          #     b <- last_df
          #   } else if(sign_strategy == 1){
          #
          #     lambda[seq(1,length(lambda),by=2)]
          #     a <- 0
          #     b <- last_df
          #   } else if(sign_strategy == 2){
          #
          #
          #   } else if (sign_strategy == 3){
          #
          #     match(TRUE, mu.R-C-D rev(odd_lambda)
          #     lambda[1:(match(TRUE,var.R < 2*cumsum(lambda^2),nomatch = length(lambda)+1L)-1L)]
          #     lambda[seq(1,length(lambda),by=2)]
          #
          #     b <- A
          #     a <- last_df - b
          #
          #   } else if (sign_strategy == 4){
          #
          #     b <- ceiling( (mu.R-B-C-D) / abs_last_lambda)
          #     a <- last_df - b
          #
          #   } else if (sign_strategy == 5){
          #     a <- last_df
          #     b <- 0
          #   }
          #
          # }



        } else {

          lambda <- abs(e_values[length(e_values)])

          a.plus.b <- floor(var.R / (2*lambda^2))

          if(a.plus.b==0){

            extended_e_values <- e_values
            mu.gauss <- mu.R
            sigma.gauss <- sqrt(var.R)

          }else{

            target_mean_in_lambda <- mu.R / lambda

            if(target_mean_in_lambda >= a.plus.b){
              a <- a.plus.b
              b <- 0
            } else if(target_mean_in_lambda <= -a.plus.b){
              a <- 0
              b <- a.plus.b
            } else {
              if(a.plus.b%%2){
                # map target_mean_in_lambda to nearest odd integer (excluding 0)
                a.minus.b <- sign(target_mean_in_lambda)*(2*pmax(1,round((abs(target_mean_in_lambda)+1)/2))-1)
              }else {
                # map target_mean_in_lambda to nearest even integer (possibly 0)
                a.minus.b <- sign(target_mean_in_lambda)*(2*round(abs(target_mean_in_lambda)/2))
              }

              a <- (a.plus.b + a.minus.b) / 2
              b <- a.plus.b - a
            }

            print(paste("Using remainder approx. with parameters",a,",",b))

            extended_e_values <- c(e_values,rep(lambda,a),rep(-lambda,b))
            sigma.gauss <- sqrt(var.R - 2*(a+b)*lambda^2)
            mu.gauss <- mu.R - lambda*(a-b)
          }

        }
      } else {
        extended_e_values <- e_values
        mu.gauss <- mu.R
        sigma.gauss <- sqrt(var.R)
      }

      if(sigma.gauss > 0){
        g <- SimpleCalcQFGauss(extended_e_values,
                               mu.gauss = mu.gauss,
                               sigma.gauss = sigma.gauss,
                               parallel.sapply)
      } else {
        g <- SimpleCalcQFGauss(e_values,
                               parallel.sapply)
      }
    }

    res$qform <- g(obs)
  })


  # keep checking if we can / need to improve the approximation
  for(j in 1:length(k)){

    e_values <- f(k[j], 0) # returned pre-sorted from largest to smallest magnitude
    sum_evalues <- sum(e_values)
    sum_evalues_2 <- sum(e_values^2)

    # check whether returned eigenvalues look incompatible with traces (numerical instability)
    if(sum_evalues_2/k[j] < traces$hsnorm2/length(traces$diag) |
       abs(traces$trace - sum_evalues) > abs(e_values[k[j]])*(n-k[j])){
      res$exit.status <- 1L
      return(res)
    }

    # accept eigenvalues
    res$prop.var <- sum_evalues_2/traces$hsnorm2
    if(k[j] > 1){res$var.ratio <- (e_values[k[j]]/e_values[k[j]-1])^2}
    res$k.qform <- k[j]

    # TEST IF WE'VE ACHIEVED A PRECISE APPROX SO WE CAN STOP
    if(res$prop.var[1] >= min.prop.var | (k[j]>1 && res$var.ratio[1] >= var.ratio.goal)){
      eval(update_res_qform)
      res$precise <- TRUE
      return(res)
    }

    # TRY TO STOP EARLY IF NO PHENOTYPE LOOKS SIGNIFICANT

    if(!is.null(stop.eval.func)){
      eval(update_res_qform)
      if((num_args==1 && stop.eval.func(10^-res$qform)) || (num_args==2 && stop.eval.func(10^-res$qform,0))){
        return(res)
      }
    }

  }

  # WE'VE HIT MAX K
  # so evaluate and exit with precise=FALSE.
  eval(update_res_qform)
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

