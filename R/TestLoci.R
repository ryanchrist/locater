#' @export
TestLoci <- function(y, pars, target.loci = 1:L(), ploidy = 2L,
                     A = NULL,
                     num.ckpts = 0L, ckpt.first.locus = FALSE,
                     verbose = FALSE,
                     use.forking = FALSE, nthreads = 1L){
  # return a list with length and names target.idx, each locater testing results

  start0 <- proc.time()[3]

  # input checks
  #####################
  if(is.null(N())){
    stop("haplotypes must be cached before running TestHaplotypes, see ?kalis::CacheHaplotypes")
  }

  if(!inherits(pars,"kalisParameters")){
    stop("pars must be a kalisParameters object, use kalis::Parameters to create one")
  }

  if(!is.vector(target.loci) || !all(as.integer(target.loci)==target.loci & target.loci > 0 & target.loci <= L()) || anyDuplicated(target.loci)){
    stop("target.loci must be a vector of non-duplicated integers in [1,L()], see ?kalis::L")
  }

  target.loci <- sort(target.loci)


  if(length(num.ckpts)!=1 || as.integer(num.ckpts)!=num.ckpts || num.ckpts <= 0){
    stop("num.ckpts must be a positive integer")
  }

  if(length(ckpt.first.locus)!=1 || !is.logical(ckpt.first.locus) || (ckpt.first.locus & num.ckpts < 2)){
    stop("ckpt.first.locus must be a logical, if TRUE num.ckpts must be >= 2")
  }

  n <- nrow(y)

  # Fit null models to y
  #######################

  h0 <- FitNull(y,A)

  # Initialize Tables / Iterator
  ###########################################

  if(verbose){print(paste("Initializing tables and checkpoints..."))}

  fwd <- MakeForwardTable(pars)
  bck <- MakeBackwardTable(pars)

  if(length(target.loci) == 1){ num.ckpts <- 0L}

  # TODO consider this option of just making the Forward go to the first locus twice if we are trying to iterate over all loci.
  # if(target.loci[1] == 1){ # what to do with first locus if trying to leave out core marker
  #   target.loci[1] <- 2
  # }


  if(length(target.loci) > 1 & num.ckpts){ # Use a ForwardIterator
    if(ckpt.first.locus){
      fwd.baseline <- MakeForwardTable(pars)
      Forward(fwd.baseline,pars,target.loci[1],nthreads)
      suppressMessages(Iter <- ForwardIterator(pars, num.ckpts - 1, target.loci, fwd.baseline,force.unif = TRUE))# we take away a checkpoint here to store the baseline
    } else {
      suppressMessages(Iter <- ForwardIterator(pars, num.ckpts, target.loci,force.unif = TRUE))# we take away a checkpoint here to store the baseline
    }
  }


  # Loop Over Loci
  ###########################################
  res <- as.list(rep(NA,length(target.loci)))
  if(verbose){print(paste("Starting loop over target loci..."))}

  for(t in length(target.loci):1){
    start1 <- proc.time()[3]
    if(num.ckpts){
      Iter(fwd,pars,target.loci[t],nthreads)
    } else {
      if(fwd$l > target.loci[t]){ResetTable(fwd)}
      Forward(fwd,pars,target.loci[t],nthreads = nthreads)
    }

    Backward(bck, pars,target.loci[t], nthreads = nthreads)

    r <- Clades(fwd, bck, pars, neighbors = TRUE, use.forking = use.forking, nthreads = nthreads)

    g <- t(Haps2Genotypes(QueryCache(target.loci[t]), ploidy = ploidy, method = "additive"))

    smt.res <- TestMarker(h0, g)

    sprigs <- Sprigs(r, use.forking = use.forking, nthreads = nthreads)

    if(length(na.omit(unique(sprigs))) < 32){
      ro.res <- list("p.value" = rep(NA_real_,ncol(y)),
                     "y" = smt.res$y)
    } else {
      ro.res <- TestSprigs(smt.res$y,sprigs)
    }

    r <- CladeMat(r, ploidy = 2L, sprigs.to.prune = sprigs, assemble = FALSE, use.forking = use.forking, nthreads = nthreads)
    r <- do.call(cbind,r)
    r <- 0.5 * (r + t(r))

    # f <- function(x){ matrix multiplication function to pass to TestClades instead of explicit matrix r
    #   # do matrix multiplication with r and smt.res$Q
    # }

    qf.res <- TestCladeMat(ro.res$y,r,smt.res$Q, other.test.pvalues = list(smt.res$p.value, ro.res$p.value),
                           use.forking = use.forking, nthreads = nthreads)

    res[[t]] <- list("smt.p.value" = smt.res$p.value,
                     "ro.p.value" = ro.res$p.value,
                     "qf.res" = qf.res)

    if(verbose){print(paste(length(target.loci) - t + 1L,"out of",length(target.loci),"loci tested.",proc.time()[3] - start1,"seconds required."))}
  }

  if(verbose){print(paste("Testing complete.",proc.time()[3] - start0,"seconds required."))}

  res
}
