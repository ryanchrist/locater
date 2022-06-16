#' @export
TestLoci <- function(y, pars, target.loci = 1:L(), ploidy = 2L,
                     A = NULL,
                     num.ckpts = 0L,
                     ckpt.first.locus = FALSE,
                     k = c(10,100),
                     neg.log10.cutoff = NULL,
                     verbose = FALSE,
                     use.bettermc = FALSE,
                     use.forking = FALSE,
                     forking.chunk.size = 100L, # forking is a very expensive operation so don't make this too small
                     mc.preschedule = FALSE,
                     nthreads = 1L){
  # return a list with length and names target.idx, each locater testing results

  initial.start <- proc.time()[3]


  # input checks
  #####################
  if(!is.matrix(y)){y <- as.matrix(y)}

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


  if(length(num.ckpts)!=1 || as.integer(num.ckpts)!=num.ckpts || num.ckpts < 0){
    stop("num.ckpts must be a positive integer")
  }

  if(length(ckpt.first.locus)!=1 || !is.logical(ckpt.first.locus) || (ckpt.first.locus & num.ckpts < 2)){
    stop("ckpt.first.locus must be a logical, if TRUE num.ckpts must be >= 2")
  }

  n <- nrow(y)

  # Although we're not sure if data.table::frank is multithreaded,
  # since kalis::Clades uses it inside a forked call, we take the extra precaution of
  # setting DTthreads to 1 here.  data.table is supposed to be able to detect when
  # it's being called inside a forked process and automatically set this to 1,
  # but while we have this dependency in kalis::Clades, we keep this here as an
  # extra precaution
  data.table::setDTthreads(1L)

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

    start0 <- start1 <- proc.time()[3]
    if(num.ckpts){
      Iter(fwd,pars,target.loci[t],nthreads)
    } else {
      if(fwd$l > target.loci[t]){ResetTable(fwd)}
      Forward(fwd,pars,target.loci[t],nthreads = nthreads)
    }
    Backward(bck, pars,target.loci[t], nthreads = nthreads)
    if(verbose){print(paste("Propogating to target",length(target.loci) - t + 1L,"took",proc.time()[3] - start1,"seconds."))}


    gc()

    start1 <- proc.time()[3]
    r <- Clades(fwd, bck, pars, neighbors = TRUE, use.bettermc = use.bettermc, use.forking = use.forking,
                forking.chunk.size = forking.chunk.size, mc.preschedule = mc.preschedule, nthreads = nthreads)
    if(verbose){print(paste("Calling clades at target",length(target.loci) - t + 1L,"took",proc.time()[3] - start1,"seconds."))}
    gc()


    start1 <- proc.time()[3]
    g <- t(Haps2Genotypes(QueryCache(target.loci[t]), ploidy = ploidy, method = "additive"))

    smt.res <- TestMarker(h0, g)

    sprigs <- Sprigs(r, use.forking = use.forking, nthreads = nthreads)

    if(length(na.omit(unique(sprigs))) < 32){
      sprigs.tested.ind <- FALSE
      ro.res <- list("p.value" = rep(NA_real_,ncol(y)),
                     "y" = smt.res$y)
    } else {
      sprigs.tested.ind <- TRUE
      ro.res <- TestSprigs(smt.res$y, sprigs,
                           ploidy = ploidy,
                           use.bettermc = use.bettermc,
                           use.forking = use.forking, nthreads = nthreads)
    }
    if(verbose){print(paste("Testing Marker and Sprigs at target",length(target.loci) - t + 1L,"took",proc.time()[3] - start1,"seconds."))}


    start1 <- proc.time()[3]
    r <- CladeMat(r, ploidy = 2L, sprigs.to.prune = if(sprigs.tested.ind){sprigs}else{NULL},
                  assemble = FALSE, use.bettermc = use.bettermc, use.forking = use.forking, forking.chunk.size = forking.chunk.size,
                  mc.preschedule = mc.preschedule, nthreads = nthreads)
    if(verbose){print(paste("Building CladeMat cols at target",length(target.loci) - t + 1L,"took",proc.time()[3] - start1,"seconds."))}
    gc()

    start1 <- proc.time()[3]
    r <- do.call(cbind,r)
    r <- 0.5 * (r + t(r))
    if(verbose){print(paste("Assembling CladeMat at target",length(target.loci) - t + 1L,"took",proc.time()[3] - start1,"seconds."))}
    gc()

    # f <- function(x){ matrix multiplication function to pass to TestClades instead of explicit matrix r
    #   # do matrix multiplication with r and smt.res$Q
    # }

    start1 <- proc.time()[3]
    qf.res <- TestCladeMat(ro.res$y,r,smt.res$Q, other.test.pvalues = list(smt.res$p.value, ro.res$p.value),
                           k = k,
                           neg.log10.cutoff = neg.log10.cutoff,
                           use.bettermc = use.bettermc,
                           use.forking = use.forking, nthreads = nthreads)
    if(verbose){print(paste("TestCladeMat at target",length(target.loci) - t + 1L,"took",proc.time()[3] - start1,"seconds."))}
    gc()


    res[[t]] <- data.frame("prop.var" = qf.res[1,],
                           "k" = qf.res[2,],
                           "smt" = -log10(smt.res$p.value),
                           "ro" = -log10(ro.res$p.value),
                           "qform.lower" = qf.res[3,],
                           "qform.upper" = qf.res[4,],
                           "qform" = qf.res[5,],
                           "sw" = qf.res[6,],
                           "ru4" = qf.res[7,],
                           "ru16" = qf.res[8,],
                           "ru64" = qf.res[9,],
                           "ri4" = qf.res[10,],
                           "ri16" = qf.res[11,],
                           "ri64" = qf.res[12,],
                           "rs4" = qf.res[13,],
                           "rs16" = qf.res[14,],
                           "rs64" = qf.res[15,]
                           )

    res[[t]]$fish <- fish(res[[t]]$smt,res[[t]]$ro,res[[t]]$qform, na.rm = TRUE)
    res[[t]]$fish.lower <- fish(res[[t]]$smt,res[[t]]$ro,res[[t]]$qform.lower, na.rm = TRUE)
    res[[t]]$fish.upper <- fish(res[[t]]$smt,res[[t]]$ro,res[[t]]$qform.upper, na.rm = TRUE)

    if(verbose){print(paste(length(target.loci) - t + 1L,"out of",length(target.loci),"loci tested.",proc.time()[3] - start0,"seconds required."))}
  }

  if(verbose){print(paste("Testing complete.",proc.time()[3] - initial.start,"seconds required."))}

  names(res) <- as.character(target.loci)

  res
}
