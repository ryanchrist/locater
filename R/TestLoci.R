
#' @import data.table

make.call.clade <- function(test.opts){

  # Add more input checks for test.opts
  if(!is.list(test.opts)){
    stop("test.opts must be a data.frame or a list")}
  test.opts <- as.data.frame(test.opts)

  # declare default values for all options; must be sorted into the three categories below
  default.opts <- list(
    # clade calling options
    "thresh" = 0.2,
    "max1var" = TRUE,
    "old.sprigs" = FALSE,
    # pre clade calling options (SMT)
    "smt.noise" = "raw",
    # clade testing options
    "max.k" = NULL,
    "sw.thresh" = 6,
    "eig.thresh" = 6.5,
    "calc.obs.T" = FALSE
  )

  if(!length(test.opts)){test.opts <- as.data.frame(default.opts)}

  if(anyDuplicated(test.opts)){stop("test.opts cannot contain duplicated rows.")}

  # some pre-processing
  if(length(setdiff(names(test.opts),names(default.opts)))){
    stop("There are some invalid options in test.opts")}

  # Fill in missing default.opts
  for(i in 1:length(default.opts)){
    if(!names(default.opts)[i]%in%names(test.opts)){
      test.opts[names(default.opts)[i]] = default.opts[[i]] }}
  # order test.opts as in default.opts and sort by options in that order
  test.opts <- test.opts[names(default.opts)]
  test.opts <- test.opts[do.call(order,test.opts),]
  test.opts["test.config"] <- as.integer(rownames(test.opts))
  test.opts <- data.table::as.data.table(test.opts)

  call.clade <- list()
  call.clade.opts <- unique(test.opts[,c("thresh","max1var","old.sprigs")])

  for(i in 1:nrow(call.clade.opts)){
    call.clade[[i]] <- list("opts" = call.clade.opts[i,],
                            "pre.clade" = list())

    # MUST UPDATE SUBSET TO MATCH clade calling options in default.opts
    call.temp.opts <- subset(test.opts,
                             thresh == call.clade.opts$thresh[i] &
                               max1var == call.clade.opts$max1var[i] &
                               old.sprigs == call.clade.opts$old.sprigs[i])
    pre.clade.opts <- unique(call.temp.opts[,c("smt.noise")])

    for(j in 1:nrow(pre.clade.opts)){

      call.clade[[i]]$pre.clade[[j]] <- list("opts" = pre.clade.opts[j,],
                                             "test.clade" = list())

      # MUST UPDATE SUBSET TO MATCH pre clade calling options in default.opts
      pre.temp.opts <- subset(call.temp.opts,
                              smt.noise == pre.clade.opts$smt.noise[j])

      test.clade.opts <- unique(pre.temp.opts[,c("max.k","sw.thresh","eig.thresh","calc.obs.T")])

      for(k in 1:nrow(test.clade.opts)){
        call.clade[[i]]$pre.clade[[j]]$test.clade[[k]] <- list(
          "opts" = test.clade.opts[k,],

          # MUST UPDATE SUBSET TO MATCH test clade options in default.opts
          "test.config" = subset(pre.temp.opts,
                                 max.k == test.clade.opts$max.k[k] &
                                   sw.thresh == test.clade.opts$sw.thresh[k] &
                                   eig.thresh == test.clade.opts$eig.thresh[k] &
                                   calc.obs.T == test.clade.opts$calc.obs.T[k])$test.config)
      }
    }
  }
  list(test.opts, call.clade)
}

# test.opts <- data.frame(
#   "smt.noise" = c(TRUE,FALSE,TRUE,TRUE,FALSE,TRUE,FALSE), # Clade-free testing options (eg: SMT, might be more complex)
#   "thresh" = c(0.2,0.2,1,0.8,0.4,0.4,1), # Clade calling options
#   "max1var" = c(rep(TRUE,3),rep(FALSE,4)),
#   "max.k" = c(128), # Clade testing options
#   "old.sprigs" = c(TRUE,TRUE,TRUE,FALSE,FALSE,FALSE,TRUE)
# )
#
# # call.clade <- make.call.clade(test.opts)

calc_k_sched <- function(x){
  if(is.null(x)){return(NULL)}
  if(x<8){return(x)}
  xx <- 2*4^(1:floor(log(x/2)/log(4)))
  if(tail(xx,1)<x){xx<-c(xx,x)}
  xx
}

#' Wrapper to Test a Chromosome with LOCATER and kalis
#'
#' A simple wrapper function designed to allow users to easliy run LOCATER in conjunction with kalis.
#'
#' @param y a \code{n} x \code{m} matrix of m quantitative phenotypes (one phenotype / column)
#' @param pars a \code{kalisParameters} object specifing HMM parameters to the Li and Stephens HMM as returned by \link{\code{kalis::Parameters}}
#' @param target.loci a vector of integers specifying index of variants among those cached with \link{\code{kalis::CacheHaplotypes}} to be screened with LOCATER
#' @param A a \code{n} x \code{q} matrix of q background covariates (an intercept -- column of 1s -- should be included). Just an intercept by default.
#' @param test.opts a list of testing options that can customized the behavior of \code{TestLoci} (further documentation coming)
#' @param verbose a logical, if TRUE print progress statements as \code{TestLoci} progresses
#' @param num.ckpts an integer, the number of checkpoints that can be stored in memory. default = 0. Increasing even to just `2L` or `3L` yields significant speed ups.
#' @param ckpt.first.locus a logical, should the first target locus be stored as a checkpoint. Default is FALSE but if the first target locus is very far from start of the cached chromosome, TRUE may yield an acceleration.
#' @param use.forking a logical, is multiprocessing by process forking allowed?  Some relatively minor acclerations are possible if TRUE, but users should verify that it is safe to launch forked processes on their compute cluster.
#' @param nthreads an integer, number of threads available
#' @return a data.table with one row per target locus including the following key fields
#' \itemize{
#' \item \code{locus.idx}: the index of corresponding variant (corresponding to \code{target.loci})
#' \item \code{phenotype}: the phenotype tested (given a column index of \code{y} )
#' \item \code{tot}: LOCATER -log10 p-value
#' \item \code{smt}: -log10 p-value from single marker testing
#' \item \code{rd}: -log10 p-value returned by running Stable Distillation
#' \item \code{qform}: -log10 p-value returned by quadratic form testing
#' }
#' @export
TestLoci <- function(y, # test phenotypes y
                     pars, # with HMM parameters pars
                     target.loci = 1:L(), # at loci target.loci
                     A = NULL, # with background covariates A
                     test.opts = list(), # testing options, may be a data.frame with more than one setting
                     # Accelerating testing
                     ############################################
                     verbose = FALSE, # if TRUE, print them. For future: If a directory as a string rather than TRUE/FALSE directory, write timings to a directory.
                     num.ckpts = 0L,
                     ckpt.first.locus = FALSE,
                     use.forking = FALSE,
                     nthreads = 1L){


  start0 <- proc.time()[3]

  # validate inputs
  #######################################
  if(!is.matrix(y)){y <- as.matrix(y)}

  if(is.null(N())){
    stop("haplotypes must be cached before running TestHaplotypes, see ?kalis::CacheHaplotypes")}

  if(!inherits(pars,"kalisParameters")){
    stop("pars must be a kalisParameters object, use kalis::Parameters to create one")}

  if(!is.vector(target.loci) || !all(as.integer(target.loci)==target.loci & target.loci > 0 & target.loci <= L()) || anyDuplicated(target.loci)){
    stop("target.loci must be a vector of non-duplicated integers in [1,L()], see ?kalis::L")}

  target.loci <- sort(target.loci)

  # Fit null models to y, validity of A is assessed here as well
  h0 <- FitNull(y,A)
  m <- ncol(y)


  # Check and arrange testing options for efficient execution
  #################################################################

  new.test.opts <- make.call.clade(test.opts)

  test.opts <- new.test.opts[[1]]
  call.clade <- new.test.opts[[2]]


  if(!is.logical(verbose) || length(verbose)!=1){stop("verbose must be TRUE or FALSE")}

  if(length(num.ckpts)!=1 || as.integer(num.ckpts)!=num.ckpts || num.ckpts < 0){
    stop("num.ckpts must be a non-negative integer")}

  if(!is.logical(ckpt.first.locus) || length(ckpt.first.locus)!=1 || (ckpt.first.locus & num.ckpts < 2)){
    stop("ckpt.first.locus must be a logical, if TRUE num.ckpts must be > 1")}

  if(!is.logical(use.forking) || length(use.forking)!=1){stop("use.forking must be TRUE or FALSE")}

  nthreads <- as.integer(nthreads)
  if(length(nthreads)!=1 || nthreads < 1){stop("nthreads must be a positive integer")}



  # Initialize Tables, Iterator, and M
  ###########################################
  if(verbose){print(paste("Initializing tables and checkpoints..."))}

  fwd <- MakeForwardTable(pars)
  bck <- MakeBackwardTable(pars)

  if(length(target.loci) == 1){num.ckpts <- 0L}

  if(length(target.loci) > 1 & num.ckpts){ # Use a ForwardIterator
    if(ckpt.first.locus){
      fwd.baseline <- MakeForwardTable(pars)
      Forward(fwd.baseline,pars,target.loci[1],nthreads)
      suppressMessages(Iter <- ForwardIterator(pars, num.ckpts - 1, target.loci, fwd.baseline,force.unif = TRUE))# we take away a checkpoint here to store the baseline
    } else {
      suppressMessages(Iter <- ForwardIterator(pars, num.ckpts, target.loci,force.unif = TRUE))# we take away a checkpoint here to store the baseline
    }}

  M <- matrix(0,N()/2,ncol(fwd$alpha)/2)


  # Loop over target loci
  ###########################################
  template.res <- test.opts
  template.res[,c("num.sprigs","k","exit.status")] <- NA_integer_
  template.res[,c("precise")] <- FALSE
  template.res[,c("obs.qform","obs.qform.T","smt","rd","qform")] <- NA_real_
  template.res <- tidyr::expand_grid(template.res,"phenotype" = if(is.null(colnames(y))){1:m}else{colnames(y)})
  res <- replicate(length(target.loci),template.res,simplify = FALSE)

  clade.details <- as.list(seq_len(length(target.loci)))

  if(verbose){print(paste("Starting loop over",length(target.loci),"target loci..."))}

  for(t in length(target.loci):1){

    # Propagate tables
    ############################
    start1 <- proc.time()[3]

    if(num.ckpts){
      Iter(fwd,pars,target.loci[t],nthreads)
    } else {
      if(fwd$l > target.loci[t]){ResetTable(fwd)}
      Forward(fwd,pars,target.loci[t],nthreads = nthreads)}
    Backward(bck, pars,target.loci[t], nthreads = nthreads)
    if(verbose){print(paste("Propagating HMM to target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start1,digits=3),"seconds."))}

    # Run tests
    ############################
    start1 <- proc.time()[3]

    for(i in 1:length(call.clade)){

      # Call Clades
      start2 <- proc.time()[3]
      neigh <- CladeMat(fwd,bck,M,unit.dist = -log(pars$pars$mu),thresh = call.clade[[i]]$opts$thresh, max1var = call.clade[[i]]$opts$max1var, nthreads = nthreads)
      if(verbose){print(paste("Calling CladeMat @ target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start2,digits=3),"seconds."))}

      sprigs <- Sprigs(neigh[[1]], old.sprigs = call.clade[[i]]$opts$old.sprigs)
      PruneCladeMat(M,neigh,sprigs,prune="singleton.info")
      PruneCladeMat(M,neigh,sprigs,prune="sprigs")

      gc()
      M <- Matrix::symmpart(M)
      gc()
      pre.clade <- call.clade[[i]]$pre.clade

      for(j in 1:length(pre.clade)){

        # Run Pre-Clade Routine (SMT)
        g <- t(Haps2Genotypes(QueryCache(target.loci[t]), ploidy = 2L, method = "additive"))
        smt.res <- TestMarker(h0, g, add.noise = pre.clade[[j]]$opts$smt.noise)

        test.clade <- pre.clade[[j]]$test.clade

        for(k in 1:length(test.clade)){

          # Test Clades: Run Renyi Distillation on Sprigs
          start2 <- proc.time()[3]
          ro.res <- TestSprigs(smt.res$y, sprigs,
                               ortho = TRUE,
                               Q = smt.res$Q,
                               use.forking = use.forking)
          if(verbose){print(paste("Call TestSprigs at target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start2,digits = 3),"seconds."))}

          # Test Clades: Test Quadratic Form
          start2 <- proc.time()[3]

          qf.res <- TestCladeMat(ro.res$y,
                                 M, # could pass function rather than M here explicitly
                                 smt.res$Q,
                                 k = calc_k_sched(test.clade[[k]]$opts$max.k),
                                 stop.eval.func = function(x, prop.var){
                                   if(prop.var==0){
                                     all(msse.test(-log10(smt.res$p.value),-log10(ro.res$p.value),-log10(x),test.1.solo = TRUE) < test.clade[[k]]$opts$sw.thresh)
                                   } else {
                                     all(msse.test(-log10(smt.res$p.value),-log10(ro.res$p.value),-log10(x),test.1.solo = TRUE) < test.clade[[k]]$opts$eig.thresh)}},
                                 calc.obs.T = test.clade[[k]]$opts$calc.obs.T,
                                 use.forking = use.forking,
                                 nthreads = 1L)

          # for now, we only store the clade details for the LAST setting of call.clade, pre.clade, and test.clade
          # that is sufficient for most use cases
          clade.details[[t]] <- attr(qf.res,"details")

          if(verbose){print(paste("Run TestCladeMat @ target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start2,digits = 3),"seconds."))}


          # Store results
          temp <- 1:m + m * (test.clade[[k]]$test.config - 1L)
          res[[t]][temp,c("num.sprigs","num.layers","obs.qform","obs.qform.T","exit.status","precise","smt","rd","qform") ] <-
            cbind(sprigs$num.sprigs,ro.res$num.layers,qf.res$obs,qf.res$obs.T,
                  qf.res$exit.status,qf.res$precise,-log10(smt.res$p.value),-log10(ro.res$p.value),qf.res$qform)

        }
      }
    }

    if(verbose){print(paste("Running tests at target",length(target.loci) - t + 1L,"out of",length(target.loci),"took",signif(proc.time()[3] - start1,digits = 3),"seconds."))}
  }

  if(verbose){print(paste("Iterating over all",length(target.loci),"target loci took",signif(proc.time()[3] - start0,digits = 3),"seconds."))}

  # merge results across loci
  names(res) <- as.character(target.loci)
  res <- data.table::rbindlist(res,idcol = "locus.idx")
  res[,locus.idx:=as.integer(locus.idx)]

  # calculate total locater signal
  res[,tot:=msse.test(smt,rd,qform,test.1.solo = TRUE)]

  # this attribute is a list where the ith element is named after the locus.idx
  names(clade.details) <- as.character(target.loci)
  attr(res,"details") <- clade.details
  res
}



TestLoci_h <- function(y, # test phenotypes y
                       pars, # with HMM parameters pars
                       target.loci = 1:L(), # at loci target.loci
                       A = NULL, # with background covariates A
                       sw.approx = FALSE, # If TRUE, use Satterthwaite Approximation for all QForm tests, critical for screening
                       test.opts = list(), # testing options, may be a data.frame with more than one setting

                       # Accelerating testing
                       ############################################
                       verbose = FALSE, # if TRUE, print them. For future: If a directory as a string rather than TRUE/FALSE directory, write timings to a directory.
                       num.ckpts = 0L,
                       ckpt.first.locus = FALSE,
                       use.forking = FALSE,
                       nthreads = 1L){


  start0 <- proc.time()[3]

  # validate inputs
  #######################################
  if(!is.matrix(y)){y <- as.matrix(y)}

  if(is.null(N())){
    stop("haplotypes must be cached before running TestHaplotypes, see ?kalis::CacheHaplotypes")}

  if(!inherits(pars,"kalisParameters")){
    stop("pars must be a kalisParameters object, use kalis::Parameters to create one")}

  if(!is.vector(target.loci) || !all(as.integer(target.loci)==target.loci & target.loci > 0 & target.loci <= L()) || anyDuplicated(target.loci)){
    stop("target.loci must be a vector of non-duplicated integers in [1,L()], see ?kalis::L")}

  target.loci <- sort(target.loci)

  # Fit null models to y, validity of A is assessed here as well
  h0 <- FitNull(y,A)
  m <- ncol(y)

  if(!is.logical(sw.approx) || length(sw.approx)!=1){stop("sw.approx must be TRUE or FALSE")}

  # Check and arrange testing options for efficient execution
  #################################################################

  new.test.opts <- make.call.clade(test.opts)

  test.opts <- new.test.opts[[1]]
  call.clade <- new.test.opts[[2]]


  if(!is.logical(verbose) || length(verbose)!=1){stop("verbose must be TRUE or FALSE")}

  if(length(num.ckpts)!=1 || as.integer(num.ckpts)!=num.ckpts || num.ckpts < 0){
    stop("num.ckpts must be a non-negative integer")}

  if(!is.logical(ckpt.first.locus) || length(ckpt.first.locus)!=1 || (ckpt.first.locus & num.ckpts < 2)){
    stop("ckpt.first.locus must be a logical, if TRUE num.ckpts must be > 1")}

  if(!is.logical(use.forking) || length(use.forking)!=1){stop("use.forking must be TRUE or FALSE")}

  nthreads <- as.integer(nthreads)
  if(length(nthreads)!=1 || nthreads < 1){stop("nthreads must be a positive integer")}



  # Initialize Tables, Iterator, and M
  ###########################################
  if(verbose){print(paste("Initializing tables and checkpoints..."))}

  fwd <- MakeForwardTable(pars)
  bck <- MakeBackwardTable(pars)

  if(length(target.loci) == 1){num.ckpts <- 0L}

  if(length(target.loci) > 1 & num.ckpts){ # Use a ForwardIterator
    if(ckpt.first.locus){
      fwd.baseline <- MakeForwardTable(pars)
      Forward(fwd.baseline,pars,target.loci[1],nthreads)
      suppressMessages(Iter <- ForwardIterator(pars, num.ckpts - 1, target.loci, fwd.baseline,force.unif = TRUE))# we take away a checkpoint here to store the baseline
    } else {
      suppressMessages(Iter <- ForwardIterator(pars, num.ckpts, target.loci,force.unif = TRUE))# we take away a checkpoint here to store the baseline
    }}

  M <- matrix(0,N(),N()) # FIXME


  # Loop over target loci
  ###########################################
  template.res <- test.opts
  #template.res[,c("num.sprigs","k")] <- NA_integer_
  #template.res[,c("prop.var","smt", "rd","qform")] <- NA_real_
  template.res <- tidyr::expand_grid(template.res,"phenotype" = if(is.null(colnames(y))){1:m}else{colnames(y)})
  res <- replicate(length(target.loci),template.res,simplify = FALSE)

  if(verbose){print(paste("Starting loop over",length(target.loci),"target loci..."))}

  for(t in length(target.loci):1){

    # Propagate tables
    ############################
    start1 <- proc.time()[3]

    if(num.ckpts){
      Iter(fwd,pars,target.loci[t],nthreads)
    } else {
      if(fwd$l > target.loci[t]){ResetTable(fwd)}
      Forward(fwd,pars,target.loci[t],nthreads = nthreads)}
    Backward(bck, pars,target.loci[t], nthreads = nthreads)
    if(verbose){print(paste("Propagating HMM to target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start1,digits=3),"seconds."))}


    kalis::DistMat(fwd,bck,type = "minus.min",M,nthreads = nthreads)
    class(M) <- "matrix"

    start2 <- proc.time()[3]
    gc()
    dM <- as.dist(Matrix::symmpart(M))
    gc()
    if(verbose){print(paste("Forming dist object @ target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start2,digits=3),"seconds."))}

    # Run tests
    ############################
    start1 <- proc.time()[3]

    for(i in 1:length(call.clade)){

      # Call Clades
      start2 <- proc.time()[3]

      d <- .Call(fastcluster:::fastcluster, N(), if(call.clade[[i]]$opts$old.sprigs){3L}else{1L}, dM, NULL)
      d$height <- d$height/(-log(pars$pars$mu)/2)

      X <-  dist2design(d)

      if(verbose){print(paste("Calling CladeMat @ target",length(target.loci) - t + 1L,"took",signif(proc.time()[3] - start2,digits=3),"seconds."))}

      pre.clade <- call.clade[[i]]$pre.clade

      for(j in 1:length(pre.clade)){

        # Run Pre-Clade Routine (SMT)
        g <- t(Haps2Genotypes(QueryCache(target.loci[t]), ploidy = 2L, method = "additive"))
        smt.res <- TestMarker(h0, g, add.noise = pre.clade[[j]]$opts$smt.noise)

        test.clade <- pre.clade[[j]]$test.clade

        for(k in 1:length(test.clade)){
          ot <- rdd(smt.res$y,X)
          Xs <- as.vector(summary(Matrix::colSums(X)))
          temp <- 1:m + m * (test.clade[[k]]$test.config - 1L)
          res[[t]][temp,c("num.clades","min_clade_size","1Q_clade_size","med_clade_size","mean_clade_size","3Q_clade_size","max_clade_size","smt","rd") ] <-
            cbind(ncol(X),Xs[1],Xs[2],Xs[3],Xs[4],Xs[5],Xs[6],-log10(smt.res$p.value),ot)
        }
      }
    }

    if(verbose){print(paste("Running tests at target",length(target.loci) - t + 1L,"out of",length(target.loci),"took",signif(proc.time()[3] - start1,digits = 3),"seconds."))}
  }

  if(verbose){print(paste("Iterating over all",length(target.loci),"target loci took",signif(proc.time()[3] - start0,digits = 3),"seconds."))}

  # merge results across loci
  names(res) <- as.character(target.loci)
  res <- data.table::rbindlist(res,idcol = "locus.idx")
  res[,locus.idx:=as.integer(locus.idx)]

  # calculate total locater signal
  res[,tot:= -pnorm((qnorm(-log(10)*smt,log.p = TRUE)+qnorm(-log(10)*rd,log.p = TRUE))/sqrt(2),log.p = TRUE)/log(10)]

  res
}







#
# if(length(k)==1 && k[1] == 0){
#   res[[t]] <- data.frame("prop.var" = qf.res[1,],
#                          "k" = qf.res[2,],
#                          "smt" = -log10(smt.res$p.value),
#                          "ro" = -log10(ro.res$p.value),
#                          "sw" = qf.res[3,])
# } else {
#   res[[t]] <- data.frame("prop.var" = qf.res[1,],
#                          "k" = qf.res[2,],
#                          "smt" = -log10(smt.res$p.value),
#                          "ro" = -log10(ro.res$p.value),
#                          "qform.lower" = qf.res[3,],
#                          "qform.upper" = qf.res[4,],
#                          "qform" = qf.res[5,],
#                          "sw" = qf.res[6,],
#                          "ru4" = qf.res[7,],
#                          "ru16" = qf.res[8,],
#                          "ru64" = qf.res[9,],
#                          "ri4" = qf.res[10,],
#                          "ri16" = qf.res[11,],
#                          "ri64" = qf.res[12,],
#                          "rs4" = qf.res[13,],
#                          "rs16" = qf.res[14,],
#                          "rs64" = qf.res[15,])
#
#   res[[t]]$fish <- fish(res[[t]]$smt,res[[t]]$ro,res[[t]]$qform, na.rm = TRUE)
#   res[[t]]$fish.lower <- fish(res[[t]]$smt,res[[t]]$ro,res[[t]]$qform.lower, na.rm = TRUE)
#   res[[t]]$fish.upper <- fish(res[[t]]$smt,res[[t]]$ro,res[[t]]$qform.upper, na.rm = TRUE)
# }

