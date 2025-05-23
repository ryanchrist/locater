

new_clade_assignments <- function(clades,sprigs.1,sprigs.2,inactive.sprigs, overwrite.single.copy.carriers = TRUE){

  inactive.sprigs.1.ind <- sprigs.1 %in% inactive.sprigs
  inactive.sprigs.2.ind <- sprigs.2 %in% inactive.sprigs

  if(overwrite.single.copy.carriers){
    res <- list("assign.to.sprigs.1" = which(inactive.sprigs.2.ind & !inactive.sprigs.1.ind),
                "assign.to.sprigs.2" = which(inactive.sprigs.1.ind & !inactive.sprigs.2.ind))
  } else {
    res <- list("assign.to.sprigs.1" = which(inactive.sprigs.2.ind & !inactive.sprigs.1.ind & !is.na(sprigs.1)),
                "assign.to.sprigs.2" = which(inactive.sprigs.1.ind & !inactive.sprigs.2.ind & !is.na(sprigs.2)))

    dual.btw.inactive.idx <- which(inactive.sprigs.1.ind & inactive.sprigs.2.ind)

    if(length(dual.btw.inactive.idx)){
      assign.vec <- sample(c(FALSE,TRUE),length(dual.btw.inactive.idx),replace = TRUE)

      temp.assign.to.sprigs.1 <- c()
      temp.assign.to.sprigs.2 <- c()

      for(i in 1:length(dual.btw.inactive.idx)){
        if(assign.vec[i]){
          temp.assign.to.sprigs.1 <- c(temp.assign.to.sprigs.1,dual.btw.inactive.idx[i])
        } else {
          temp.assign.to.sprigs.2 <- c(temp.assign.to.sprigs.2,dual.btw.inactive.idx[i])
        }
      }
      res$assign.to.sprigs.1 <- c(res$assign.to.sprigs.1,temp.assign.to.sprigs.1)
      res$assign.to.sprigs.2 <- c(res$assign.to.sprigs.2,temp.assign.to.sprigs.2)
    }
  }
  res
}

solve_active_nodes <- function(x, weights){

  if(nrow(x) >= 16){
    # as a fail safe for now in case the graph component is large, we just set all nodes to active
    # and then inactivate them one at a time from highest degree to lowest degree nodes until we find a solution
    config <- rep(1L,nrow(x))
    priority.to.remove <- nrow(x) - rank(colSums(x),ties.method = "random") + 1

    for(i in 1:nrow(x)){
      config[which(priority.to.remove == i)] <- 0L
      if(sum(config * (x%*%config)) == 0){break}
    }
    return(as.logical(config))
  }


  x <- as.matrix(x)
  storage.mode(x) <- "integer"
  configs <- as.matrix(expand.grid(replicate(n = nrow(x),c(0L,1L),simplify = FALSE)))

  valid.config.ind <- rowSums(configs * (configs%*%x)) == 0

  weights <- (1:nrow(x))^4

  config.scores <- configs[valid.config.ind,] %*% weights

  as.logical(configs[which(valid.config.ind)[which.max(config.scores)],])
}

calc_sprig_phenotypes <- function(y.resids,sprigs,n.unique.sprigs,ploidy){

  if(all(is.na(sprigs))){stop("all provided sprigs are NA")}

  if(ploidy == 1){

    sprigs.to.remove.sizes <- tabulate(sprigs, nbins = n.unique.sprigs)
    sprig.coefficient <- 1/sqrt(sprigs.to.remove.sizes)


    return(list("skip.renyi" = FALSE,
                "y.sprigs" = c(tapply(y.resids,
                                      sprigs,
                                      function(x){sum(x)/sqrt(length(x))})),
                "renyi.sprigs" = sprigs,
                "non.renyi.sprigs" = integer(),
                "sprigs.to.remove.sizes" = sprigs.to.remove.sizes,
                "sprig.coefficient" = sprig.coefficient))

  }

  if(ploidy!=2){
    stop("currently only ploidy = 1 or 2 supported")
  }

  ###########################################
  ###########################################
  # Make Unambiguous Clade Assignments
  ###########################################
  ###########################################


  sprigs.1 <- sprigs[seq.int(1,length(sprigs),2)]
  sprigs.2 <- sprigs[seq.int(2,length(sprigs),2)]

  # build clade assignments based on unambiguous sprig assignment:
  # either both alleles are assigned to the same sprig or only one allele is assigned to a sprig
  renyi.sprigs <- pmax(sprigs.1,sprigs.2,na.rm = T)
  dual.samples.idx <- which(sprigs.1 != sprigs.2)
  renyi.sprigs[dual.samples.idx] <- NA_integer_

  # the resulting renyi.sprigs here may be missing some nodes/sprigs, meaning that
  # there are no samples unambiguously assigned to those nodes/sprigs (these are nodes in our initial network with weight 0)
  # instead of trying to keep those nodes in the Renyi test by randomly assigning some of their edges to them, at least for now
  # we just remove those nodes.
  # If an edge in our graph connects a zero weight node to a positively weighted node, we assign that edge (that dual sample)
  # to the positively weighted node. In the rare case that an edge in our graph connects two zero weight nodes, we simply
  # say that edge doesn't belong to a clade for Renyi testing and test those nodes/clades/sprigs in the quadratic form
  # This procedure ensures that the resulting graph will only consist of positively weighted nodes that are still connected via
  # edges to some other nodes



  # these are the sprigs/nodes without anyone unambiguously assigned to them which will need to be tested in the quadratic form
  zero.weighted.sprigs <- which(tabulate(renyi.sprigs,n.unique.sprigs) == 0)

  if(length(zero.weighted.sprigs)==n.unique.sprigs){
    warning("All sprigs had zero weight")
    return(list("skip.renyi" = TRUE))
  }



  temp.assign.list <- new_clade_assignments(renyi.sprigs,sprigs.1,sprigs.2,inactive.sprigs = zero.weighted.sprigs)

  renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]

  # Remove samples from dual.samples.idx that have now been unambiguously assigned
  dual.samples.idx <- setdiff(dual.samples.idx,union(temp.assign.list$assign.to.sprigs.1,
                                                     temp.assign.list$assign.to.sprigs.2))
  #length(dual.samples.idx)
  #dual.samples.idx[(sprigs.1[dual.samples.idx] %in% zero.weighted.sprigs & sprigs.2[dual.samples.idx] %in% zero.weighted.sprigs)]

  # Remove samples from dual.samples.idx that were edges between two zero.weighted.sprigs
  dual.samples.idx <- dual.samples.idx[!(sprigs.1[dual.samples.idx] %in% zero.weighted.sprigs & sprigs.2[dual.samples.idx] %in% zero.weighted.sprigs)]
  #length(dual.samples.idx)

  # dual.samples.idx now indicates samples that form edges between two positively weighted nodes and need to be assigned/resolved

  # In summary:
  # Q: what do you do when a sprig is in zero.weighted.sprigs (ie: all samples under a sprig are ambiguously assigned) ?
  # A: we assign NA to that sprig's putative.sprig.y, immediately remove it from Renyi testing, and test it only in the quadratic form
  # Q: what do you do when a doubleton has one homozygous sample under it?
  # A: we just keep it and test it in the Renyi


  # this should always be true: renyi.sprigs at this stage includes all sprigs that will considered as candidates from Renyi screening
  #all.equal(sort(union(na.omit(unique(renyi.sprigs)),zero.weighted.sprigs)),1:200)



  ###########################################
  ###########################################
  # Sort Dual Samples Into Clades
  ###########################################
  ###########################################

  # generate graph
  sprig.dag <- as(Matrix::sparseMatrix(
    sprigs.1[dual.samples.idx],
    sprigs.2[dual.samples.idx],dims = rep(n.unique.sprigs,2),symmetric = T),"ngCMatrix")


  # calculate node weights
  sample.weights <- rep(1L,length(renyi.sprigs))
  sample.weights[which(sprigs.1 == sprigs.2)] <- 2L

  temp.sprig.sd <- c(tapply(sample.weights^2, renyi.sprigs, function(x){sqrt(sum(x))}))
  putative.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  putative.sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  temp.sprig.y <- c(tapply(y.resids * sample.weights, renyi.sprigs, sum))
  putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  putative.sprig.y[as.integer(names(temp.sprig.y))] <- temp.sprig.y
  putative.sprig.y <- putative.sprig.y / putative.sprig.sd


  if(!setequal(zero.weighted.sprigs, which(is.na(putative.sprig.y)))){
    stop("error: zero weighted sprigs did not line up with those missing putative sprig phenotypes")
  }

  message(paste("# of putative sprigs was ",length(putative.sprig.y)))
  message(paste("# of zero weighted sprigs: ",length(zero.weighted.sprigs)))


  w <- rep(NA_real_,n.unique.sprigs)
  if(length(zero.weighted.sprigs)){
    r <- renyi(putative.sprig.y[-zero.weighted.sprigs]) #
    message(paste("# of ranks in initial Renyi transform was ",length(r$ranks)))
    w[-zero.weighted.sprigs] <- rank2gauss(rank(r$signs * r$ranks))^4
  } else {
    r <- renyi(putative.sprig.y) #
    message(paste("# of ranks in initial Renyi transform was ",length(r$ranks)))
    w<- rank2gauss(rank(r$signs * r$ranks))^4
  }

  old.y.sprigs <- putative.sprig.y
  old.y.sprigs[zero.weighted.sprigs] <- NA_real_
  #plot(putative.sprig.y,w)


  components <- find_irreducible_components(sprig.dag,keep.isolated.nodes = FALSE)
  inactive.sprigs <- active.sprigs <- replicate(length(components),integer(),simplify=FALSE)

  # Loop through components, solving for optimal assignments, making those assignments and recording nodes/sprigs left for the quadratic form to test
  for(i in 1:length(components)){
    if(length(components[[i]])==1){next}
    optimal.config <- solve_active_nodes(sprig.dag[components[[i]],components[[i]]], w[components[[i]]])
    active.sprigs[[i]] <- components[[i]][optimal.config]
    inactive.sprigs[[i]] <- components[[i]][!optimal.config]
  }

  active.sprigs <- unlist(active.sprigs)
  inactive.sprigs <- unlist(inactive.sprigs)
  non.renyi.sprigs <- c(inactive.sprigs,zero.weighted.sprigs)



  temp.assign.list <- new_clade_assignments(renyi.sprigs, sprigs.1, sprigs.2, inactive.sprigs = inactive.sprigs, overwrite.single.copy.carriers = FALSE)
  renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]

  dual.samples.idx <- setdiff(dual.samples.idx,union(temp.assign.list$assign.to.sprigs.1,
                                                     temp.assign.list$assign.to.sprigs.2))

  dual.samples.idx <- dual.samples.idx[!(sprigs.1[dual.samples.idx] %in% inactive.sprigs &
                                           sprigs.2[dual.samples.idx] %in% inactive.sprigs)]

  if(length(dual.samples.idx)){stop("Some dual samples that should be assigned/resolved remain unassigned/unresolved.")}



  temp.assign.list <- new_clade_assignments(renyi.sprigs, sprigs.1, sprigs.2, inactive.sprigs = inactive.sprigs, overwrite.single.copy.carriers = TRUE)
  new.assign.renyi.sprigs <- rep(NA_integer_,length(renyi.sprigs))
  new.assign.renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  new.assign.renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]



  ###########################################
  ###########################################
  # Update sprig-level phenotypes
  # (removing sprigs that are now inactive
  # and adding newly assigned clade mamebers)
  ###########################################
  ###########################################


  # caution: we performed the Renyi on putative.sprig.y[-zero.weighted.sprigs] but all definitions of
  # inactive.sprigs were calculated in the original reference frame (same sprig labels as sprigs.1 etc.)
  # so we need to convert the indexing before inverting the Renyi



  updated.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  temp.exclude.ind <- rep(FALSE,n.unique.sprigs)
  temp.exclude.ind[inactive.sprigs] <- TRUE

  if(length(zero.weighted.sprigs)){
    temp.exclude.ind <- temp.exclude.ind[-zero.weighted.sprigs]
    updated.putative.sprig.y[-zero.weighted.sprigs] <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind)) #
  } else {
    updated.putative.sprig.y <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind)) #
  }


  inactive.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  temp.exclude.ind <- rep(FALSE,n.unique.sprigs)
  temp.exclude.ind[active.sprigs] <- TRUE

  if(length(zero.weighted.sprigs)){
    temp.exclude.ind <- temp.exclude.ind[-zero.weighted.sprigs]
    inactive.putative.sprig.y[-zero.weighted.sprigs] <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind)) #
  } else {
    inactive.putative.sprig.y <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind)) #
  }

  updated.putative.sprig.y <- pmax(updated.putative.sprig.y,inactive.putative.sprig.y,na.rm = TRUE)



  if(any(is.na(updated.putative.sprig.y[-c(zero.weighted.sprigs)]))){
    stop("There is at least one NA in updated.putative.sprig.y that does not correspond to a sprig listed in zero.weighted.sprigs!")
  }

  # show that the Renyi procedure made on subtle changes to putative.sprig.y
  # plot(putative.sprig.y,updated.putative.sprig.y); abline(0,1)

  new.putative.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  new.temp.sprig.sd <- c(tapply(sample.weights^2, new.assign.renyi.sprigs, function(x){sqrt(sum(x))}))
  new.putative.sprig.sd[as.integer(names(new.temp.sprig.sd))] <- new.temp.sprig.sd

  new.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  new.temp.sprig.y <- c(tapply(y.resids * sample.weights, new.assign.renyi.sprigs, sum))
  new.putative.sprig.y[as.integer(names(new.temp.sprig.y))] <- new.temp.sprig.y
  new.putative.sprig.y <- new.putative.sprig.y / new.putative.sprig.sd

  if(!setequal(which(!is.na(new.putative.sprig.y)),na.omit(unique(new.assign.renyi.sprigs))) ){
    stop("The sprigs for which newly-assigned sprig level phenotypes were calculated do not match the sprigs present in new.assign.renyi.sprigs.")
  }


  new.putative.sprig.sd[is.na(new.putative.sprig.sd)] <- 0
  new.putative.sprig.y[is.na(new.putative.sprig.y)] <- 0


  y.sprigs <- (putative.sprig.sd * updated.putative.sprig.y + new.putative.sprig.sd * new.putative.sprig.y) /
    sqrt(putative.sprig.sd^2 + new.putative.sprig.sd^2)

  if(!setequal(which(is.na(y.sprigs)), c(zero.weighted.sprigs))){
    stop("There is at least one NA in sprig.y that does not correspond to a sprig listed in zero.weighted.sprigs!")
  }



  sprigs.to.remove.sizes <- tabulate(sprigs, nbins = n.unique.sprigs)
  sprigs.to.remove.sizes[non.renyi.sprigs] <- NA_integer_


  sprig.sd <- rep(NA_real_,n.unique.sprigs)
  temp.sprig.sd <- c(tapply(sample.weights^2, renyi.sprigs, function(x){sqrt(sum(x))}))
  sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  sprig.weight <- rep(NA_real_,n.unique.sprigs)
  temp.sprig.weight <- c(tapply(sample.weights, renyi.sprigs, sum))
  sprig.weight[as.integer(names(temp.sprig.weight))] <- temp.sprig.weight

  sprig.coefficient <- sprig.sd / sprig.weight



  return(list("skip.renyi" = FALSE,
              "y.sprigs" = y.sprigs,
              "renyi.sprigs" = renyi.sprigs,
              "old.y.sprigs" = old.y.sprigs,
              "non.renyi.sprigs" = non.renyi.sprigs,
              "sprigs.to.remove.sizes" = sprigs.to.remove.sizes,
              "sprig.coefficient" = sprig.coefficient))

  # demo:
  #plot(putative.sprig.y,y.sprigs); abline(0,1)
  #qqnorm(y.sprigs);abline(0,1)

}



old_test_sprigs <- function(y,sprigs, k.max, ploidy = 2L){

  if(all(is.na(sprigs))){
    return(list("p.value" = NA_real_,
                "y" = y,
                "outliers.idx" = NA_integer_,
                "outliers.exps" = NA_real_))
  }

  ###########################################
  ###########################################
  # Calc Sprig Level Phenotypes
  ###########################################
  ###########################################

  # sprigs here to start must NOT be a factor

  y.resids.ranks <- rank(y,ties.method = "random")
  y.resids <- rank2gauss(y.resids.ranks)

  n.unique.sprigs <- length(na.omit(unique(sprigs)))


  # if(ploidy == 1){
  #
  #   sprigs.to.remove.sizes <- tabulate(sprigs, nbins = n.unique.sprigs)
  #   sprig.coefficient <- 1/sqrt(sprigs.to.remove.sizes)
  #
  # }

  if(ploidy!=2){
    stop("currently only ploidy = 2 is supported")
  }

  ###########################################
  ###########################################
  # Make Unambiguous Clade Assignments
  ###########################################
  ###########################################


  sprigs.1 <- sprigs[seq.int(1,length(sprigs),2)]
  sprigs.2 <- sprigs[seq.int(2,length(sprigs),2)]

  # build clade assignments based on unambiguous sprig assignment:
  # either both alleles are assigned to the same sprig or only one allele is assigned to a sprig
  renyi.sprigs <- pmax(sprigs.1,sprigs.2,na.rm = T)
  dual.samples.idx <- which(sprigs.1 != sprigs.2)
  renyi.sprigs[dual.samples.idx] <- NA_integer_

  # the resulting renyi.sprigs here may be missing some nodes/sprigs, meaning that
  # there are no samples unambiguously assigned to those nodes/sprigs (these are nodes in our initial network with weight 0)
  # instead of trying to keep those nodes in the Renyi test by randomly assigning some of their edges to them, at least for now
  # we just remove those nodes.
  # If an edge in our graph connects a zero weight node to a positively weighted node, we assign that edge (that dual sample)
  # to the positively weighted node. In the rare case that an edge in our graph connects two zero weight nodes, we simply
  # say that edge doesn't belong to a clade for Renyi testing and test those nodes/clades/sprigs in the quadratic form
  # This procedure ensures that the resulting graph will only consist of positively weighted nodes that are still connected via
  # edges to some other nodes



  # these are the sprigs/nodes without anyone unambiguously assigned to them which will need to be tested in the quadratic form
  zero.weighted.sprigs <- which(tabulate(renyi.sprigs,n.unique.sprigs) == 0)


  temp.assign.list <- new_clade_assignments(renyi.sprigs,sprigs.1,sprigs.2,inactive.sprigs = zero.weighted.sprigs)

  renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]

  # Remove samples from dual.samples.idx that have now been unambiguously assigned
  dual.samples.idx <- setdiff(dual.samples.idx,union(temp.assign.list$assign.to.sprigs.1,
                                                     temp.assign.list$assign.to.sprigs.2))
  #length(dual.samples.idx)
  #dual.samples.idx[(sprigs.1[dual.samples.idx] %in% zero.weighted.sprigs & sprigs.2[dual.samples.idx] %in% zero.weighted.sprigs)]

  # Remove samples from dual.samples.idx that were edges between two zero.weighted.sprigs
  dual.samples.idx <- dual.samples.idx[!(sprigs.1[dual.samples.idx] %in% zero.weighted.sprigs & sprigs.2[dual.samples.idx] %in% zero.weighted.sprigs)]
  #length(dual.samples.idx)

  # dual.samples.idx now indicates samples that form edges between two positively weighted nodes and need to be assigned/resolved

  # In summary:
  # Q: what do you do when a sprig is in zero.weighted.sprigs (ie: all samples under a sprig are ambiguously assigned) ?
  # A: we assign NA to that sprig's putative.sprig.y, immediately remove it from Renyi testing, and test it only in the quadratic form
  # Q: what do you do when a doubleton has one homozygous sample under it?
  # A: we just keep it and test it in the Renyi


  # this should always be true: renyi.sprigs at this stage includes all sprigs that will considered as candidates from Renyi screening
  #all.equal(sort(union(na.omit(unique(renyi.sprigs)),zero.weighted.sprigs)),1:200)


  #########################################################
  #########################################################
  # Calculative Putative Sprig (Vertex) Level Phenotypes
  #########################################################
  #########################################################


  # calculate node weights
  sample.weights <- rep(1L,length(renyi.sprigs))
  sample.weights[which(sprigs.1 == sprigs.2)] <- 2L

  temp.sprig.sd <- c(tapply(sample.weights^2, renyi.sprigs, function(x){sqrt(sum(x))}))
  putative.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  putative.sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  temp.sprig.y <- c(tapply(y.resids * sample.weights, renyi.sprigs, sum))
  putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  putative.sprig.y[as.integer(names(temp.sprig.y))] <- temp.sprig.y
  putative.sprig.y <- putative.sprig.y / putative.sprig.sd


  if(!setequal(zero.weighted.sprigs, which(is.na(putative.sprig.y)))){
    stop("error: zero weighted sprigs did not line up with those missing putative sprig phenotypes")
  }


  ###########################################
  ###########################################
  # Sort Dual Samples Into Clades
  ###########################################
  ###########################################

  # generate graph
  sprig.dag <- as(Matrix::sparseMatrix(
    sprigs.1[dual.samples.idx],
    sprigs.2[dual.samples.idx],dims = rep(n.unique.sprigs,2),symmetric = T),"ngCMatrix")



  rank.order <- order(abs(putative.sprig.y), na.last = TRUE, decreasing=TRUE)
  active.sprig.ind <- rep(0L,n.unique.sprigs)
  active.sprig.ind[zero.weighted.sprigs] <- NA_integer_
  for(i in seq_len(n.unique.sprigs - length(zero.weighted.sprigs))){

    current.sprig.idx <- rank.order[i]

    if(active.sprig.ind[current.sprig.idx]==-1L){next}

    neighbors <- which(sprig.dag[current.sprig.idx,])
    neighbor.status <- active.sprig.ind[neighbors]
    if(any(neighbor.status==1L)){stop("a vertex that should have been marked as inactive (-1) was left unmarked")}

    active.sprig.ind[current.sprig.idx] <- 1L
    active.sprig.ind[neighbors] <- -1L
  }
  active.sprig.ind[active.sprig.ind == -1] <- 0L

  # check solution
  # active.sprig.ind.check <- active.sprig.ind
  # active.sprig.ind.check[is.na(active.sprig.ind.check)] <- 0L
  # if((active.sprig.ind.check %*% sprig.dag %*% active.sprig.ind.check)[1,1] != 0){
  #   stop("not a valid solution!")
  # }

  active.sprigs <- which(active.sprig.ind==1)
  inactive.sprigs <- which(!active.sprig.ind)
  non.renyi.sprigs <- c(inactive.sprigs,zero.weighted.sprigs)



  alt.renyi.sprigs <- renyi.sprigs

  temp.assign.list <- new_clade_assignments(renyi.sprigs, sprigs.1, sprigs.2, inactive.sprigs = inactive.sprigs, overwrite.single.copy.carriers = FALSE)
  renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.1[temp.assign.list$assign.to.sprigs.1]
  renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.2[temp.assign.list$assign.to.sprigs.2]

  alt.renyi.sprigs[temp.assign.list$assign.to.sprigs.1] <- sprigs.2[temp.assign.list$assign.to.sprigs.1]
  alt.renyi.sprigs[temp.assign.list$assign.to.sprigs.2] <- sprigs.1[temp.assign.list$assign.to.sprigs.2]


  dual.samples.idx <- setdiff(dual.samples.idx,union(temp.assign.list$assign.to.sprigs.1,
                                                     temp.assign.list$assign.to.sprigs.2))

  # dual.samples.idx <- dual.samples.idx[!(sprigs.1[dual.samples.idx] %in% inactive.sprigs &
  #                                          sprigs.2[dual.samples.idx] %in% inactive.sprigs)]
  #
  if(length(dual.samples.idx)){stop("Some dual samples that should be assigned/resolved remain unassigned/unresolved.")}

  # See there are still some sprigs that have only one member -- that's OK, they'll be "tested" here but have little to no chance of being an outlier
  # may want to remove them in the future.
  #barplot(table(table(renyi.sprigs)))


  temp.sprig.sd <- c(tapply(sample.weights^2, renyi.sprigs, function(x){sqrt(sum(x))}))
  renyi.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  renyi.sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  temp.sprig.y <- c(tapply(y.resids * sample.weights, renyi.sprigs, sum))
  renyi.sprig.y <- rep(NA_real_,n.unique.sprigs)
  renyi.sprig.y[as.integer(names(temp.sprig.y))] <- temp.sprig.y
  renyi.sprig.y <- renyi.sprig.y / renyi.sprig.sd

  temp.sprig.sd <- c(tapply(sample.weights^2, alt.renyi.sprigs, function(x){sqrt(sum(x))}))
  alt.renyi.sprig.sd <- rep(NA_real_,n.unique.sprigs)
  alt.renyi.sprig.sd[as.integer(names(temp.sprig.sd))] <- temp.sprig.sd

  alt.y.resids <- ifelse(renyi.sprigs==alt.renyi.sprigs,y.resids,rnorm(length(y.resids)))

  temp.sprig.y <- c(tapply(alt.y.resids * sample.weights, alt.renyi.sprigs, sum))
  alt.renyi.sprig.y <- rep(NA_real_,n.unique.sprigs)
  alt.renyi.sprig.y[as.integer(names(temp.sprig.y))] <- temp.sprig.y
  alt.renyi.sprig.y <- alt.renyi.sprig.y / alt.renyi.sprig.sd

  tilde.z <- ifelse(active.sprig.ind,renyi.sprig.y,alt.renyi.sprig.y)


  if(length(zero.weighted.sprigs)){
    r <- rdistill::renyi(tilde.z[-zero.weighted.sprigs]) #
  } else {
    r <- rdistill::renyi(tilde.z) #
  }

  outliers.exps <- tail(r$exps,k.max)

  if(length(outliers.exps) < k.max){
    p_value = NA_real_
  } else {
    p_value <- rdistill::mpse.test(rev(outliers.exps))
    r$exps[(r$n - k.max + 1L):r$n] <- rexp(k.max)
  }

  star.z <- rep(NA_real_,n.unique.sprigs)

  if(length(zero.weighted.sprigs)){
    star.z[-zero.weighted.sprigs] <- rdistill::inverse.renyi(r) #
  } else {
    star.z <- rdistill::inverse.renyi(r) #
  }

  sprig.sd <- ifelse(active.sprig.ind,renyi.sprig.sd,alt.renyi.sprig.sd)

  y.updates <- sample.weights * ((star.z - tilde.z) / sprig.sd)[renyi.sprigs]
  y.updates[is.na(y.updates)] <- 0

  # Update y.resids
  ###########################################
  y.resids <- y.resids + y.updates

  list("p.value" = p_value,
       "y" = y.resids,
       #"outliers.idx" = outlier.idx,
       "outliers.exps" = outliers.exps
  )
}

old_TestSprigs <- function(y, sprigs, k.max = 2^min(max(4,ceiling(log2(sprigs$num.sprigs/100))),7) , ploidy = 2L, use.forking = FALSE, nthreads = 1L){

  if(sprigs$num.sprigs < 64){ # require there to be at least 64 sprigs to run
    return(list("p.value" = rep(NA_real_,ncol(y)),"y" = y))}

  if(use.forking){
    res <- parallel::mclapply(as.list(as.data.frame(y)),test_sprigs,sprigs = sprigs$assignments, k.max = k.max, ploidy = ploidy)
  } else {
    res <- apply(y,2,test_sprigs, sprigs = sprigs$assignments, k.max = k.max, ploidy = ploidy, simplify = FALSE)
  }

  list("p.value" = unlist(lapply(res,function(z){getElement(z,"p.value")})),
       "y" = do.call(cbind,lapply(res,function(z){getElement(z,"y")})))
}


design2layer <- function(X, col.idx = NULL){

  if(!inherits(X,"sparseMatrix")){
    stop("X must be a sparseMatrix from Matrix package, see Matrix::sparseMatrix")}

  if(!is.integer(col.idx)){col.idx <- as.integer(col.idx)}
  if(!(is.atomic(col.idx) && is.integer(col.idx))){
    stop("col.idx must be a vector of integers")}

  # if(is.null(col.idx)){
  #   a.levels <- seq_len(ncol(X))
  # }else{
  #   a.levels <- col.idx}

  p <- ncol(X)
  X <- X[,col.idx]
  #var.total <- Matrix::colSums(X^2)

  X <- X * (Matrix::rowSums(X) == 1)
  #var.assigned <- Matrix::colSums(X^2)
  #prop.var <- var.assigned / var.total

  a <- rep(NA_integer_,nrow(X))
  X <- Matrix::mat2triplet(X)
  a[X$i] <- col.idx[X$j]

  rdistill::rdlayer(a = a,p = p,ortho = FALSE)
}


part.clades <- function(x1, x2, p = max(c(max(x1,na.rm=T),max(x2,na.rm=T)),na.rm=T)){

  n <- length(x1)
  if(length(x2) != n){stop("x1 and x2 must have equal length")}

  x1[is.na(x1)] <- p+1L
  x2[is.na(x2)] <- p+1L

  # Declare adjacency matrix A between clades
  # note redundant i,j entries are ignored here
  A <- Matrix::sparseMatrix(x1,x2,
                            dims=c(p+1,p+1))
  A <- A | Matrix::t(A)
  csA <- Matrix::colSums(A)

  # we set A[p+1,p+1] <- TRUE to simplify our downstream code
  # this ensures that p+1 clade is always treated as have a homozygous sample and automatically set to the second layer.
  # even in case there are not samples with homozygous (NA,NA) alleles
  A[p+1,p+1] <- TRUE

  # now assign each clade to a first or second layer
  in.first.layer <- rep(FALSE,p+1)
  assigned <- rep(FALSE,p+1)
  edge.first.layer <- rep(FALSE,p+1)
  edge.second.layer <-  rep(FALSE,p+1)

  # Assign all homozygous samples and the p+1 clade to the second layer
  # (Samples with two NAs for their sprigs (homozygous for the p+1 clade) end up
  # assigned to second layer along with the p+1 clade.)
  assigned[Matrix::diag(A)] <- TRUE

  # initialize edge.second.layer
  edge.second.layer <-  assigned

  # see which samples are hard assigned to the second
  while(any(!assigned)){

    # if no edge.second.layer assign a high degree clade to the second layer
    if(!any(edge.second.layer)){
      temp.idx <- which.max(csA[!assigned])
      edge.second.layer[!assigned][temp.idx] <- TRUE
      assigned[!assigned][temp.idx] <- TRUE}

    # Assign any unassigned nodes (clades) that border the edge of the second layer
    # to be in the first layer
    edge.first.layer <- !assigned & Matrix::rowSums(A[,edge.second.layer,drop=FALSE]) > 0
    assigned[edge.first.layer] <- TRUE
    in.first.layer[edge.first.layer] <- TRUE

    if(any(edge.first.layer)){
      # Assign any unassigned nodes (clades) that border the edge of the first layer
      # to be in the second layer
      edge.second.layer <- !assigned & Matrix::rowSums(A[,edge.first.layer,drop=FALSE]) > 0
      assigned[edge.second.layer] <- TRUE
    } else {
      edge.second.layer <-  rep(FALSE,p+1)
    }
  }

  if(!all(assigned)){stop("Error occurred during clade partitioning: not all clades were assigned to layer 1 or layer 2")}

  # Excluding clades with a homozygote carrier and the p+1 clade, (done with !diag(A))
  # all clades assigned to one layer must have at least one neighbor in the other layer
  if(!all(Matrix::colSums(A * (!in.first.layer))[in.first.layer & !Matrix::diag(A)] > 0)){
    stop("Error occurred during clade partitioning: some clade in the 1st layer that needs a neighbor in 2nd layer doesn't have one")}
  #
  if(!all(Matrix::colSums(A * (in.first.layer))[!in.first.layer & !Matrix::diag(A)] > 0)){
    stop("Error occurred during clade partitioning: some clade in the 2nd layer that needs a neighbor in 1st layer doesn't have one")}

  return(in.first.layer[-(p+1)])
}



ortho_part_clades <- function(X){

  n <- nrow(X)
  p <- ncol(X)

  r <- order(Matrix::colSums(X))
  layers <- list(rep(NA_integer_,n))
  L <- 1

  # num_assigned <- 0L
  for(j in r){
    assigned <- FALSE
    i_range <-X@p[c(j,j+1L)]
    support_j <- X@i[seq.int(i_range[1]+1L,i_range[2])]+1L # non-zero rows of X[,j]

    for(l in 1:L){
      if(all(is.na(layers[[l]][support_j]))){
        layers[[l]][support_j] <- j
        assigned <- TRUE
        break
      }
    }

    if(!assigned){
      L <- L + 1
      layers[[L]] <- rep(NA_integer_,n)
      layers[[L]][support_j] <- j
    }
    # num_assigned <- num_assigned + 1L
    # print(num_assigned)
  }

  # pack solutions into RD layers
  for(l in 1:L){
    layers[[l]] <- rdistill::rdlayer(layers[[l]],p,ortho=TRUE)}

  # sort layers from largest to smallest
  layers[order(sapply(layers,rdistill::size),decreasing = TRUE)]
}


##load sprigs
# sprigs <- readRDS("~/Downloads/temp_sprigs.rds")
# p <- sprigs$num.sprigs
# x <- sprigs$assignments
#
# # form design matix X
# x[is.na(x)] <- p+1L
# n <- length(x)/2
# X <- Matrix::sparseMatrix(i = rep(1:n,each=2), j = x, x = 1L, dims=c(n,p+1))
# X <- X[,-(p+1)]
#
# layers <- ortho_part_clades(X)


new_test_sprigs <- function(y, x, layers, sprigs, ortho = FALSE, ...){

  k <- function(x){2^pmin(pmax(4,ceiling(log2(x/100))),7)}

  if(ortho){
    rdistill::rdistill(y = y, x = x, l = layers,
                       filt_opts = list("method" = "thresh", "t" = qbeta(0.05,8,ncol(x)-8+1)),
                       test_opts = list("k" = k))
  } else {
    rdistill::rdistill(y = y, x = x, l = layers,
                       filt_opts = list("method" = "topk","k" = k),
                       test_opts = list("k" = k))
  }

}


#' Test Sprigs
#'
#' Test Sprigs object as returned by \link{\code{kalis::Sprigs}} using Stable Distillation.
#'
#' @param y a \code{n} x \code{m} matrix of phenotype residuals, one phenotype per column
#' @param sprigs a list as returned by \link{\code{kalis::Sprigs}}
#' @param ortho a logical, to be deprecated, currently kept for background compatibility
#' @param Q a \code{n} x \code{q} orthogonal matrix with columns spanning the column space of the background covariates
#' @param use.forking a logical, is multiprocessing by process forking allowed?  Some relatively minor acclerations are possible if TRUE, but users should verify that it is safe to launch forked processes on their compute cluster.
#' @export
TestSprigs <- function(y, sprigs, ortho = FALSE, Q = matrix(1/sqrt(nrow(y)),nrow(y),1), use.forking = FALSE){

  # min number of predictors required to run Renyi Distillation at all
  if(sprigs$num.sprigs < 8){ # require there to be at least 8 sprigs to run
    return(list("p.value" = rep(NA_real_,ncol(y)),"y" = y),"num.layers"=0L)}

  p <- sprigs$num.sprigs
  x <- sprigs$assignments

  # Convert sprig assignments x into a design matrix x for clades (n samples x p predictors)
  x[is.na(x)] <- p+1L
  n <- length(x)/2
  X <- Matrix::sparseMatrix(i = rep(1:n,each=2), j = x, x = 1L, dims=c(n,p+1))
  X <- X[,-(p+1)]

  # if(!is.null(Q)){

    # if(use.forking){
    #   res <- parallel::mclapply(as.list(as.data.frame(y)),
    #                             rdistill::rdistill_pivot,x = X, Q = Q, max_num_causal = 16, ...)
    # } else {
    #   res <- apply(y, 2, rdistill::rdistill_pivot, x = X, Q = Q, max_num_causal = 16, ..., simplify = FALSE)
    # }

    # return(list("p.value" = unlist(lapply(res,function(z){getElement(z,"p_value")})),
    #             "y" = do.call(cbind,lapply(res,function(z){getElement(z,"y")})),
    #             "num.layers" = ncol(X)))

    res <- distill_pivot_par(y = y, x = X, Q = Q, max_num_causal = 16)

    return(list("p.value" = res$p_value,
                "y" = res$y,
                "num.layers" = nrow(res$u) - c(Matrix::colSums(is.na(res$u))) ))

  # } else {
  #
  #   if(ortho){
  #     layers  <- ortho_part_clades(X)
  #   } else {
  #     # set function for calculating k based on number of clades in a layer
  #     in.first.layer <- part.clades(x1 = x[seq.int(1,length(x),2)],
  #                                   x2 = x[seq.int(2,length(x),2)],
  #                                   p = p)
  #     layers <- list(design2layer(X,which(in.first.layer)),
  #                    design2layer(X,which(!in.first.layer)))
  #     # sort layers so they're tested from largest to smallest
  #     layers <- layers[order(sapply(layers,rdistill::size),decreasing = TRUE)]
  #   }
  #
  #   if(use.forking){
  #     res <- parallel::mclapply(as.list(as.data.frame(y)),new_test_sprigs, x = X, layers = layers, sprigs = sprigs, ortho = ortho, ...)
  #   } else {
  #     res <- apply(y, 2, new_test_sprigs, x = X, layers = layers, sprigs = sprigs, ortho = ortho, ..., simplify = FALSE)
  #   }
  #
  #   return(list("p.value" = unlist(lapply(res,function(z){getElement(z,"gpval_layers")})),
  #                 "y" = do.call(cbind,lapply(res,function(z){getElement(z,"y_new")})),
  #                 "num.layers" = length(layers)))
  # }


}




