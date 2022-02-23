

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
    r <- renyi(putative.sprig.y[-zero.weighted.sprigs])
    message(paste("# of ranks in initial Renyi transform was ",length(r$ranks)))
    w[-zero.weighted.sprigs] <- rank2gauss(rank(r$signs * r$ranks))^4
  } else {
    r <- renyi(putative.sprig.y)
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
    updated.putative.sprig.y[-zero.weighted.sprigs] <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  } else {
    updated.putative.sprig.y <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  }


  inactive.putative.sprig.y <- rep(NA_real_,n.unique.sprigs)
  temp.exclude.ind <- rep(FALSE,n.unique.sprigs)
  temp.exclude.ind[active.sprigs] <- TRUE

  if(length(zero.weighted.sprigs)){
    temp.exclude.ind <- temp.exclude.ind[-zero.weighted.sprigs]
    inactive.putative.sprig.y[-zero.weighted.sprigs] <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
  } else {
    inactive.putative.sprig.y <- inverse.renyi(r, exclude.idx = which(temp.exclude.ind))
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



test_sprigs <- function(y,sprigs,ploidy = 2L){

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


  if(ploidy == 1){

    sprigs.to.remove.sizes <- tabulate(sprigs, nbins = n.unique.sprigs)
    sprig.coefficient <- 1/sqrt(sprigs.to.remove.sizes)

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


  outlier.test.points <- 2^(0:6)

  if(length(zero.weighted.sprigs)){
    r <- ro::renyi(tilde.z[-zero.weighted.sprigs])
  } else {
    r <- ro::renyi(tilde.z)
  }

  outliers.exps <- tail(r$exps,32)
  if(length(outliers.exps) < 32){
    p_value = NA_real_
  } else {
    p_value <- ro::mpse.test(rev(outliers.exps))
  }
  temp.exclude.idx <- which(r$ranks > r$n - 32)


  star.z <- rep(NA_real_,n.unique.sprigs)

  if(length(zero.weighted.sprigs)){
    star.z[-zero.weighted.sprigs] <- ro::inverse.renyi(r, exclude.idx = temp.exclude.idx)
  } else {
    star.z <- ro::inverse.renyi(r, exclude.idx = temp.exclude.idx)
  }

  outlier.sprig.idx <- which(is.na(star.z) & !is.na(tilde.z))

  sprig.sd <- ifelse(active.sprig.ind,renyi.sprig.sd,alt.renyi.sprig.sd)

  y.updates <- sample.weights * ((star.z - tilde.z) / sprig.sd)[renyi.sprigs]
  y.updates[is.na(y.updates)] <- 0


  # Update y.resids
  ###########################################
  y.resids <- y.resids + y.updates

  # Replace outlier y.resids with iid Gaussians
  ###############################################
  outlier.idx <- which(renyi.sprigs %in% outlier.sprig.idx)
  y.resids[outlier.idx] <- rnorm(length(outlier.idx))

  list("p.value" = p_value,
       "y" = y.resids,
       "outliers.idx" = outlier.idx,
       "outliers.exps" = outliers.exps
  )
}

#' @export
TestSprigs <- function(y, sprigs, ploidy = 2L, use.forking = FALSE, nthreads = 1L){

  if(use.forking){
    res <- parallel::mclapply(as.list(as.data.frame(y)),test_sprigs,sprigs = sprigs, ploidy = ploidy)
  } else {
    res <- apply(y,2,test_sprigs, sprigs = sprigs, ploidy = ploidy, simplify = FALSE)
  }

  list("p.value" = unlist(lapply(res,function(z){getElement(z,"p.value")})),
       "y" = do.call(cbind,lapply(res,function(z){getElement(z,"y")})))
}

