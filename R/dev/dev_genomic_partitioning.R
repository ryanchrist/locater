
rm(list=ls())


design2layer <- function(X, col.idx = NULL){

  if(!inherits(X,"sparseMatrix")){
    stop("X must be a sparseMatrix from Matrix package, see Matrix::sparseMatrix")}

  if(!is.integer(col.idx)){col.idx <- as.integer(col.idx)}
  if(!(is.atomic(col.idx) && is.integer(col.idx))){
    stop("col.idx must be a vector of integers")}

  if(is.null(col.idx)){
    a.levels <- seq_len(ncol(X))
  }else{
    a.levels <- col.idx}

  X <- X[,col.idx]
  var.total <- Matrix::colSums(X^2)

  X <- X * (Matrix::rowSums(X) == 1)
  var.assigned <- Matrix::colSums(X^2)
  prop.var <- var.assigned / var.total

  a <- rep(NA_integer_,nrow(X))
  X <- Matrix::mat2triplet(X)
  a[X$i] <- col.idx[X$j]

  list("a" = a,
       "a.levels" = a.levels,
       "prop.var" = prop.var,
       "ortho" = FALSE,
       "a.anyNA"= anyNA(a))
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

#saveRDS(sprigs,"~/Downloads/temp_sprigs.rds")





part.clades2 <- function(x1, x2,
                         p = max(c(max(x1,na.rm=T),max(x2,na.rm=T)),na.rm=T),
                         min_pv_assigned = 2/3,
                         min_pv_available = 4/5){

  water_fill <- function(x,w){
    if(!is.integer(w) || any(w)<0){
      stop("w must be an integer vector with non-negative entries")}
    W <- sum(w)
    if(!is.integer(x) || length(x) != 1 || x < 0 || x> W){
      stop("x must be an integer in [0,sum(w)]")}

    w <- (x/W) * w
    fw <- floor(w)
    csr <- x - sum(fw)
    if(csr > 0){
      oo <- order(w - fw,decreasing = TRUE)[1:csr]
      fw[oo] <- fw[oo] + 1L}
    fw
  }

  hae <- function(s,e,v,d){ # hard assign edges function

    m <- length(e)
    if(length(v)!=m || length(d)!=m){
      stop("the length of e, v, and d must be equal")}

    e <- pmin(e,v)
    Se <- sum(e)
    if(s > Se){stop("s cannot be greater than pmin(e,v)")}
    if(s == Se){ return(e) }

    a <- rep(0,m)

    # hard assign any edges from nodes that have no other unassigned neighbors (free variance to hard assign)
    a[d==0] <- e[d==0]
    d[d==0] <- Inf
    v <- v - a
    e <- e - a
    s <- s - sum(a)

    while(s > 0){
      r <- (v-1)/d
      r[e==0] <- -Inf
      i <- which.max(r)
      a[i] <- a[i] + 1L
      v[i] <- v[i] - 1L
      e[i] <- e[i] - 1L
      s <- s - 1L
    }

    a
  }

  make_hard_assignments_based_on_demand_for_available_samples <- expression({

    if(any(l_neighbors)){

      neighbors_in_l <- neighbors[l_neighbors]

      num_open_edges <- A@x[idx[l_neighbors]] - D[[l]][js,neighbors_in_l]

      # var_budget[,l] =  floor(ac * (1-min_pv_available)) - rowSums(D)
      ad <- pmin(var_budget[neighbors_in_l,l],num_open_edges) # number of edges that can be directed towards js

      t <- t + sum(ad)

      if(t < tau[js]){
        # we don't have enough hard variance available so we make these assignments and move on
        proposed_D_col[neighbors_in_l,l] <- ad
        proposed_var_budget[neighbors_in_l,l] <- proposed_var_budget[neighbors_in_l,l] - ad
      }

      if(t == tau[js]){
        # we have EXACTLY enough hard variance assigned and it doesn't take any edges from unassigned nodes, no need to consider other potential layer assignments
        proposed_D_col[neighbors_in_l,l] <- ad
        proposed_var_budget[neighbors_in_l,l] <- proposed_var_budget[neighbors_in_l,l] - ad

        solution_found[l] <- TRUE
        break
      }

      if(t > tau[js]){
        # we have MORE THAN enough hard variance assigned and it doesn't take any edges from unassigned nodes, no need to consider other potential layer assignments

        # quickly figure out which edges to direct towards js

        temp_a  <- hae(shortfall,
                       num_open_edges,
                       proposed_var_budget[neighbors_in_l,l],
                       # number of undirected edges from nodes i ~ j to unassigned nodes that are not j (we subtract num_open_edges to account for this "that are not j").
                       colSums(A[!assigned,neighbors_in_l,drop=FALSE]-D[!assigned,neighbors_in_l,drop=FALSE])- num_open_edges)

        proposed_D_col[neighbors_in_l,l] <- temp_a
        proposed_var_budget[neighbors_in_l,l] <- proposed_var_budget[neighbors_in_l,l] - temp_a

        solution_found[l] <- TRUE
        break
      }

    }
  })

  find_assignment_option_for_js <- expression({

    #l <- 5 # FIXME

    if(var_budget[js,l] < 0){next}

    # ASSESS AMOUNT OF "FREE" VARIANCE FROM PREDICTORS THAT ARE EITHER
    # A) ALLOCATED TO A LAYER BESIDES L
    # B) UASSIGNED AND INELIGIBLE FOR ASSIGNMENT TO LAYER L
    # WE DONT NEED TO MAKE HARD ASSIGNMENTS FROM THESE PREDICTORS BECAUSE
    # THOSE HARD ASSIGNMENTS WILL BE PICKED UP WHILE FORMING THE DESIGN MATRIX
    # FOR THE LAYER (THE ROWS WITH ONLY ONE NON-ZERO ENTRY)
    ################################################################################
    ################################################################################

    i_range <-A@p[c(js,js+1L)]
    idx <- (i_range[1]+1L):i_range[2]
    neighbors <- A@i[idx]+1L # non-zero rows of A[,js]

    # update idx and neighbors by removing self from neighborhood (if present)
    idx <- idx[neighbors!=js]
    neighbors <- neighbors[neighbors!=js]

    a_neighbors <- a[neighbors]

    d <- A@x[idx[ (a_neighbors!=0 & a_neighbors!=l) | # neighbors that are assigned to a layer beside l
                    (a_neighbors==0 & var_budget[neighbors,l] < 0) ]] # neighbors that are not assigned and are ineligible for layer l


    t <- A[js,js] + sum(d)

    if(t>=tau[js]){
      # we have enough variance hard assigned to js without making any hard assignments, no need to consider other potential layer assignments
      # assign js to layer l
      solution_found[l] <- TRUE
      break
    }


    # TRY TO GET VARIANCE FROM OTHER PREDICTORS IN LAYER l
    # WHILE RESPECTING THEIR MIN_VAR_AVAILABLE CONSTRAINT
    ########################################################
    ########################################################

    shortfall <- tau[js] - t
    l_neighbors <- a_neighbors==l
    eval(make_hard_assignments_based_on_demand_for_available_samples)

    # TRY TO GET VARIANCE FROM UNASSIGNED PREDICTORS
    # WITHOUT MAKING ANY INELIGIBLE FOR LAYER l
    ########################################################
    ########################################################

    shortfall <- tau[js] - t
    l_neighbors <- a_neighbors==0 & var_budget[neighbors,l] >= 0
    eval(make_hard_assignments_based_on_demand_for_available_samples)


    # TRY TO GET VARIANCE FROM UNASSIGNED PREDICTORS
    # BY MAKING SOME OR ALL INELIGIBLE FOR LAYER l
    ########################################################
    ########################################################

    shortfall <- tau[js] - t
    l_neighbors <- a_neighbors==0 & var_budget[neighbors,l] >= 0
    neighbors_in_l <- neighbors[l_neighbors]
    num_open_edges <- A@x[idx[l_neighbors]] - D[[l]][js,neighbors_in_l]

    # Look to pull required hard assigned variance from unassigned nodes.
    num_options <- rowSums(var_budget[neighbors_in_l,,drop=FALSE] >= 0)
    oo <- order(num_options,decreasing = TRUE)

    num_options <- num_options[oo]
    neighbors_in_l <- neighbors_in_l[oo]
    num_open_edges <- num_open_edges[oo]

    # this is the set of samples whose edges need to be assigned to js
    select_idx <- 1:match(TRUE,cumsum(num_open_edges) >= shortfall)
    connections_to_fully_assign <- neighbors_in_l[select_idx]

    cost[l] <- sum(1/num_options[select_idx])
    proposed_D_col[connections_to_fully_assign,l] <- A[connections_to_fully_assign,js]
    proposed_var_budget[connections_to_fully_assign,l] <- -1L # make the budget simply negative (ineligible)

    solution_found[l] <- TRUE
  })


  min_pv_assigned = 2/3
  min_pv_available = 4/5

  sprigs <- readRDS("~/Downloads/temp_sprigs.rds")
  p <- sprigs$num.sprigs
  x <- sprigs$assignments
  x1 = x[seq.int(1,length(x),2)]
  x2 = x[seq.int(2,length(x),2)]

  n <- length(x1)
  if(length(x2) != n){stop("x1 and x2 must have equal length")}

  x1[is.na(x1)] <- p+1L
  x2[is.na(x2)] <- p+1L

  L <- 1
  D <- list(Matrix::sparseMatrix(1,1,x=0.0,dims=c(p,p)))
  assigned <- rep(FALSE,p)
  a <- rep(0L,p)


  # Declare adjacency matrix A between clades
  # note redundant i,j entries added together
  A <- Matrix::sparseMatrix(x1,x2, x = 1,
                            dims=c(p+1,p+1))
  A <- 2 * Matrix::symmpart(A)

  diag(A) <- 2 * diag(A) + A[,p+1]
  A <- A[-(p+1),-(p+1)]
  A <- as(A,Class = "generalMatrix")

  # quick sanity check:
  #  ac <- tabulate(x1,p) + tabulate(x2,p) # allele count of each clade
  # all.equal(ac,colSums(A))

  ac <- colSums(A)
  tau <- ceiling(min_pv_assigned * ac)
  var_budget <- matrix(floor(ac * (1-min_pv_available)),p,L)   # count of alleles available from each clade

  while(any(!assigned)){

    # select js, this is slow but OK for now.
    js <- which(!assigned)[order(rowSums(var_budget[!assigned,,drop=FALSE]>=0),ac[!assigned],decreasing = TRUE)[1]]

    #js <- 5 # FIXME

    cost <- rep(0,L)
    solution_found <- rep(FALSE,L)
    proposed_D_col <- Matrix::sparseMatrix(1,1,x=0.0,dims=c(p,L))
    proposed_var_budget <- var_budget

    # LOOP OVER EXISTING LAYERS
    for(l in 1:L){
      eval(find_assignment_option_for_js)}

    if(any(solution_found)){
      l <- which(solution_found)[which.min(cost[solution_found])]
    } else {
      # make new layer
      l <- L <- L+1
      D[[L]] <- Matrix::sparseMatrix(1,1,x=0.0,dims=c(p,p))
      var_budget <- cbind(var_budget,rep(floor(ac * (1-min_pv_available)),p))
      proposed_D_col <- cbind(proposed_D_col,Matrix::sparseMatrix(integer(),integer(),dims=c(p,1)))
      proposed_var_budget <- var_budget
      eval(find_assignment_option_for_js)
    }

    a[js] <- l
    D[[l]][,js] <- proposed_D_col[,l]
    var_budget[,l] <- proposed_var_budget[,l]

    # add in validation function here to make sure new assignment is OK.
    assigned[js] <- TRUE
  }

  list(a,D)

}



# important! :
# if an edge is directed (assigned to a js in one layer), then it's
# available to be hard assigned to the node on the other side of the edge
# in EVERY other layer (past and future layers)

#sort entries of x1 and x2 so that x1 <= x2.

# sufficient to find a solution:
# the layer assignment of each row, a vector a in {NA,1,...,L}^n
# a vector in {-1,0,1}^|edges between predictors in same layer| saying
# 0 if an edge is directed in none of the layers
# or a direction it's in if it's hard assigned: towards x1 or x2. Of course that only works
# if we assume that all edges are directed the same way when there are multiple edges
# but that seems reasonable for now...



# to match format expected by rdistill, within a given layer we need
# an assignment vector a giving the hard assignment of each row.
# from above, that's the same as subsetting X to only the columns in
# a given layer, hard assigning all rows with rowSum 1 or an entry of 2
# to the corresponding predictor column.
# then just using the hard assignment edge directions to pick out which
# of the columns to use for the rows with two 1s in them.


# maybe we need a second type of matrix where each entry gives the sample number that underlies
# a given edge -- only problem is what do you do when there are multiple edges between the same two nodes.
# maybe we make an extended key system. in our new symmetric encoding matrix M, if there's only
# one edge between nodes (predictors) i and j (i <= j), then M[i,j] gives the row number (1,...,n)
# if there are multiple edges then it gives a number starting from n+1 (or maybe even a negative number)
# which gives the list element in a list of vectors where each vector is a set of row numbers.

# It could also be that the matrix encodes a number like the p vector in the sparse matrix representation
# such that xx[p:(p+yy[p])] gives the set of row numbers.  yy here can be a sparse vector perhaps.
#  ... wait, this is just another sparse matrix structure could just be that we make M[i,j] give the row number
# of another sparse matrix, M2, such that M2[M[i,j],] has 1s indicating the samples that connect.
# M2 is # of pairs of connected edges by n (a LARGE sparse matrix)...but that's OK...same as other list structure
# but probably faster and easier to work with? M and M2 together give another means of encoding the matrix.
# we can also keep another vector b = rowSums(M2) which allow us to quickly see how m




# maybe we just end up with the sufficient info described above, and when
# we go and actually build a submatrix X for a layer l, we use it to complete
# our assignment vector a on the spot, looping through pairs of nodes that
# we know are connected with multiple edges and just finding those rows
# and hard assigning the required number of them -- won't be too expensive
# unless there are many nodes connected by multiple edges.

# so yeah lets go back and do this later,,,,





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









sprigs <- readRDS("~/Downloads/temp_sprigs.rds")
p <- sprigs$num.sprigs
x <- sprigs$assignments
x1 = x[seq.int(1,length(x),2)]
x2 = x[seq.int(2,length(x),2)]

in.first.layer <- part.clades(x1 = x[seq.int(1,length(x),2)],
                              x2 = x[seq.int(2,length(x),2)],
                              p = p)

# Need to understand if/when we have clades (eg: doubletons) where their
# entire variance is not available
x[c(2493,2494,2603,2604)]



# These entries are sample that carry at least one doubletons where only half of that doubleton's variance is available:
which(((in.first.layer[x][seq.int(1,length(x),2)] & in.first.layer[x][seq.int(2,length(x),2)]) |
         (!in.first.layer[x][seq.int(1,length(x),2)] & !in.first.layer[x][seq.int(2,length(x),2)])) &
        (tabulate(x)[x][seq.int(1,length(x),2)]==2 | tabulate(x)[x][seq.int(2,length(x),2)]==2))


in.first.layer[833]

# Convert sprig assignments x into a design matrix x for clades (n samples x p predictors)
x[is.na(x)] <- p+1L
n <- length(x)/2
x <- Matrix::sparseMatrix(i = rep(1:n,each=2), j = x, x = 1L, dims=c(n,p+1))
x <- x[,-(p+1)]

layer <- list(design2layer(x,which(in.first.layer)),
              design2layer(x,which(!in.first.layer)))


y <- matrix(rnorm(nrow(x)),ncol=1)

test <- rdistill::rdistill(y,x,layers=layer)



###############################################################
# Breadcrumb ideas below
###############################################################

#
# tapply(initial.sprigs.1,initial.sprigs.2,unique)
#
#
# while(anyNA(first.layer)){
#
#   # if there are any edges to follow (or from root, follow them)
#   # otherwise pick highest degree remaining node.
#
#
#
# }
#
#
# #
#
#
# # start line
# ####################
#
# rm(list=ls())
# #saveRDS(sprigs,"~/Downloads/temp_sprigs.rds")
# sprigs <- readRDS("~/Downloads/temp_sprigs.rds")
# p <- sprigs$num.sprigs
# x <- sprigs$assignments
#
# n <- length(x)/2
#
# assigned <- now.assigned <- rep(FALSE,p)
#
# initial.sprigs.1 <- x[seq.int(1,length(x),2)]
# initial.sprigs.2 <- x[seq.int(2,length(x),2)]
#
#
# sprigs.1 <- initial.sprigs.1
# sprigs.2 <- initial.sprigs.2
#
# done <- FALSE
#
#
# # FIXME maybe add a function here to remove any redundant phenotypes (eg: if there exists AB AB genotypes for two doubletons A and B)
#
#
# # Loop over to make layers
# #   Loop up allele counts (don't need to reconsider doubletons once we've checked them b/c still blocked)
# #     Assign all free and unambigous samples recording assignment vector and samples that are still free
# #        note that the assignment vector a is a subset of non-free samples because samples become non-free
# #        even if they just carry an allele that has had some sample assigned to it
#
# assigned.sprig <- rep(FALSE,p)
# ac <- tabulate(x,p)
# max.ac <- max(ac)
#
# layers <- list()
# num.layers <- 0L
# done <- FALSE
# old.num.assigned.sprig <- -1
# num.assigned.sprig <- 0
#
# iterate.over.afs <- TRUE
#
# assign_and_update_exp <- expression({
#   # Assign all unambiguously associated samples
#   assign.sprigs.1 <- which(!assigned.sample & !is.na(sprigs.1) & is.na(sprigs.2))
#   if(length(assign.sprigs.1)){a[assign.sprigs.1] <- sprigs.1[assign.sprigs.1]}
#   assign.sprigs.2 <- which(!assigned.sample & is.na(sprigs.1) & !is.na(sprigs.2))
#   if(length(assign.sprigs.1)){a[assign.sprigs.2] <- sprigs.2[assign.sprigs.2]}
#
#   layer.assigned.sprig <- tabulate(a,p) > 0 # here we only want to track assigned sprigs in this layer?
#   as1 <- layer.assigned.sprig[initial.sprigs.1]
#   as2 <- layer.assigned.sprig[initial.sprigs.2]
#   assigned.sample <-  (as1 & !is.na(as1)) | (as2 & !is.na(as2))
#   table(assigned.sample,useNA="always")
# })
#
#
# while(old.num.assigned.sprig < num.assigned.sprig & num.assigned.sprig < p){
#
#   old.num.assigned.sprig <- num.assigned.sprig
#
#   layer.assigned.sprig <- rep(FALSE,p)
#
#   num.layers <- num.layers + 1
#   assigned.sample <- rep(FALSE,n) # a superset of is.na(a) b/c includes carriers of clades that have had a sample hard assigned to them
#   a <- rep(NA_integer_,n)
#
#   # Remove from layer all previously assigned sprigs
#   layer.sprigs.1 <- initial.sprigs.1
#   layer.sprigs.1[assigned.sprig[initial.sprigs.1]] <- NA_integer_
#   layer.sprigs.2 <- initial.sprigs.2
#   layer.sprigs.2[assigned.sprig[initial.sprigs.2]] <- NA_integer_
#
#
#   if(iterate.over.afs){
#     for(cac in 2:max.ac){
#
#       # remove entries for samples that have already been assigned in this layer or sprigs that do not have the current allele count
#       sprigs.1 <- layer.sprigs.1
#       sprigs.2 <- layer.sprigs.2
#       sprigs.not.in.cacs <- ac!=cac
#       sprigs.1[assigned.sample | sprigs.not.in.cacs[sprigs.1]] <- NA_integer_
#       sprigs.2[assigned.sample | sprigs.not.in.cacs[sprigs.2]] <- NA_integer_
#
#       eval(assign_and_update_exp)
#       #print(paste(sum(layer.assigned.sprig,na.rm = T),sum(assigned.sample,na.rm = T)))
#     }
#   } else {
#     # remove entries for samples that have already been assigned in this layer or sprigs that do not have the current allele count
#     sprigs.1 <- layer.sprigs.1
#     sprigs.2 <- layer.sprigs.2
#
#     sprigs.1[assigned.sample] <- NA_integer_
#     sprigs.2[assigned.sample] <- NA_integer_
#
#     eval(assign_and_update_exp)
#     #print(paste(sum(layer.assigned.sprig,na.rm = T),sum(assigned.sample,na.rm = T)))
#   }
#
#   # samples that could still be assigned to a sprig that's already in the layer
#   to.assign <- which(is.na(a) & (layer.assigned.sprig[initial.sprigs.1] | layer.assigned.sprig[initial.sprigs.2]))
#
#   if(anyNA(initial.sprigs.1[to.assign]) | anyNA(initial.sprigs.2[to.assign])){
#     stop("there are some samples that should have been assigned to a sprig immediately because they are carriers of a sprig and a NA (no sprig) in to.assign")}
#
#   cac <- tabulate(a,p)
#
#   if(length(to.assign)){
#     for(i in 1:length(to.assign)){
#       A <- initial.sprigs.1[to.assign[i]]
#       B <- initial.sprigs.2[to.assign[i]]
#
#       if((cac[A] / ac[A]) > (cac[B] / ac[B])){
#         a[to.assign[i]] <- A
#         cac[A] <- cac[A] + 1L
#       } else {
#         a[to.assign[i]] <- B
#         cac[B] <- cac[B] + 1L
#       }
#     }
#   }
#
#   # confirm new to.assign has length zero
#   if(any(is.na(a) & (layer.assigned.sprig[initial.sprigs.1] | layer.assigned.sprig[initial.sprigs.2]),na.rm = T)){
#     stop("some to.assign were not assigned and should have been")}
#
#   layers[[num.layers]] <- list("a" = a, "n.sprigs"= sum(tabulate(a,p) > 0))
#
#   assigned.sprig <- assigned.sprig | layer.assigned.sprig
#
#   num.assigned.sprig <- sum(assigned.sprig)
#   print(num.assigned.sprig)
# }
#
# # Need to do a bunch of checks for this and figure out why we're getting these really small trailing layers
# # Looks like setting iterate.over.afs to TRUE increases size of early layers, but both versions have about same number of overall layers
#
#
# # interesting case:
#
# cbind(initial.sprigs.1,initial.sprigs.2)[initial.sprigs.1 %in% c(91,552,559) | initial.sprigs.2 %in% c(91,552,559),]
# # where the first entry 17 ends up in
# head(to.assign)
#
# cbind(sprigs.1,sprigs.2)[c(17,339,451),]
#
#
#
#
#
# res1 <- c(length(to.assign),
#           sum(assigned.sprig),
#           sum(is.na(a)),
#           sum(assigned.sample,na.rm=T))
# res1
# #
#
#
#
#
#
#
#
# layers <- list()
# num.layers <- 0L
#
# done <- FALSE
#
# layer.sprigs.1 <- initial.sprigs.1
# layer.sprigs.2 <- initial.sprigs.2
#
#
# # reset sprigs
# sprigs.1 <- layer.sprigs.1
# sprigs.2 <- layer.sprigs.2
#
#
#
#
# # I don't think we need to worry about allele counts -- using fact that within this clade epoch,
# # each sample can only be connecting two clades.  So we have a type of network where clades are the nodes
# # and samples form the edges.  Samples associated with only one clade are free hanging edges (without a second node)
# # All we're doing is assigning free edges to their one connected node, then removing those nodes from the graph and iterating.
# # It cleans up well.
#
# # If we can't get a clade to slot in here b/c it's prop variance explained is not sufficient to out rank it's neighbors with high probability
# # , we just don't test it here in this layer and don't prune it and save for the next epoch of clades.  We also I guess don't test and higher
# # clade including that clade yet?
#
#
#
# # do sprig assignment for unambiguously associated samples
# old.a <- integer()
#
# done <- FALSE
#
# while(!done){
#
#   if(exists("a")){
#     sprigs.1[sprigs.1 %in% unique(a)] <- NA_integer_
#     sprigs.2[sprigs.2 %in% unique(a)] <- NA_integer_
#     a <- pmax(a,sprigs.1,sprigs.2,na.rm = T)
#   } else {
#     a <- pmax(sprigs.1,sprigs.2,na.rm = T)}
#   dual.samples.idx <- which(sprigs.1 != sprigs.2)
#   a[dual.samples.idx] <- NA_integer_
#
#   print(length(na.omit(unique(a))))
#   if(length(na.omit(unique(a))) == length(na.omit(unique(old.a)))){done <- TRUE}
#   old.a <- a
# }
#
# # Then maybe finish up by going back and assigning more samples to any sprigs with the lowest prop.var.
# ac <- tabulate(x,p)
# cac <- tabulate(a)
#
# # samples that could still be assigned to someone
# to.assign <- which(is.na(a) & (!is.na(initial.sprigs.1) | !is.na(initial.sprigs.1)))
#
# if(length(to.assign)){
#
#   for(i in 1:length(to.assign)){
#     A <- initial.sprigs.1[to.assign[i]]
#     B <- initial.sprigs.2[to.assign[i]]
#
#     if((cac[A] / ac[A]) > (cac[B] / ac[B])){
#       a[to.assign[i]] <- A
#       cac[A] <- cac[A] + 1L
#     } else {
#       a[to.assign[i]] <- B
#       cac[B] <- cac[B] + 1L
#     }
#   }
# }
#
# # confirm new to.assign has length zero
# to.assign <- which(is.na(a) & (!is.na(initial.sprigs.1) | !is.na(initial.sprigs.1)))
# length(to.assign)
# #
#
# min(cac/ac)
#
# cac[4]
# cbind(initial.sprigs.1,initial.sprigs.2)[which(initial.sprigs.1==4 | initial.sprigs.2==4),]
# #
#
# initial.sprigs.1
#
#
#
#
#
# sum(table(table(x)))
#
#
# # do sprig assignment for unambiguously associated samples
# if(exists("a")){
#   a <- pmax(a,sprigs.1,sprigs.2,na.rm = T)
# } else {
#   a <- pmax(sprigs.1,sprigs.2,na.rm = T)}
# dual.samples.idx <- which(sprigs.1 != sprigs.2)
# a[dual.samples.idx] <- NA_integer_
# length(unique(a))
#
#
# sprigs.1[sprigs.1 %in% unique(a)] <- NA_integer_
# sprigs.2[sprigs.2 %in% unique(a)] <- NA_integer_
#
# # do sprig assignment for unambiguously associated samples
# if(exists("a")){
#   a <- pmax(a,sprigs.1,sprigs.2,na.rm = T)
# } else {
#   a <- pmax(sprigs.1,sprigs.2,na.rm = T)}
# dual.samples.idx <- which(sprigs.1 != sprigs.2)
# a[dual.samples.idx] <- NA_integer_
# length(unique(a))
#
#
#
#
#
# # for now we just require that at least one sample is assigned to a predictor
# layer.now.assigned <- tabulate(a,p) > 0
#
#
#
#
#
# #### Start line of new section
#
# # FIXME need to add function here to remove any redundant phenotypes (eg: if there exists AB AB genotypes for two doubletons A and B)
# free.sample <- rep(TRUE,n)
#
# ac <- table(x)
# max.ac <- max(ac)
#
# layers <- list()
# num.layers <- 0L
#
# done <- FALSE
#
#
# layer.sprigs.1 <- initial.sprigs.1
# layer.sprigs.2 <- initial.sprigs.2
#
#
#
#
#
#
# while(!done){
#
#   num.layers <- num.layers + 1
#
#   # update assigned and remove already assigned clades/sprigs/predictors
#   assigned[now.assigned] <- TRUE
#   layer.sprigs.1[layer.sprigs.1 %in% which(assigned) ] <- NA_integer_
#   layer.sprigs.2[layer.sprigs.2 %in% which(assigned) ] <- NA_integer_
#
#   layer.assigned <- layer.now.assigned <- rep(FALSE,p)
#
#   for(i in 2:max.ac){
#     #i <- 2L
#     cacs <- 2:i # current allele counts
#
#     # reset sprigs
#     sprigs.1 <- layer.sprigs.1
#     sprigs.2 <- layer.sprigs.2
#
#     # remove already assigned clades/sprigs/predictors AND any samples carrying
#     # this is the key difference between this algorithm and the other.
#     layer.assigned[layer.now.assigned] <- TRUE
#     # FIXME we need a different layer.assigned in the line below because this removing samples
#     # that should not be removed b/c it's also removing samples removed in other layers.
#     to.remove <- sprigs.1 %in% which(layer.assigned) | sprigs.2 %in% which(layer.assigned)
#     sprigs.1[to.remove] <- NA_integer_
#     sprigs.2[to.remove] <- NA_integer_
#
#     # remove sprigs not in current allele counts
#     sprigs.not.in.cacs <- names(ac)[!(ac %in% cacs)]
#     sprigs.1[sprigs.1 %in% sprigs.not.in.cacs] <- NA_integer_
#     sprigs.2[sprigs.2 %in% sprigs.not.in.cacs] <- NA_integer_
#
#     # do sprig assignment for unambiguously associated samples
#     if(exists("a")){
#       a <- pmax(a,sprigs.1,sprigs.2,na.rm = T)
#     } else {
#       a <- pmax(sprigs.1,sprigs.2,na.rm = T)}
#     dual.samples.idx <- which(sprigs.1 != sprigs.2)
#     a[dual.samples.idx] <- NA_integer_
#
#     # for now we just require that at least one sample is assigned to a predictor
#     layer.now.assigned <- tabulate(a,p) > 0
#   }
#
#   layers[[num.layers]] <- list("a" = a, "num.clades" = sum(layer.now.assigned))
#
#   now.assigned <- assigned | layer.now.assigned
#
#   if(!any(now.assigned & !assigned)){done <- TRUE}
#   print(sum(now.assigned))
# }
#
# # Hmm...this is not working either...
#
# sapply(layers,function(x){length(x$num.clades)})
#
# # I think the last layer here is one layer, where as the other
# # algorithm pumps out a different layer with each iteration.
# # I think this is definitely a larger layer...check, but then iterate this
#
#
# #
#
#
#
#
# # Iterate assigning unambiguously associated samples
# while(!done){
#   assigned[now.assigned] <- TRUE
#   sprigs.1[sprigs.1 %in% which(assigned)] <- NA_integer_
#   sprigs.2[sprigs.2 %in% which(assigned)] <- NA_integer_
#
#   if(exists("a")){
#     a <- pmax(a,sprigs.1,sprigs.2,na.rm = T)
#   } else {
#     a <- pmax(sprigs.1,sprigs.2,na.rm = T)}
#   dual.samples.idx <- which(sprigs.1 != sprigs.2)
#   a[dual.samples.idx] <- NA_integer_
#
#   # for now we just require that at least one sample is assigned to a predictor
#   now.assigned <- tabulate(a,p) > 0
#
#   if(!any(now.assigned & !assigned)){done <- TRUE}
#   print(sum(now.assigned))
# }
#
#
#
#
#
#
# # Wait, we only need to run sprig calling to know which neighborhoods to extend and to get higher level clades.  So why don't we run blobby and get
# # all of our aprigs at all of the epochs, and put all of those sprigs into one big design matrix
#
# # then we take that design matrix and run Renyi distillation on it with classic partitioning etc.  Any straggler layers (below a minimum size) we form into
# # a chi square test statistic and weight appropriately.  Additionally, we calculate the QForm out of the rest of the clade matrix, estimate it's number of
# # clades and then add that to this mix (weighted appropriately).  Then we just do the weighted Renyi with two p-values with weight > 1.
#
# # I need to revisit current partitioning algorithm because it seems that it would run into this potential problem where some of the samples I'd like
# # to draw on for a selected predictor have already been assigned elsewhere for calculating W', so many predictors would not be able to obtain 100%
# # prop.var even during extraction.  We could get around this by modifying partitioning by saying that after a predictor is allocated, no other samples  else in it's image
# # can be hard allocated...I need to think this through more...and start by revisiting allocation.  It seems to me here that even simple allocation could
# # easily leave predictors without full prop.var being tested.
#
#
#
# # Hmmm...maybe this isn't the answer b/c pulling in more sprigs in this way means that many sprigs will not be eligible to
# # gather their full sample size later in RD
# # Maybe we do the initial assignment like this but actually break the problem into multiple layers...?
#
# # Well, we could actually allocate some dually-assigned (ambiguously associated) samples, but then we'd just have to sacrifice them in that
# # those samples could not be used in the calculation of the final Ws.  We'd have to renormalize each column of X  so that some samples could
# # be tie-breakers and some could not....but this is almost a different form of Renyi distillation where we have a min
#
# # Maybe for now we stick with Renyi distillation we know and use the standard partitioning method we already have to build a set of distillation layers.
# # once the distillation layers are smaller than some threshold, we don't include any more predictors/clades in that clade round.  Then we go up to the
# # next neighborhoods for all of the allocated sprigs....then we run partitioning again yielding more layers, etc....
#
# # When we partition, prioritize clades from prior epochs, then rare vars, then samples that are not in high demand?
#
# # this way we end up with a bunch of renyi layers with some layers mixing clades from across different epochs of clade calling.
#
#
# # any sprigs that do not have any samples assigned to them at this point must be in competition with other samples
#
# # now let's go through and just assign samples even if there are conflicts.
# # we allocate more samples to sprigs starting with the rarest to the most common variants (breaking ties by prioritizing variants with sprigs with the least var explained)
# # in this first roun, we just go on up allocating, trying to get up to half of variance explained.
#
# # Then we loop back and assign samples
#
#
#
# cbind(initial.sprigs.1,initial.sprigs.2)[is.na(a),]
#
#
# # check that every assignment made in a is valid (corresponds to one of the sprigs carried by that sample)
# sum(!is.na(a)) == length(unique(sort(c(which(a==initial.sprigs.1),
#                                        which(a==initial.sprigs.2)))))
#
#
#
#
#
#
#
#
#
# # temp.free.samples <- rep(FALSE,n)
# # temp.free.samples[dual.samples.idx] <- TRUE
# #
# # # NA assignments for samples that have already been assigned
# # sprigs.1[!temp.free.samples] <- NA_integer_
# # sprigs.2[!temp.free.samples] <- NA_integer_

