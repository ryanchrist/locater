
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
                     Matrix::colSums(A[!assigned,neighbors_in_l,drop=FALSE]-D[!assigned,neighbors_in_l,drop=FALSE])- num_open_edges)

      proposed_D_col[neighbors_in_l,l] <- temp_a
      proposed_var_budget[neighbors_in_l,l] <- proposed_var_budget[neighbors_in_l,l] - temp_a

      solution_found[l] <- TRUE
      break
    }

  }
})

find_assignment_option_for_js <- expression({

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
  num_open_edges <- A@x[idx[l_neighbors]] - D[[l]][js,neighbors_in_l] # don't need to take pmin here with var_budget
  # because we've given up on caring if taking the variance makes the neighbor ineligible for layer l by this step

  if(sum(num_open_edges) >= shortfall){
    # if there are enough open edges, Look to pull required hard assigned variance from unassigned nodes.
    num_options <- Matrix::rowSums(var_budget[neighbors_in_l,,drop=FALSE] >= 0)
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
  }

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
sparse_mat_1_col <- Matrix::sparseMatrix(1,1,x=0.0,dims=c(p,1))
proposed_D_call_template <- sparse_mat_1_col

assigned <- rep(FALSE,p)
a <- rep(0L,p)


# Declare adjacency matrix A between clades
# note redundant i,j entries added together
A <- Matrix::sparseMatrix(x1,x2, x = 1,
                          dims=c(p+1,p+1))
A <- 2 * Matrix::symmpart(A)

Matrix::diag(A) <- Matrix::diag(A) + A[,p+1]
A <- A[-(p+1),-(p+1)]
A <- as(A,Class = "generalMatrix")
# quick sanity check:
#  ac <- tabulate(x1,p) + tabulate(x2,p) # allele count of each clade
# all.equal(ac,colSums(A))

ac <- Matrix::colSums(A)
tau <- ceiling(min_pv_assigned * ac)
var_budget <- matrix(floor(ac * (1-min_pv_available)),p,L)   # count of alleles available from each clade

while(any(!assigned)){

  # select js, this is slow but OK for now.
  js <- which(!assigned)[order(Matrix::rowSums(var_budget[!assigned,,drop=FALSE]>=0),ac[!assigned],decreasing = TRUE)[1]]

  #if(js==7){stop("yolo")}

  cost <- rep(0,L)
  solution_found <- rep(FALSE,L)
  proposed_D_col <- proposed_D_call_template
  proposed_var_budget <- var_budget

  # LOOP OVER EXISTING LAYERS
  for(l in 1:L){
    if(var_budget[js,l] < 0){next}
    eval(find_assignment_option_for_js)}

  if(any(solution_found)){
    l <- which(solution_found)[which.min(cost[solution_found])]
  } else {
    # make new layer
    l <- L <- L+1
    proposed_D_call_template <- cbind(proposed_D_call_template,sparse_mat_1_col)
    D[[L]] <- Matrix::sparseMatrix(1,1,x=0.0,dims=c(p,p))

    var_budget <- cbind(var_budget,floor(ac * (1-min_pv_available)))

    proposed_D_col <- cbind(proposed_D_col,sparse_mat_1_col)

    proposed_var_budget <- var_budget
    for(l in L){ # done for not to get around problem that break can break out of the outer while loop
      eval(find_assignment_option_for_js)
    }
  }

  a[js] <- l
  D[[l]][,js] <- proposed_D_col[,l]
  var_budget[,l] <- proposed_var_budget[,l]

  # add in validation function here to make sure new assignment is OK.
  assigned[js] <- TRUE

  print(paste(js,"assigned to layer",l,"out of",L,"layers."))
}

all(assigned)
table(a)


# check solution:

# form design matix X
x[is.na(x)] <- p+1L
n <- length(x)/2
X <- Matrix::sparseMatrix(i = rep(1:n,each=2), j = x, x = 1L, dims=c(n,p+1))
X <- X[,-(p+1)]

passed <- rep(FALSE,L)

layers <- as.list(1:L)

for(l in 1:L){

  b <- rep(NA_integer_,n)
  X1 <- X[,a==l]
  rsX1 <- rowSums(X1)==1
  b[rsX1] <- which(a==l)[mat2triplet(t(X1[rsX1,]))$i]
  range(b,na.rm=TRUE)

  # Looks like we have the rest of the needed variance ready to be assigned in D[[1]]!
  # temp_ans <- tapply(X[cbind(seq_len(n),b)],b,sum)
  # cbind(temp_ans,tau[as.integer(names(temp_ans))])[temp_ans < tau[as.integer(names(temp_ans))],]
  # colSums(D[[1]][,c(132,247,294)])

  temp_D <- D[[l]]
  temp_D[a!=l,] <- 0
  temp_D[,a!=l] <- 0

  tD <- mat2triplet(temp_D)

  if(length(tD$i)){
    for(s in 1:length(tD$i)){

      #temp_row_idx <- which(X[,tD$i[s]]!=0 & X[,tD$j[s]]!=0)

      temp_row_idx <- intersect(X[,tD$i[s],drop=FALSE]@i, X[,tD$j[s],drop=FALSE]@i)+1L

      if(!length(temp_row_idx)){stop("hard assignment encoded in D does not correspond to a valid row of X")}

      if(length(temp_row_idx)==1){
        b[temp_row_idx] <-  tD$j[s]
      } else {
        b[sample(temp_row_idx,tD$x[s])] <- tD$j[s]
      }
    }
  }

  temp_ans <- tapply(X[cbind(seq_len(n),b)],b,sum)
  temp_levels <- as.integer(names(temp_ans))
  passed[l] <- all(temp_ans >= tau[temp_levels]) & all(temp_ans + colSums(X[is.na(b),temp_levels]) >= ceiling(ac[temp_levels]*min_pv_available))

  layers[[l]] <- rdistill:::new_rdlayer(a = b,p = p,ortho = FALSE,a_anyNA = anyNA(b),a_levels = temp_levels)

}

all(passed)



################################
################################
# Below is an attempt to more efficiently
# unpack the assignments learned using the
# graph representation. In other words,
# to unpack the hard assignments
# encoded in each D[[l]]
# more efficiently by leveraging
# the Matrix package to replace
# the for loop.
################################
################################

# this is after x1 and x2 NAs have been assigned to p+1
# here we declare any heterozygous NA samples as homozygous samples
# just for assigning samples


# temp solution:
mean(x1==p+1)
x1 <- ifelse(x1%in%which(a==1),x1,p+1)
mean(x1==p+1)

mean(x2==p+1)
x2 <- ifelse(x2%in%which(a==1),x2,p+1)
mean(x2==p+1)

temp_ind <- x1==p+1 & x2!=p+1
x1[temp_ind] <- x2[temp_ind]
temp_ind <- x1!=p+1 & x2==p+1
x2[temp_ind] <- x1[temp_ind]

temp_ind <- x1 <= x2
# need to fix this so that we're not assigning j where we don't have to...pull this from A below using matrix triplet?
# iX <- pmin(x1,x2,na.rm = TRUE)
# jX <- pmax(x1,x2,na.rm = TRUE)
iX <- ifelse(temp_ind,x1,x2)
jX <- ifelse(temp_ind,x2,x1)
rm(temp_ind)
sX <- order(jX,iX)
iX <- iX[sX]
jX <- jX[sX]

#all(iX<=jX)





A <- Matrix::sparseMatrix(x1,x2, x = 1,
                          dims=c(p+1,p+1))[-(p+1),-(p+1)]
A <- Matrix::symmpart(A)
A <- triu(A)
A <- 2 * A
diag(A) <- diag(A)/2


# check that the number of samples to assign in A match the number of samples that there are to assign
# (exclude samples that have NA,NA alleles)
sum(A) == n-sum(x1==p+1 & x2==p+1)

num_to_assign <- sum(A)

iX <- iX[1:num_to_assign]
jX <- jX[1:num_to_assign]
sX <- sX[1:num_to_assign]


# quick check that A is consistent with iX and jX
tA <- mat2triplet(A)
all.equal(iX,rep(tA$i,times=tA$x))
all.equal(jX,rep(tA$j,times=tA$x))

layers <- replicate(L,list("a"=rep(NA_integer_,p)),simplify = FALSE)

temp_D <- D[[1]]
diag(temp_D) <- diag(A) # assign anything on the diagonal (homozygous carriers and heterozygous-NULL carriers to j)
temp_D <- triu(temp_D)

A@x <- c(0L,cumsum(A@x)[-length(A@x)])+1L

temp_A <- A * as(temp_D!=0,"dMatrix")

temp_idx <- sequence(temp_D@x,temp_A@x)

a_temp <- rep(NA_integer_,n)
a_temp[sX[temp_idx]]  <- jX[temp_idx]

all(is.na(a_temp) || a_temp==x1 || a_temp==x2)
all(tabulate(a_temp,p)[a!=1]==0)




new_a_temp <- a_temp

temp_D <- D[[1]]
temp_D <- t(tril(temp_D))
temp_A <- A * as(temp_D!=0,"dMatrix")
temp_idx <- sequence(temp_D@x,temp_A@x)
new_a_temp[sX[temp_idx]] <- iX[temp_idx]

# end of solution: new_a_temp


all(is.na(new_a_temp) || new_a_temp==x1 || new_a_temp==x2)
all(tabulate(new_a_temp,p)[a!=1]==0)

all(X[cbind(1:n,new_a_temp)[-which(is.na(new_a_temp)),]]>0)

t_ans <- tapply(X[cbind(1:n,new_a_temp)],new_a_temp,sum)
cbind(t_ans,tau[as.integer(names(t_ans))])[t_ans < tau[as.integer(names(t_ans))],]

rowSums(X[a==1,])

# write simple loop to check solution
# try simple toy example



layers[[1]]$a <- new_a_temp

#

# temp_D <- cbind(temp_D,A[,p+1][-(p+1)])
# temp_D <- rbind(temp_D,A[,p+1])
# diag(temp_D) <- diag(A) # this will allow us to automatically assign homozygous carriers downstream
# temp_D[p+1,p+1] <- 0
















design2layer_new <- function(X, col.idx = NULL, D = NULL){

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

  a <- rep(NA_integer_,nrow(X))

  if(!is.null(D)){
    for(j in col.idx){
      D[,j]
    }
  }

  X <- X * (Matrix::rowSums(X) == 1)
  var.assigned <- Matrix::colSums(X^2)


  prop.var <- var.assigned / var.total


  X <- Matrix::mat2triplet(X)
  a[X$i] <- col.idx[X$j]

  list("a" = a,
       "a.levels" = a.levels,
       "prop.var" = prop.var,
       "ortho" = FALSE,
       "a.anyNA"= anyNA(a))
}


l <- 1

# Convert sprig assignments x into a design matrix x for clades (n samples x p predictors)

test <- design2layer(x,which(a==1))

layers <- list(),
design2layer(x,which(!in.first.layer)))




#list(a,D)

