set.seed(27)

colpal <- c(palette.colors(10,"Tableau10",alpha = 0.75),"#000000","#E41A1C")

library(locater)
library(kalis)

# Load simulated haplotypes and recombination map
CacheHaplotypes("~/Downloads/msprime_haps_0.h5",transpose=T)
raw.map <- read.table("~/Downloads/map_0.txt")
N()
L()
nthreads <- 8L
# CacheHaplotypes("/home/docker/msprime_haps_0.h5",transpose=T)
# raw.map <- read.table("/home/docker/map_0.txt")
pos <- raw.map[,1]
map <- raw.map[,2]

ploidy <- 2L
n <- N() / ploidy

# running this with high recombination penalty --> lots of big distance jumps
pars <- Parameters(CalcRho(diff(map),10000),mu = 10^-4,use.speidel = FALSE,check.rho = TRUE)

fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)
t <- floor(L()/2)
Forward(fwd,pars,t,nthreads = nthreads)
Backward(bck,pars,t,nthreads=nthreads)


M <- matrix(0,N(),N())
system.time({DistMat(fwd,bck,type = "minus.min",M,nthreads = nthreads)})
class(M) <- "matrix"
system.time({M <- as.dist(Matrix::symmpart(M))})

METHODS <- c("single", "complete", "average", "mcquitty",
             "ward.D", "centroid", "median", "ward.D2")


yy <- matrix(0,length(METHODS),6)
rownames(yy) <- METHODS



for(i in 1:length(METHODS)){
  system.time({d <- .Call(fastcluster:::fastcluster, N(), i, M, NULL)}) # 2L = complete liknage

  # we divide by 2 here because symmpart above divides by 2
  d$height <- d$height/(-log(pars$pars$mu)/2)

  # class(d) <- c("hclust","list")
  # dend <- as.dendrogram(d)
  # cairo_pdf("~/Desktop/large_toy_hap_tree_zoomed_1k.pdf",width = 8,height = 5)
  # plot(dend,xlim=c(1,1000))
  # dev.off()
  #


  n <- length(d$height) + 1L
  d_height <- d$height
  d_height[n] <- 0
  d_merge <- d$merge
  d_merge[d_merge < 0] <- n
  d_merge <- rbind(d_merge,n)
  jumps1 <- d_height - d_height[d_merge[,1]]
  jumps2 <- d_height - d_height[d_merge[,2]]
  yy[i,1] <- mean(jumps1>=1)
  yy[i,2] <- mean(jumps2>=1)
  yy[i,3] <- mean(jumps1>=1 & jumps2 >= 1)
  yy[i,4] <- mean(jumps1>=2)
  yy[i,5] <- mean(jumps2>=2)
  yy[i,6] <- mean(jumps1>=2 & jumps2 >= 2)
  print(i)
}

yy
#yy_4 <- yy
# look at centroid, median and single linkage.



plot(ecdf(c(d_height - d_height[d_merge[,1]],d_height - d_height[d_merge[,2]])))
abline(v=1)






system.time({d <- .Call(fastcluster:::fastcluster, N(), 3L, M, NULL)}) # 2L = complete liknage
d$height <- d$height/(-log(pars$pars$mu)/2)

n <- length(d$height) + 1L
d_height <- d$height
d_height[n] <- 0
d_merge <- d$merge
d_merge[d_merge < 0] <- n
d_merge <- rbind(d_merge,n)
jumps1 <- d_height - d_height[d_merge[,1]]
jumps2 <- d_height - d_height[d_merge[,2]]

# jumps 2 <= jumps 1
# require both to get large jumps

cutoff <- 1
selected_splits <- which(jumps1 >= cutoff | jumps2 >= cutoff)

cl1 <- d_merge[selected_splits,1]
if(any(jumps1[selected_splits] < cutoff)){
  cl1[jumps1[selected_splits] < cutoff] <- n
}
#cl1 <- cl1[cl1!=n]
cl2 <- d_merge[selected_splits,2]
if(any(jumps2[selected_splits] < cutoff)){
  cl2[jumps2[selected_splits] < cutoff] <- n
}
#cl2 <- cl2[cl2!=n]

extract_carriers <- function(x,i){
  if(i == length(d$order)){return(integer())} # i == n indicates i is a leaf
  i1 <- i2 <- i
  while(i1 > 0){i1 <- x$merge[i1,1]}
  while(i2 > 0){i2 <- x$merge[i2,2]}
  x$order[match(-i1,x$order):match(-i2,x$order)]
}

make_design_matrix <- function(x,cl){

  i_list <- j_list <- as.list(1:length(cl))

  for(k in 1:length(cl)){
    i_list[[k]] <- extract_carriers(x,cl[k])
    j_list[[k]] <- rep(k,length(i_list[[k]]))
  }

  Matrix::sparseMatrix(unlist(i_list),unlist(j_list),dims = c(n,length(cl)))
}

X1 <- make_design_matrix(d,cl1)
X2 <- make_design_matrix(d,cl2)

sc1 <- colSums(X1)
sc2 <- colSums(X2)

found1 <- as.list(rep(NA_integer_,L()))
found2 <- as.list(rep(NA_integer_,L()))
found12 <- as.list(rep(NA_integer_,L()))
for(i in 1:L()){
  carrier_haps <- c(QueryCache(loci.idx = i))==1
  ac <- sum(carrier_haps)

  clades_with_matching_ac_1 <- which(sc1 == ac)
  found1[[i]] <- clades_with_matching_ac_1[colSums(X1[carrier_haps,clades_with_matching_ac_1,drop = FALSE])==ac]

  clades_with_matching_ac_2 <- which(sc2 == ac)
  found2[[i]] <- clades_with_matching_ac_2[colSums(X2[carrier_haps,clades_with_matching_ac_2,drop = FALSE])==ac]

  found12[[i]] <- c(found1[[i]],length(cl1) + found2[[i]])

  print(i)
}

found1 <- unlist(found1)
found2 <- unlist(found2)
found12 <- unlist(found12)
# how many variants matched each called clade
tc1 <- tabulate(found1,length(cl1))
plot(ecdf(tc1))
# percentage of called clades that were matched by an observed variant
sum(tc1[cl1!=n]>0)

sum(cl1!=n)

tc2 <- tabulate(found2,length(cl2))
plot(ecdf(tc2))
# percentage of called clades that were matched by an observed variant
sum(tc2[cl2!=n]>0)


sum(tc1[cl1!=n]>0) + sum(tc2[cl2!=n]>0)
sum(cl1!=n) + sum(cl2!=n)



tc12 <- tabulate(found12,length(cl1) + length(cl2))
tc12 <- tc12[1:length(cl1)] + tc12[length(cl1)+1:length(cl2)]
plot(ecdf(tc12))

# percentage of significant jumps where at least one of two immediate children corresponded to an observed variant
mean(tc12[!(cl1==n & cl2 == n)]>0)

sum(tc12[!(cl1==n & cl2 == n)]>0)

#

mean(X[,ncol(X)])

sum(found>0)


#




# Look at concordance between edges called long by single linkage vs. mean or centroid linkage!

# This solution so far counts many positions where one of the children is a leaf, so we're making a big jump when
# we incorporate the loan wolf.  Try this again where we eliminate those...may want to try some other approaches here as well
# around the clustering or symmetrization


# Find mergers of interest


# Find the two branches under each merge of interest and their respective leaves

# I think we can find the children using the order field.  Find the left most descendent of one branch
# and the left most descendent of the next branch to get the contiguous set of children

# for each merger, use distances to figure out which branch carries the mutation
# and select the leaves each of those selected branches to get the predictors to test


# Do the above and try clique calling.


# run RD.

dim(d$merge)
length(d$height)

fc <- fastcluster::hclust(M,method="average")
fc$height <- fc$height /-log(pars$pars$mu)
plot(fc)

hist(fc$height[fc$merge[,1]>0] - fc$height[fc$merge[fc$merge[,1]>0,1]])
hist(fc$height[fc$merge[,2]>0] - fc$height[fc$merge[fc$merge[,2]>0,2]])
plot(fc$height)
which(diff(fc$height/-log(pars$pars$mu))>0.01)
str(fc)


str(M)
M <- matrix(0,N()/2,N()/2)

neigh <- CladeMat(fwd,bck,M,unit.dist = -log(pars$pars$mu),thresh = 0.2,max1var = FALSE)
mean(neigh[[2]][[2]] == neigh[[2]][[3]])
