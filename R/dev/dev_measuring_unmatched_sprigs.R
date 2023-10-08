set.seed(27)

library(locater)
library(kalis)
verbose <- TRUE
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
pars <- Parameters(CalcRho(diff(map),10^-4),mu = 10^-4,use.speidel = FALSE,check.rho = TRUE)

fwd <- MakeForwardTable(pars)
bck <- MakeBackwardTable(pars)
t <- floor(L()/2)
Forward(fwd,pars,t,nthreads = nthreads)
Backward(bck,pars,t,nthreads=nthreads)

M <- matrix(0,N()/2,N()/2)
neigh <- CladeMat(fwd,bck,M,unit.dist = -log(pars$pars$mu),thresh = 0.5, max1var = TRUE, nthreads = nthreads)

sprigs <- Sprigs(neigh[[1]])

sc <- tabulate(sprigs$assignments,nbins = sprigs$num.sprigs)

require(Matrix)
haps <- Matrix(QueryCache(),sparse = TRUE)
haps <- t(haps)
ac <- colSums(haps)
hist(ac)
a <- sprigs$assignments
a[is.na(a)] <- 0L

num_sprigs <- sprigs$num.sprigs
#
# a <- rep(0L,N())
# a[haps[,1]==1] <- 1L
# num_sprigs <- 1L


found <- rep(0L,num_sprigs)
for(i in 1:num_sprigs){
  carriers <- which(a == i)
  candidate_vars <- which(ac==length(carriers))
  for(j in 1:length(candidate_vars)){
    if(identical(carriers,haps[,candidate_vars[j],drop=FALSE]@i+1L)){
      found[i] <- candidate_vars[j]
      break
    }
  }
  print(i)
}

mean(found>0)

# Allele count distribution abount variants that are found
table(ac[found[found>0]])

# Allele count distribution in general
colpal <- c(palette.colors(10,"Tableau10",alpha = 0.75),"#000000","#E41A1C")

barplot(rbind(table(ac[ac<=6]),
              table(ac[found[found>0]]),
              table(sc[!found])[1:5]),las=1,beside = T,col = colpal[1:3])
legend("topright",fill = colpal[1:3],
       legend=c("obs variants","matched sprigs","unmatched sprigs"))


M <- matrix(0,N(),N())
system.time({DistMat(fwd,bck,type = "minus.min",M,nthreads = nthreads)})
class(M) <- "matrix"
system.time({M <- as.dist(Matrix::symmpart(M))})
system.time({d <- .Call(fastcluster:::fastcluster, N(), 2L, M, NULL)}) # 2L = complete liknage
#

dd <- diff(d$height)/(-log(pars$pars$mu))

# Find mergers of interest
which(dd>=0.5)

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
