
# demonstration of how our projection-based phenotype simulation approach
# stabalizes the effects induced for causal variants


require(Matrix)
require(rdistill)
N.haps <- rep(2e4,3)
n <- sum(N.haps)/2
target_effect_per_variant <- 4
A <- cbind(1,as(t(fac2sparse(rep(1:3,times=N.haps/2))),"denseMatrix"),rbinom(n,1,0.5),rnorm(n))
rsA <- rowSums(A)
qA <- qr(A)

n.reps <- 500

x1 <- x2 <- x3 <- matrix(0,n.reps,2)

sign_vec <- c(-1,-1)

for(i in 1:n.reps){

  G <- cbind(runif(n) < 0.01,runif(n) < 0.01)
  storage.mode(G) <- "integer"
  G <- as(G,"sparseMatrix")
  table(G[,1],G[,2])
  X <- std_sparse_matrix(G)

  qX_causal <- qr(X)
  Q <- qr.Q(qX_causal)

  causal_idx <- c()
  causal_r2 <- c()

  for(j in sample(ncol(X))){
    r2 <- sum(qr.fitted(qA,X[,j])^2)
    causal_r2  <- c(causal_r2,r2)
  }

  qAX <- qr(cbind(A,X))

  zz <- rnorm(n)

  y1 <- rank2gauss(rank(
    qr.resid(qA, zz + as.vector(X %*% (sign_vec * qnorm(10^(-target_effect_per_variant)/2))))
    , ties.method = "random"))

  x1[i,] <- -pnorm(tail(qr.coef(qAX,y = y1),2),lower.tail = FALSE, log.p = TRUE)/log(10)

  y2 <- rank2gauss(rank(
    qr.resid(qA, qr.resid(qX_causal,zz) + as.vector(X %*% (sign_vec * qnorm(10^(-target_effect_per_variant)/2)/sqrt(1-causal_r2))))
    , ties.method = "random"))

  x2[i,] <- -pnorm(tail(qr.coef(qAX,y = y2),2),lower.tail = FALSE, log.p = TRUE)/log(10)

  y3 <- rank2gauss(rank(
    qr.resid(qA, qr.resid(qX_causal,zz) + as.vector(X %*% (sign_vec * qnorm(10^(-target_effect_per_variant)/2))))
    , ties.method = "random"))

  x3[i,] <- -pnorm(tail(qr.coef(qAX,y = y3),2),lower.tail = FALSE, log.p = TRUE)/log(10)


  print(i)
}


layout(matrix(1:4,2,2))
plot(x1,ylim = c(0,12),xlim=c(0,12))
abline(h=4)
abline(v=4)
plot(x2,ylim = c(0,12),xlim=c(0,12))
abline(h=4)
abline(v=4)

plot(rowMeans(x1),abs(x1[,1] - x1[,2]))
plot(rowMeans(x2),abs(x2[,1] - x2[,2]))
plot(rowMeans(x3),abs(x3[,1] - x3[,2]))

plot(x2[,1],x3[,1])
abline(0,1)

plot(x2[,2]/x3[,2])
abline(0,1)

