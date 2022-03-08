FitNull <- function(y, A = NULL){

  if(!is.matrix(y)){y <- as.matrix(y)}
  if(is.null(A)){A = matrix(1,nrow(y),1)}

  if(nrow(y) != nrow(A)){
    stop("nrow(y) must equal nrow(A)")
  }

  if(!any(apply(A,2,function(x){all(x == x[1])},simplify = T))){
    A <- cbind(1,A) # add intercept to A if no intercept in input
  }

  Q <- qr.Q(qr(A))
  y.resids <- y - Q %*% crossprod(Q,y)
  sumsq <- colSums(y.resids^2)

  list("y" = y,
       "Q" = Q,
       "sumsq" = sumsq)
}

TestMarker <- function(x, g){
  if(all(g==g[1])){
    return(list("p.value" = rep(NA_real_,ncol(x$y)),
                "y" = x$y, # iid Gaussians under null
                "Q" = x$Q))}

  n <- nrow(x$Q)
  g <- g - x$Q %*% crossprod(x$Q,g)
  g <- g / sqrt(sum(g^2))
  Q.local <- cbind(g,x$Q)
  temp <- crossprod(Q.local,x$y)
  Z2 <- c(temp[1,])^2
  y.resids.local <- scale(x$y - Q.local %*% temp, center = FALSE) + c(Q.local %*% crossprod(Q.local,rnorm(n)))

  # use statmod::qresiduals (randomized quantile residuals to extend this to GLMs)

  nu <- n-ncol(x$Q)-1
  list("p.value" = exp(pf(q = Z2*nu / (x$sumsq - Z2) , df1 = 1, df2 = nu, lower.tail = FALSE, log.p = TRUE)),
       "y" = y.resids.local, # iid Gaussians under null
       "Q" = Q.local)
}

#' Test All Markers
#' Fast method for performing SMT on a genotype matrix
#'
#' @export
TestCachedMarkers <- function(y, A = NULL, from = 1L, to = kalis::L(), ploidy = 2L, model = "additive"){

  if(is.null(kalis::N())){stop("haplotypes must be cached with kalis::CacheHaplotypes before running FindTargetVars")}
  if(length(from)!=1 || from > kalis::L() || from < 1){
    stop("from must be an integer in [1,kalis::L()]")}
  if(length(to)!=1 || to > kalis::L() || to < 1){
    stop("to must be an integer in [1,kalis::L()]")}
  if(from>to){
    stop("from must be <= to")}

  if(!is.matrix(y)){y <- as.matrix(y)}
  if(is.null(A)){A = matrix(1,nrow(y),1)}

  if(nrow(y) != nrow(A)){
    stop("nrow(y) must equal nrow(A)")
  }

  X <- Haps2Genotypes(QueryCache(loci.idx = from:to),
                      ploidy = ploidy,
                      method = model)
  res <- fast_smt(y,X,qr.Q(qr(A)))
  if(!is.null(colnames(y))){colnames(res) <- colnames(y)}
  rownames(res) <- from:to

  as.data.frame(res)
}

