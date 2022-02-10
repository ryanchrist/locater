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

  list("y.resids" = y,
       "Q" = Q,
       "sumsq" = sumsq)
}

TestMarker <- function(x, g){
  if(all(g==g[1])){
    return(list("p.value" = rep(NA_real_,ncol(x$y.resids)),
                "y" = x$y.resids, # iid Gaussians under null
                "Q" = x$Q))}

  g <- g - x$Q %*% crossprod(x$Q,g)
  g <- g / sqrt(sum(g^2))
  Q.local <- cbind(g,x$Q)
  temp <- crossprod(Q.local,x$y.resids)
  Z2 <- c(temp[1,])^2
  y.resids.local <- scale(x$y.resids - Q.local %*% temp, center = FALSE) + c(Q.local %*% crossprod(Q.local,rnorm(n)))

  nu <- n-ncol(x$Q)-1

  list("p.value" = exp(pf(q = Z2*nu / (x$sumsq - Z2) , df1 = 1, df2 = nu, lower.tail = FALSE, log.p = TRUE)),
       "y" = y.resids.local, # iid Gaussians under null
       "Q" = Q.local)
}

TestAllMarkers <- function(y,x,Q){
  fast_smt(y,x,Q)
}

