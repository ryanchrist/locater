#' Test a Clade Matrix for association with phenotypes
#'
#' Test a local phylogeny for association with phenotypes based on a matrix prodced by \code{kalis::CladesMat}
#
#' @param y a matrix of phenotypes
#' @param M a matrix as produced by \code{kalis::CladesMat}
#' @param Q an orthogonal matrix whose columns span the column space of the background covariates
#' @param traces an optional traces argument
#' @return
#'   With with p-values for Clade Testing
#' @seealso
#'   \code{\link{Clades}} to generate \code{kalisClades} object
#' @examples
#' \dontrun{
#' }
#'
#' @export
TestCladeMat <- function(y,M,Q,traces = NULL, other.test.pvalues = NULL, point.est = FALSE, k = NULL, use.forking = FALSE, nthreads = 1L){

  # process other.test.pvalues dropping any tests that have any missing values
  if(!is.null(other.test.pvalues)){
    if(is.numeric(other.test.pvalues)){
      if(any(is.na(other.test.pvalues))){
        other.test.res <- NULL
      } else {
        other.test.res <- as.list(-log10(other.test.pvalues))
      }
    } else if(is.list(other.test.pvalues)) {
      other.test.pvalues <- other.test.pvalues[sapply(other.test.pvalues,function(x){!any(is.na(x))})]
      if(length(other.test.pvalues)){
        other.test.res <- lapply(other.test.pvalues,function(x){-log10(x)})
        other.test.res <- data.table::transpose(other.test.pvalues)
      } else {
        other.test.res <- NULL
      }
    } else {
      stop("other.test.pvalues must be a vector or a list")
    }
  }

  if(!is.null(k)){k <- as.integer(k)}
  n <- nrow(y)

  traces <- calc_traces(M, Q = Q, nthreads = nthreads)

  if(!traces$hsnorm2){
    warning("PMP had all zero entries -- perhaps no clades were called.")
    if(boughs.branches){return(rep(NA_real_,2))}else{return(NA_real_)}
  }

  matmul <- function(x,args){
    x <- x - Q %*% crossprod(Q,x)
    x <- M %*% x
    x - Q %*% crossprod(Q,x)
  }

  if(use.forking){
    parallel.sapply <- function(x,FUN,...){unlist(parallel::mclapply(x,FUN,...,mc.cores = nthreads))}
  } else {
    parallel.sapply <- sapply
  }

  ###     Test Overall Quadratic Form         ###
  ###############################################

  obs <- c(colSums(y * matmul(y,1)))

  CalcBounds(f = function(k,args){RSpectra::eigs_sym(matmul,
                                                     k = k, n = n, args = args,
                                                     opts = list("ncv" = min(n, max( 4*((2*k+1)%/%4+1), 20)) ))},
             # ensures that we'll always be using a multiple of 4 unless we request all of the eigenvalues.
             args = 0,
             obs = obs,
             traces,
             k = if(!is.null(k)){k}else{if(point.est){min(n-1,200)}else{unique(pmin(n-1,c(4,20,200)))}}, # need n-1 instead of n here because if n,
             # RSpectra falls back to eigen which cannot handle functional definition of matrix
             only.point.est = if(point.est){TRUE}else{FALSE}, # if TRUE, then the first element of k is used for calculated trunc part
             only.bounds.est = if(!point.est){TRUE}else{FALSE},
             lower.tail = FALSE, #if(j <= num.dist.mats){TRUE}else{FALSE}
             other.test.res = other.test.res,
             parallel.sapply = parallel.sapply)
}
