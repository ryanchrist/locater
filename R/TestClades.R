#' Test a Clade Matrix for association with phenotypes
#'
#' Test a local phylogeny for association with phenotypes based on a matrix prodced by \code{kalis::CladesMat}
#
#' @param y a matrix of phenotypes
#' @param M a matrix as produced by \code{kalis::CladesMat}
#' @param Q an orthogonal matrix whose columns span the column space of the background covariates
#' @param traces an optional traces argument
#' @param point.est a logical, if TRUE evaluate eigenvalues according to schedule k until at least 99% of the variation in M is captured and return the resulting p-value.
#' If FALSE (default), stop evaluation of eigenvalues once bounds can rule out that locus is significant
#' @return
#'   With with p-values for Clade Testing
#' @seealso
#'   \code{\link{Clades}} to generate \code{kalisClades} object
#' @examples
#' \dontrun{
#' }
#'
#' @export
TestCladeMat <- function(y, M, Q, traces = NULL,
                         other.test.pvalues = NULL,
                         min.prop.var = 0.98,
                         k = c(10,100,200,400,800),
                         neg.log10.cutoff = NULL,
                         use.bettermc = FALSE, use.forking = FALSE, nthreads = 1L){

  # if any of the other.test.pvalues are missing/non.finite, they are just ignored
  # when they are combined with the QF bounds to assess whether to continue eigendecomposition
  if(!is.null(other.test.pvalues)){
    if(is.numeric(other.test.pvalues)){
        other.test.res <- as.list(-log10(other.test.pvalues))
    } else if(is.list(other.test.pvalues)) {
        other.test.res <- lapply(other.test.pvalues,function(x){-log10(x)})
        other.test.res <- data.table::transpose(other.test.res)
    } else {
      stop("other.test.pvalues must be a vector or a list")
    }
  }

  if(!is.null(k)){k <- as.integer(k)}
  n <- nrow(y)

  traces <- calc_traces(M, Q = Q, nthreads = nthreads)

  matmul <- function(x,args){
    x <- x - Q %*% crossprod(Q,x)
    x <- M %*% x
    x - Q %*% crossprod(Q,x)
  }

  if(use.forking){
    if(use.bettermc){
      parallel.sapply <- function(x,FUN,...){unlist(bettermc::mclapply(x,FUN,...,mc.cores = nthreads))}
    } else {
      parallel.sapply <- function(x,FUN,...){unlist(parallel::mclapply(x,FUN,...,mc.cores = nthreads))}
    }
  } else {
    parallel.sapply <- sapply
  }

  ###     Test Overall Quadratic Form         ###
  ###############################################

  SimpleCalcBounds(y,
                   matmul,
                   traces,
                   min.prop.var = min.prop.var,
                   k = k,
                   neg.log10.cutoff = neg.log10.cutoff,
                   other.test.res = other.test.res,
                   lower.tail = FALSE,
                   parallel.sapply = parallel.sapply)
}
