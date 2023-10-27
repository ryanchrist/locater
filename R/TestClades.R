#' Test a Clade Matrix for association with phenotypes
#'
#' Test a local phylogeny for association with phenotypes based on a matrix prodced by \code{kalis::CladesMat}
#'
#' Concerning the specification of k, k=0 represents evaluating no eigenvalues and simply calculating the Satterthwaite (SW) approximation based on the traces of the matrix.
#' If k is NULL, k is set to 0. If the k provided does not include 0, 0 is appended to k so the SW approximation is always calculated.
#' stop.eval.function is a function that can be used to prematurely stop the procedure if the results are almost guaranteed to be insignificant. This function is called
#' after each number of k eigenvalues is calculated (including after k=0).
#' It always takes one or two arguments. The first argument must always be a vector of estimated QForm p-values in [0,1]. If a second argument is given in the function definition,
#' the proportion of variance explained by the current set of evaluated eigenvalues (a number in [0,1]) will be passed to it. This allows the user to use different early stopping rules
#' based on the p-values of other test statistics, the current estimated set of QForm p-values, and the proportion of variance captured by the current set of eigenvalues.
#' This allows the user to do an initial screen with just the SW approximation (leaving k and stop.eval.func as NULL). If we always want to bipass the SW approximation and
#' try more accurate approximations, we can define stop.eval.func to return FALSE whenever the second argument (prop.var) is 0.
#'
#' @param y a matrix of phenotypes
#' @param M a matrix as produced by \code{kalis::CladesMat}
#' @param Q an orthogonal matrix whose columns span the column space of the background covariates
#' @param k a vector of non-negative integers, number of eigenvalues to calculate to try to approximate the null distribution, see Details.
#' @param min.prop.var a double in [0,1], if we've obtained enough eigenvalues to explain at least this minimum proportion of variance, then we can trust our approximation and stop calculating eigenvalues, see Details.
#' @param var.ratio.goal a double in [0,1], if lambda_k / lambda_{k-1} >= var.ratio.goal, the spectrum has plateaued, then we can trust our approximation and stop calculating eigenvalues, see Details.
#' @param stop.eval.func a function that returns TRUE if all associations can be ruled insignificant and eigendecomposition stopped early, see Details.
#' @param cs.approx a logical, should a difference of chi-squares be used to try to approximate remainder or should all of the remainder be left to a Gaussian approximation (the default)
#' @param use.forking a logical, should forking be used during underlying FFT?
#' @param nthreads a positive integer, how many cores in multithreaded routines?
#' @return
#'   With with p-values for Clade Testing
#' @seealso
#'   \code{\link{Clades}} to generate \code{kalisClades} object
#' @examples
#' \dontrun{
#' }
#'
#' @export
TestCladeMat <- function(y, M, Q,
                         k = NULL, # c(8,32,128,512)
                         min.prop.var = 0.98,
                         var.ratio.goal = 0.9,
                         stop.eval.func = NULL,
                         cs.approx = FALSE,
                         use.forking = FALSE,
                         nthreads = 1L){

  if(!is.null(k)){k <- as.integer(k)}
  n <- nrow(y)

  traces <- calc_traces(M, Q = Q, nthreads = nthreads)

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

  SimpleCalcBounds(y,
                   matmul,
                   traces,
                   min.prop.var = 0.98,
                   var.ratio.goal = 0.98,
                   k = k, # vector of positive integers, if null, we just do SW approx taking k=0.
                   stop.eval.func = stop.eval.func,
                   cs.approx = cs.approx,
                   parallel.sapply = base::sapply)
}
