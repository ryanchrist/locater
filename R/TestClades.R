#' Test a Clade Matrix (local relatedness matrix) for association with phenotypes
#'
#' Test a local relatedness matrix for association with phenotype(s).
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
#' @param y a \code{n} by \code{m} matrix of quantitative phenotypes for \code{n} samples and \code{m} phenotypes
#' @param M a local \code{n} by \code{n} relatedness/clade matrix, such as one produced by \code{kalis::CladesMat}
#' @param Q a \code{n} by \code{q} orthogonal matrix whose columns span the column space of the background covariates
#' @param k a vector of non-negative integers, number of eigenvalues to calculate to try to approximate the null distribution, see Details.
#' @param min.prop.var a double in [0,1], if we've obtained enough eigenvalues to explain at least this minimum proportion of variance, then we can trust our approximation and stop calculating eigenvalues, see Details.
#' @param var.ratio.goal a double in [0,1], if lambda_k / lambda_{k-1} >= var.ratio.goal, the spectrum has plateaued, then we can trust our approximation and stop calculating eigenvalues, see Details.
#' @param stop.eval.func a function that returns TRUE if all associations can be ruled insignificant and eigendecomposition stopped early, see Details.
#' @param calc.obs.T experimental feature coming soon that will allow return of truncated quadratic form statistics \code{T}
#' @param use.forking a logical, should forking be used during underlying FFT?
#' @param nthreads a positive integer, how many cores in multithreaded routines?
#' @return A data.frame with \code{m} rows (the first row corresponding to the first phenotype/column in \code{y}). The data frame has fields:
#' \itemize{
#' \item \code{qform}: the quadratic form -log10 p-value
#' \item \code{obs}: the observed value of the quadratic form test statistic
#' \item \code{obs.T}: the observed value of the truncated quadratic form test statistic (experimental feature, coming soon, all NAs for now)
#' \item \code{precise}: a logical, if TRUE the value of \code{qform} can be trusted based on our robust approximation of the remainder term -- the eigendecomposition stopping criteria was reached.
#' If FALSE, the largest number of eigenvalues in the \code{k} provided was not sufficient to satisfy the stopping criteria so the approximation may be unreliable -- consider retrying with a larger eigendecomposition budget in \code{k}.
#' \item \code{exit.status}: an integer, if 0 no errors were encountered.
#' }
#' @seealso
#'   \code{\link{Clades}} to generate \code{kalisClades} object
#' @examples
#' \dontrun{
#' }
#'
#' @export
TestCladeMat <- function(y, M, Q,
                         k = NULL, # c(8,32,128,512)
                         prop.var.goal = 0.95,
                         var.ratio.goal = 0.95,
                         stop.eval.func = NULL,
                         calc.obs.T = FALSE,
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
                   prop.var.goal = prop.var.goal,
                   var.ratio.goal = var.ratio.goal,
                   k = k, # vector of positive integers, if null, we just do SW approx taking k=0.
                   stop.eval.func = stop.eval.func,
                   calc.obs.T = calc.obs.T,
                   parallel.sapply = base::sapply)
}
