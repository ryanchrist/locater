divide_interval <- function(x,min.gap){
  # x be map along interval including end points which are evaluated
  if(length(x)==2){return(integer(0))} # no points in interval to add

  idx <- rep(NA_integer_,ceiling((tail(x,1) - x[1]) / min.gap))
  j <- 0

  lead.x <- x[1]

  for(i in seq.int(2,length(x)-1)){
    if(x[i] - lead.x >= min.gap){
      j <- j+1
      idx[j] <- i
      lead.x <- x[i]
    }
  }
  na.omit(idx)
}



# map <- cumsum(runif(1e3))
# targets0 <- sort(sample.int(length(map),20,replace = FALSE))
#
# FillGaps(targets0,map,min.gap=1)


#' Helper for Finding Target Variants for Screening
#'
#' Use recombination distance (and preliminary single marker test results) to select target variants for screening
#
#' @param map a vector specifying the recombination map 'CDF' at all cached variants
#' @param min.cM the desired minimum recombination distance between target loci in cM
#' @param from the first variant index to consider as a target
#' @param to the last variant index to consider as a target
#' @param initial.targets an integer vector or data.frame.
#' If an integer vector, the indices of variants that should be included in the returned targets.
#' If a \code{data.frame}, single marker test results in the format returned by \code{\link{TestCachedMarkers}} to be used
#' in combination with \code{smt.thresh} to implicitly define initial target variants
#' @param smt.thresh when \code{initial.targets} is provided as a \code{data.frame} of single marker test results,
#' a -log10 p-value threshold such that variants with any signal above that threshold should be included in initial targets
#' @return
#'   integer vector of target variant indices
#' @seealso
#'   \code{\link{TestCachedMarkers}} to select variants based on initial single marker test results
#' @examples
#' \dontrun{
#' }
#'
#' @export
FindTargetVars <- function(map, min.cM = 0.1, from = 1L, to = length(map),
                           initial.targets = NULL, # have the precalculated by previous function?
                           smt.thresh = NULL){

  if(!is.numeric(map) | length(map) == 0){stop("map must be a numeric vector")}

  if(length(map)==1){return(1L)}
  if(length(map)==2){return(1:2)}

  if(!is.null(min.cM) && (!is.finite(min.cM) || length(min.cM) !=1 || min.cM < 0)){
    stop("min.cM must be a non-negative number")}

  if(length(from)!=1 || from > length(map) || from < 1){
    stop("from must be an integer in [1,length(map)]")}
  if(length(to)!=1 || to > length(map) || to < 1){
    stop("to must be an integer in [1,length(map)]")}
  if(from>to){
    stop("from must be <= to")}
  if(from==to){
    return(as.integer(from))}
  if(from + 1L == to){
    return(as.integer(c(from,to)))}

  if(is.null(initial.targets)){
    initial.targets <- integer()
  } else if(is.integer(initial.targets)){
    # move on to next step
  } else if(is.data.frame(initial.targets)){
    if(is.null(smt.thresh)){
      stop("smt.thresh must be provided if initial.targets is defined implicitly with a data.frame of SMT results, see Details")}
    if(!is.finite(smt.thresh) || length(smt.thresh) !=1 || smt.thresh < 0){
      stop("smt.thresh must be a non-negative number")}

    initial.targets <- as.integer(rownames(initial.targets)[which(c(rowSums(initial.targets > smt.thresh))>0)])
  } else {
    stop("initial.targets must be NULL, an integer vector, or defined implicitly with a data.frame of SMT results, see Details")
  }

  if(any(initial.targets <= 0) | any(initial.targets > length(map))){
    stop("all initial.targets must be in [1,length(map)]")}

  idx <- initial.targets

  idx <- idx[ idx >= from & idx <= to]
  if(!length(idx)){idx <- c(from,to)}

  long.intervals <- which(diff(map[idx]) >= min.cM)
  if(!length(long.intervals)){return(idx)} # no gaps to fill

  l <- lapply(long.intervals,FUN = function(x){idx[x] - 1L + divide_interval(map[idx[x]:idx[x+1]],min.cM)})
  sort(Reduce(union,c(l,list(idx))))
}
