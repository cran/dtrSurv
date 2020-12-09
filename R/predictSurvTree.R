# Internal function to predict value of a tree
#
# Function is not exported
#
# @param x A data.frame object {n x p}. The model.frame for covariates to be
#   considered in splitting
#
# @param params A Parameters object. All information that regulates tree and 
#   specifies analysis preferences.
#
# @param nCat An integer object. Number of categories in each covariate. If 
#   covariate i is continuous, nCat[i] = 0; if covariate i is an ordered factor,
#   nCat[i] = 1; if covariate i is an unordered factor, nCat[i] = levels()
#
# @param nodes A list object. The nodes of a tree.
#
#' @include class_Parameters.R
#
.predictSurvTree <- function(..., 
                             x, 
                             params,
                             nCat,
                             nodes) {

  if (is.null(x = x)) return( NULL )

  x <- data.matrix(frame = x)

  usedNodes <- !is.na(x = nodes$nodes[,1L])
  nUsed <- sum(usedNodes)

  n <- nrow(x = x)
  np <- ncol(x = x)
  nTimes <- .NTimes(params)

  nodes$nodes[is.na(nodes$nodes)] <- 0.0

  res <- .Fortran("predictSurvTree",
                  n = as.integer(n),
                  np = as.integer(np),
                  xt = as.double(x = x),
                  nCat = as.integer(x = nCat),
                  nt = as.integer(x = nrow(x = nodes$survFunc)),
                  nNodes = as.integer(x = nUsed),
                  tsurvFunc = as.double(x = nodes$survFunc[,usedNodes]),
                  mean = as.double(nodes$mean[usedNodes]),
                  survProb = as.double(nodes$survProb[usedNodes]),
                  nCols = as.integer(x = ncol(x = nodes$nodes)),
                  tnodes = as.double(x = nodes$nodes[usedNodes,]),
                  predSurvFunc = as.double(numeric(length = n*nTimes)),
                  predMean = as.double(numeric(length = n)),
                  predSurvProb = as.double(numeric(length = n)),
                  PACKAGE = "dtrSurv")

  isSurv <- .CriticalValueCriterion(params) %in% c("surv.mean", "surv.prob")
  valueObj <- list()
  valueObj[[ "survFunc" ]] <- matrix(data = res$predSurvFunc, nrow = nTimes, ncol = n)
  valueObj[[ "mean" ]] <- res$predMean
  if (isSurv) valueObj[[ "survProb" ]] <- res$predSurvProb

  return( valueObj )

}
