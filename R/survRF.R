# Internal function to generate forest
#
# Function is not exported
#
# @param x A data.frame object {n x p}. The model.frame for covariates to be
#   considered in splitting
#
# @param delta An integer vector object {n}. Vector of indicators for censoring 
#   (1 = not censored; 0 = censored)
#
# @param pr A matrix object {nt x n}. Probability mass vector of survival 
#   function
#
# @param params A Parameters object. All information that regulates tree and 
#   specifies analysis preferences.
#
# @param mTry An integer object. The maximum number of covariates to use 
#   for splitting algorithm.
#
# @params sampleSize An integer object. The number of samples to draw for each
#    tree
#
#' @include class_Parameters.R
#' @include class_CriticalValue.R class_CriticalValueMean.R
#' @include class_CriticalValueSurvival.R
#' @include class_SurvRF.R
#' @import parallel
.survRF <- function(..., x, delta, pr, params, mTry, sampleSize) {

  # if x_i is an unordered factor, nCat_i is the number of levels 
  # if x_i is not a factor, nCat is 0
  # if x_i is an ordered factor, nCat is 1
  xLevels <- lapply(X = x, FUN = levels)

  nCat <- sapply(X = x, FUN = nlevels)

  nCat <- ifelse(test = sapply(X = x, FUN = is.ordered), yes = 1L, no = nCat)

  # number of individuals in training data
  nSamples <- nrow(x = x)

  # number of time points
  nTimes <- nrow(x = pr)

  # total number of trees to be grown in the forest
  # .NTree() is a getter method defined for Parameters objects
  nTree <- .NTree(object = params)

  # determine the number of samples to include in each tree
  sampleSize <- ceiling(x = sampleSize * nSamples)

  message("sampleSize: ", sampleSize)

  # maximum number of nodes in a tree
  maxNodes <- 2L * sampleSize + 1L

  # convert factors to integers
  x = data.matrix(frame = x)

  nr = nrow(x = x)

  # send step specific x, pr, delta, mTry, nCat to Fortran
  res = .Fortran("setUpInners",
                 t_n = as.integer(x = nrow(x = x)),
                 t_np = as.integer(x = ncol(x = x)),
                 t_x = as.double(x = x),
                 t_pr = as.double(x = t(x = pr)),
                 t_delta = as.integer(x = delta),
                 t_mTry = as.integer(x = mTry),
                 t_nCat = as.integer(x = nCat),
                 t_sampleSize = as.integer(x = sampleSize),
                 t_nTree = as.integer(x = params@nTree),
                 t_nrNodes = as.integer(x = maxNodes),
                 PACKAGE = "dtrSurv")

  survTree <- .Fortran("survTree",
                       forestSurvFunc = as.double(numeric(nTimes*nr)),
                       forestMean = as.double(numeric(nr)),
                       forestSurvProb = as.double(numeric(nr)),
                       PACKAGE = "dtrSurv")

  # retrieve trees from Fortran 
  trees <- list()
  for (i in 1L:nTree) {

    trees[[ i ]] <- list()
    treeDims <- .Fortran("treeDim", 
                         iTree = as.integer(x = i), 
                         nr = as.integer(x = 1L), 
                         nc = as.integer(x = 1L),
                         PACKAGE = "dtrSurv")

    temp <- .Fortran("getTree", 
                     iTree = as.integer(x = i), 
                     nr = as.integer(x = treeDims$nr), 
                     nc = as.integer(x = treeDims$nc), 
                     nodes = as.double(x = numeric(length = treeDims$nc*treeDims$nr)),
                     survFunc = as.double(x = numeric(length = treeDims$nr*nTimes)),
                     mean = as.double(x = numeric(length = treeDims$nr)),
                     survProb = as.double(x = numeric(length = treeDims$nr)),
                     PACKAGE = "dtrSurv")

    trees[[ i ]]$nodes <- matrix(data = temp$nodes, 
                                 nrow = treeDims$nr,  
                                 ncol = treeDims$nc)
    trees[[ i ]]$survFunc <- matrix(data = temp$survFunc, 
                                    nrow = nTimes, 
                                    ncol = treeDims$nr)
    trees[[ i ]]$mean <- temp$mean

    trees[[ i ]]$survProb <- temp$survProb
  }

  forest <- list()
  forest[[ "survFunc" ]] <- matrix(data = survTree$forestSurvFunc,
                                   nrow = nTimes, ncol = nr)
  forest[[ "mean" ]] <- survTree$forestMean

  crit <- .CriticalValueCriterion(params)
  if (crit %in% c("surv.mean", "surv.prob")) {
    forest[[ "survProb" ]] <- survTree$forestSurvProb
  }

  return( new(Class = "SurvRF",
              "trees" = trees,
              "forest" = forest,
              "variables" = colnames(x = x),
              "mTry" = mTry,
              "nCat" = nCat,
              "xLevels" = xLevels) )

}
