# Virtual class to denote objects are arising from survRF step
#
# Methods
#   .Predict(object, newdata, ...) {new; not allowed}
setClass(Class = "SurvRFObject",
         contains = c("VIRTUAL"))

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
setGeneric(name = ".Predict",
           def = function(object, newdata, ...) { standardGeneric(".Predict") })

setMethod(f = ".Predict",
          signature = c(object = "ANY",
                        newdata = "ANY"),
          definition = function(object, newdata, ...) { stop("not allowed") })

#-------------------------------------------------------------------------------
# method to make predictions for new data at each tx level
#-------------------------------------------------------------------------------
setGeneric(name = ".PredictAll",
           def = function(object, ...) { 
                   standardGeneric(".PredictAll") 
                 })

setMethod(f = ".PredictAll",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

# Class for storing survRF results for pooled analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot trees A list object. The results of the tree building algorithm for
#    each tree in the forest
#
#  @slot forest A list object. The values averaged across all trees in the
#    forest
#
#  @slot variable A character vector. The variables considered in the
#    analysis
#
#  @slot mTry An integer. The maximum number of covariates considered for 
#    splitting
#
#  @slot nCat An integer vector. The number of categories for each covariate
#    considered. >=2 unordered factor, 1 ordered factor, 0 continuous
#
#  @slot xLevels A list object. The categories in each covariate considered.
#
# Methods
#   .Predict(object, newdata, ...) {defined}
#
setClass(Class = "SurvRF",
         slots = c("trees" = "list",
                   "forest" = "list",
                   "variables" = "character",
                   "mTry" = "integer",
                   "nCat" = "integer",
                   "xLevels" = "list"),
         contains = c("SurvRFObject"))

.stageSurvRF <- function(object) {
  res <- list()
  if (is(object = object, class2 = "SurvRFStratified")) {
    for (i in 1L:length(x = object@strat)) {
      res[[ i ]] <- list()
      res[[ i ]][[ "trees" ]] <- object@strat[[ i ]]@trees
      res[[ i ]][[ "forest" ]] <- object@strat[[ i ]]@forest
    }
  } else {
    res[[ "trees" ]] <- object@trees
    res[[ "forest" ]] <- object@forest
  }
}

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with NULL retrieves fitted values from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing survFunc, mean, and? survProb
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "SurvRF",
                        newdata = NULL),
          definition = function(object, newdata, ...) {
              return( object@forest )
            })


#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculates value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing survFunc, mean, and? survProb
#-------------------------------------------------------------------------------
#' @include predictSurvTree.R
setMethod(f = ".Predict",
          signature = c(object = "SurvRF",
                        newdata = "data.frame"),
          definition = function(object, 
                                newdata, 
                                ..., 
                                params) {

              # verify that there are not new levels in the the data
              # this assumes that newdata has been passed in with
              # covariates in the order used in the analysis. This is
              # guaranteed if the data.frame is created from the model.
              xLevels <- lapply(X = newdata, FUN = levels)

              for (i in length(x = xLevels)) {
                if (is.null(x = xLevels[[ i ]]) && 
                    is.null(x = object@xLevels[[ i ]])) next
                if (any(! {xLevels[[ i ]] %in% object@xLevels[[ i ]]})) {
                  stop("new factor levels not present in the training data",
                       call. = FALSE)
                }
              }

              # verify that type of data is the same as the training data
              # type means numeric (nCat = 0), ordered factor (nCat = 1), or 
              # factor (nCat = length(levels(x)))
              nCat <- sapply(X = xLevels, FUN = length)

              nCat <- ifelse(test = sapply(X = newdata, FUN = is.ordered), 
                             yes = 1L, 
                             no = nCat)

              if (any(unlist(x = object@nCat) != unlist(x = nCat))) {
                stop("type of predictors in newdata do not match the training data",
                       call. = FALSE)
              }

              nTree <- length(x = object@trees)

              # predict for first tree
              # .predictSurvTree() is an internal function defined in predictSurvTree.R

              newResult <- .predictSurvTree(x = newdata, 
                                            params = params, 
                                            nCat = nCat,
                                            nodes = object@trees[[ 1L ]])

              iTree <- 2L
              while (iTree <= nTree) {

                # predict for tree iTree; sum result
                # .predictSurvTree() is an internal function defined in predictSurvTree.R
    
                tmp <- .predictSurvTree(x = newdata, 
                                        params = params, 
                                        nCat = nCat,
                                        nodes = object@trees[[ iTree ]])

                newResult$survFunc <- newResult$survFunc + tmp$survFunc
                newResult$mean <- newResult$mean + tmp$mean
                if (!is.null(x = newResult$survProb)) {
                  newResult$survProb <- newResult$survProb + tmp$survProb
                }

                iTree <- iTree + 1L
              }

              newResult$survFunc <- newResult$survFunc / nTree
              newResult$mean <- newResult$mean / nTree
              if (!is.null(x = newResult$survProb)) {
                newResult$survProb <- newResult$survProb / nTree
              }

              return( newResult )

            })

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculates value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object containing survFunc, mean, and? survProb
#-------------------------------------------------------------------------------
setMethod(f = ".PredictAll",
          signature = c(object = "SurvRF"),
          definition = function(object, ..., newdata, model, txLevels, txName, params) {

              # set all cases to receive first treatment
              newdata[[ txName ]][] <- txLevels[1L]

              # extract new model frame
              x <- stats::model.frame(formula = model, data = newdata)

              # remove response from x
              if (attr(x = terms(x = model), which = "response") == 1L) {
                x <- x[,-1L,drop = FALSE]
              }

              # calculate the estimated values for this treatment level
              # .Predict() is a method; called here for objects of class SurvRF, which
              # is defined in this file
              res <- .Predict(object = object, newdata = x, params = params, ...)

              res$survFunc <- list(res$survFunc)

              i <- 2L

              # repeat this process for each tx level
              while (i <= length(x = txLevels)) {

                # set all cases to receive the ith treatment
                newdata[[ txName ]][] <- txLevels[i]

                # extract new model frame
                x <- stats::model.frame(formula = model, data = newdata)

                # remove response from x
                if (attr(x = terms(x = model), which = "response") == 1L) {
                  x <- x[,-1L,drop=FALSE]
                }

                # calculate the estimated values for this treatment level
                # .Predict() is a method; called here for objects of class SurvRF, which
                # is defined in file class_SurvRF.R
                tt <-  .Predict(object = object, 
                                newdata = x, 
                                params = params, ...)

                res[[ "survFunc" ]][[ i ]] <- tt$survFunc
                res[[ "mean" ]] <- cbind(res$mean, tt$mean)
                res[[ "survProb" ]] <- cbind(res$survProb, tt$survProb)

                i <- i + 1L
              }

              opt <- .optimal(params = params, predicted = res, txLevels = txLevels)

              return( list("predicted" = res, "optimal" = opt) )
            })

# Class for storing survRF results for stratified analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot strat A list object. The results of the tree building algorithm for
#   each tx group. List will be a list of SurvRF objects
#
# Methods
#   .Predict(object, newdata, ...) {defined}
#
setClass(Class = "SurvRFStratified",
         slots = c("strat" = "list"),
         contains = c("SurvRFObject"))

#-------------------------------------------------------------------------------
# method to make predictions for new data or to retrieve fitted values
#-------------------------------------------------------------------------------
# method with NULL retrieves fitted values from SurvRFStratified object
#-------------------------------------------------------------------------------
# method returns a list object each element a list containing survFunc, mean, 
#   and? survProb
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "SurvRFStratified",
                        newdata = NULL),
          definition = function(object, newdata, ...) {
              res <- list()
              for (i in 1L:length(x = object@strat)) {
                res[[ i ]] <- .Predict(object@strat[[ i ]], newdata = NULL, ...)
              }
              return( res )
            })

#-------------------------------------------------------------------------------
# method makes predictions for new data or retrieves fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculated value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object each element a list containing survFunc, mean, 
#   and? survProb
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "SurvRFStratified",
                        newdata = "data.frame"),
          definition = function(object, newdata, ..., params) {

              res <- list()

              res[[ 1L ]] <- .Predict(object = object@strat[[ 1L ]], 
                                      newdata = newdata,  
                                      params = params, ...)

              i <- 2L
              while (i <= length(x = object@strat)) {
                res[[ i ]] <- .Predict(object = object@strat[[ i ]], 
                                       newdata = newdata, 
                                       params = params, ...)
                i <- i + 1L
              }

              return( res )

            })


#-------------------------------------------------------------------------------
# method makes predictions for new data or retrieves fitted values
#-------------------------------------------------------------------------------
# method with data.frame calculated value for new data from SurvRF object
#-------------------------------------------------------------------------------
# method returns a list object each element a list containing survFunc, mean, 
#   and? survProb
#-------------------------------------------------------------------------------
setMethod(f = ".PredictAll",
          signature = c(object = "SurvRFStratified"),
          definition = function(object, ..., newdata, model, params, txLevels) {

              # extract new model frame
              x <- stats::model.frame(formula = model, data = newdata)

              # remove response from x
              if (attr(x = terms(x = model), which = "response") == 1L) {
                x <- x[,-1L,drop=FALSE]
              }

              res <- .Predict(object = object@strat[[ 1L ]], 
                              newdata = x, 
                              params = params, ...)

              res$survFunc <- list(res$survFunc)

              i <- 2L
              while (i <= length(x = object@strat)) {
                tt <- .Predict(object = object@strat[[ i ]], 
                               newdata = x, 
                               params = params, ...)

                res[[ "survFunc" ]][[ i ]] <- tt$survFunc
                res[[ "mean" ]] <- cbind(res$mean, tt$mean)
                res[[ "survProb" ]] <- cbind(res$survProb, tt$survProb)

                i <- i + 1L
              }

              opt <- .optimal(params = params, predicted = res, txLevels = txLevels)

              return( list("predicted" = res, "optimal" = opt) )

            })

.optimal <- function(params, predicted, txLevels) {

  crit <- .CriticalValueCriterion(params)

  if (crit == "mean") {

    # identify which element contains the maximum expected survival time
    optTx <- apply(X = predicted$mean,
                   MARGIN = 1L,
                   FUN = .whichMax,
                   tieMethod = params@tieMethod)

  } else if (crit == "surv.mean") {
    # identify which element
    # contains the maximum expected survival time
    optTx <- apply(X = predicted$mean,
                   MARGIN = 1L,
                   FUN = .whichMax,
                   tieMethod = params@tieMethod)

    # index for which the survival probability is max
    tmp <- apply(X = predicted$survProb,
                 MARGIN = 1L,
                 FUN = .whichMax,
                 tieMethod = "NA")

    # for those that are not tied, replace mean survival time
    # with mean survival probability
    isna <- is.na(x = tmp)
    optTx[!isna] <- tmp[!isna]

  } else if (crit == "surv.prob") {

    # identify which element contains the maximum expected survival probability
    optTx <- apply(X = predicted$survProb,
                   MARGIN = 1L,
                   FUN = .whichMax,
                   tieMethod = params@tieMethod)
  }

  # extract the survival function at optimal tx
  optSv <- matrix(data = 0.0, nrow = length(x = optTx), ncol = .NTimes(params))

  for (i in 1L:length(optTx)) {
    optSv[i,] <- predicted$survFunc[[optTx[i]]][,i]
  }    

  return( new(Class = "Optimal",
              "optimalTx" = txLevels[optTx], 
              "optimalY" = optSv,
              "type" = crit) )
}

#-------------------------------------------------------------------------------
# internal function to identify the maximum value of input vector x
#  @param x a vector of values 
#  @param tieMethod a character indicating the method to be used to 
#    breaks {first, random, NA}
#-------------------------------------------------------------------------------
# function returns a single numeric or NA
#-------------------------------------------------------------------------------
.whichMax <- function(x, tieMethod) {

  ind <- which(x >= {max(x) - 1e-8})

  if (length(x = ind) == 1L) return( ind )

  if (tieMethod == "first") return( ind[1L] )

  if (tieMethod == "random") return( resample(x = ind, size = 1L) )

  return( NA )
}
