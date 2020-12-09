# Class for storing primary results from a single stage of the Q-learning
#   survival analysis
#
# Class is not exported and is for internal convenience only
#
#  @slot txName A character object. The name of the treatment variable for the
#    Q-learning step
#
#  @slot txLevels A vector. The treatment options available.
#
#  @slot model A formula object. The model for extracting the covariates to be
#    considered for splitting
#
#  @slot survRF A SurvRF object. The primary results of the tree building
#    algorithm
#
#  @slot eligibility A logical vector object. TRUE indicates that the case was
#    eligible to be included in analysis
#
#  @slot valueAllTx A list object. The value of the tree for each tx level.
#
#  @slot optimal An Optimal object. The estimated optimal tx and optimal value.

# Methods
#
# Functions
#
#' @include class_Optimal.R class_SurvRF.R
#' @include class_Parameters.R

setClass(Class = "DTRSurvStep",
         slots = c("txName" = "character",
                   "txLevels" = "vector",
                   "model" = "formula",
                   "survRF" = "ANY",
                   "eligibility" = "logical",
                   "valueAllTx" = "list",
                   "optimal" = "Optimal"))

#-------------------------------------------------------------------------------
# function to return key stage results as a list
#-------------------------------------------------------------------------------
# function is not exported
#-------------------------------------------------------------------------------
.stage <- function(object, ...) {
  result <- list()

  result[[ "txName" ]] <- object@txName
  result[[ "txLevels" ]] <- object@txLevels
  result[[ "sampleSize" ]] <- sum(object@eligibility)
  result[[ "valueAllTx" ]] <- object@valueAllTx
  result[[ "optimal" ]] <- .OptimalAsList(object = object)
  result[[ "stages" ]] <- .stageSurvRF(object = object@survRF)

  return( result )
}


.meanValue <- function(object, ...) {
  res <- list()

  # returns the mean of the expected survival times
  res[[ "Et" ]] <-  mean(x = apply(X = object@valueAllTx$mean, 
                                   MARGIN = 1L,  
                                   FUN = max))
  if (!is.null(object@valueAllTx$survProb)) {
    res[[ "St" ]] <-  mean(x = apply(X = object@valueAllTx$survProb, 
                                     MARGIN = 1L,  
                                     FUN = max))
  }
  return( res )
}

# generic defined in class_Optimal.R
setMethod(f = ".OptimalAsList",
          signature = c(object = "DTRSurvStep"),
          definition = function(object, ...) {
                return( .OptimalAsList(object = object@optimal) )
              })

#-------------------------------------------------------------------------------
# method to retrieve predicted values
#-------------------------------------------------------------------------------
# if findOptimal is TRUE, method stops with error
# if findOptimal is FALSE, method returns a Value object
#-------------------------------------------------------------------------------
setMethod(f = ".Predict",
          signature = c(object = "DTRSurvStep",
                        newdata = NULL),
          definition = function(object, newdata, ...) {

              return( .Predict(object = object@survRF, newdata = NULL, ...) )

            })

#-------------------------------------------------------------------------------
# method to predict value for new data
#-------------------------------------------------------------------------------
# if findOptimal is TRUE, method returns a list containing a Value object and 
#   an Optimal object
# if findOptimal is FALSE, method returns a Value object
#-------------------------------------------------------------------------------
#' @include class_Optimal.R
#' @importFrom stats model.frame
setMethod(f = ".Predict",
          signature = c(object = "DTRSurvStep",
                        newdata = "data.frame"),
          definition = function(object, 
                                newdata, 
                                ...,  
                                params,
                                findOptimal) {

              # update model to remove response
              mod <- update(object@model, NULL ~ .)

              # ensure data contains all model covariates
              x <- tryCatch(expr = stats::model.frame(formula = mod, 
                                                      data = newdata),
                            error = function(e) {
                                      stop("variables in the training data missing in newdata",
                                           call. = FALSE)
                                     })

              # remove response from x
              # this should no longer happen, but keeping anyway
              if (attr(x = terms(x = mod), which = "response") == 1L) {
                x <- x[,-1L,drop = FALSE]
              }

              if (findOptimal) {
                # if optimal is requested make predictions for all possible
                # treatment options

                resV <- .PredictAll(object = object@survRF, 
                                    newdata = newdata, 
                                    params = params,
                                    model = mod,
                                    txName = object@txName,
                                    txLevels = object@txLevels)

                return( list("value" = resV$predicted, "optimal" = resV$optimal) )

              } else {
                return( .Predict(object = object@survRF, 
                                 newdata = x, 
                                 params = params, ...) )
              }
          })


# Internal function to create random forest
#
# Function is not exported
#
# @param model A survival formula object, the rhs of which specifies the 
#   covariates to be considered in the splitting algorithm
#
# @params data A data.frame object containing covariate and treatment histories
#
# @params priorStep A DTRSurvStep object. The analysis from a previous step
#
# @params params A Parameters object.
#
# @params txName A character object or a character vector object. The names
#   of the treatment variables in data
#
# @params mTry An integer object or NULL. If integer, the maximum number of 
#    covariates to consider for each split. If NULL, mTry is sqrt(np)
#
# @params sampleSize An integer object. The number of samples to draw for each
#    tree
#
#' @importFrom stats na.pass
#' @importFrom stats update
#' @importFrom stats terms
#' @importFrom stats complete.cases
#' @importFrom stats model.frame
#' @importFrom stats model.response
#' @include shiftMat.R survRF.R
.dtrSurvStep <- function(..., 
                         model, 
                         data, 
                         priorStep, 
                         params,
                         txName,
                         mTry,
                         sampleSize) {

  mod <- model

  # identify order 1 terms in formula
  order1 <- attr(x = stats::terms(x = mod), which = "order") == 1L
  if (any(order1)) {
    stageCov <- attr(x = stats::terms(x = mod), which = "term.labels")[order1]
  } else {
    stop("problem in identifying covariates, verify formula\n", call. = FALSE)
  }

  # warn about order > 1
  orderHigh <- attr(x = stats::terms(x = mod), which = "order") > 1L
  if (any(orderHigh)) message("interaction terms are ignored")

  # extract model frame
  x <- stats::model.frame(formula = mod, data = data, na.action = na.pass)

  message("model ", appendLF = FALSE)
  tm <- as.character(mod)
  message(tm[2], " ~ ", tm[3])

  # identify individuals with complete data
  elig <- stats::complete.cases(x)

  # extract response and delta from model frame
  response <- stats::model.response(data = x)
  delta <- response[,2L]
  response <- response[,1L]

  # remove response from x
  if (attr(x = terms(x = mod), which = "response") == 1L) {
    x <- x[,-1L,drop = FALSE]
  }

  # responses that are zero indicate censored at a previous stage
  zeroed <- abs(x = response) < 1e-8

  elig <- elig & !zeroed

  if (sum(elig) == 0L) stop("no cases have complete data", call. = FALSE)

  message("cases in stage: ", sum(elig))

  # maximum number of covariates to try
  if (is.null(x = mTry)) {
    mTry <- as.integer(x = ceiling(x = sqrt(x = ncol(x = x))))
    message("maximum # of covariates considered for splitting set to ", mTry)
  } else if (mTry > ncol(x = x)) {
    message("mTry reset as it is larger than the # of available covariates")
    mTry <- as.integer(x = ceiling(x = sqrt(x = ncol(x = x))))
    message("maximum # of covariates considered for splitting set to ", mTry)
  } else {
    mTry <- as.integer(x = mTry)
  }

  if (is.null(x = priorStep)) {
    # priorStep is NULL for first step of the analysis.
    # transform the time variable to a probability mass vector

    # identify time points <= response
    tSurv <- sapply(X = response[elig], 
                    FUN = function(s, tp) { as.integer(x = {s < tp}) },
                    tp = .TimePoints(object = params))

    # time point nearest the response without going over
    # {nTimes x nElig}

    pr <- {rbind(tSurv[-1L,],1)-tSurv}
  } else {
    # priorStep is not NULL when q < Q
    #
    # the number of timepoints
    # .NTimes() is a getter method defined for Parameters objects
    nTimes <- .NTimes(object = params)

    # create an empty matrix for all previously eligible cases
    survMatrix <- matrix(data = 0.0, 
                         nrow = nTimes, 
                         ncol = nrow(x = x))

    survMatrix[1L,] <- 1.0

    # retrieve estimated survival function from previous step
    survMatrix[,priorStep@eligibility] <- t(.OptimalY(object = priorStep@optimal))

    # shift the survival function down in time (T_i - Tq) and 
    # transform to a probability mass vector for only those
    # eligible for this stage
    # .shiftMat is an internal function defined in shiftMat.R

    pr <- .shiftMat(timePoints = .TimePoints(object = params), 
                    survMatrix = survMatrix[,elig,drop = FALSE], 
                    shiftVector = response[elig],  
                    surv2prob = TRUE)

    pr[abs(pr)<1e-8] <- 0.0

  }

  if (any(is.na(x = pr))) stop("NA not permitted in pr -- contact maintainer", 
                               call. = FALSE)

  if (any(pr > 1.0) || any(pr < 0.0)) {
    stop("pr must obey 0 <= pr <= 1 -- contact maintainer", call. = FALSE)
  }

  # identify tx levels in limited data
  if (is.factor(x = data[,txName])) {
    txLevels <- levels(x = factor(x = data[elig,txName]))
  } else {
    txLevels <- sort(x = unique(x = data[elig, txName]))
  }

  if (length(x = txLevels) == 1L) {
    message("***only one treatment level in data***")
  }

  if (.Pooled(object = params)) {

      message("pooled analysis; treatments ", paste(txLevels,collapse=" "))

      # this will be a SurvRF object
      result <- .survRF(x = x[elig,,drop=FALSE], 
                        y = response[elig], 
                        pr = pr, 
                        delta = delta[elig], 
                        params = params,
                        mTry = mTry,
                        txLevels = txLevels,
                        model = mod,
                        sampleSize = sampleSize)

  } else {

    message("stratified analysis")

    # result will be a list of SurvRF objects
    result <- list()

    for (i in 1L:length(x = txLevels)) {

      message("  treatment level ", txLevels[i])

      nms <- as.character(x = txLevels[i])

      di <- {data[elig,txName] == txLevels[i]}

      use <- elig & {data[,txName] == txLevels[i]}

      result[[ nms ]] <- .survRF(x = x[use,,drop=FALSE], 
                                 y = response[use], 
                                 pr = pr[,di], 
                                 delta = delta[use], 
                                 params = params,
                                 mTry = mTry,
                                 txLevels = txLevels[i],
                                 model = mod,
                                 sampleSize = sampleSize)

    }

    result <- new(Class = "SurvRFStratified", "strat" = result)

  }

  # calculate the estimated values for all treatment levels
  # .PredictAll() is a method; called here for objects of class SurvRF, which
  # is defined in file class_SurvRF.R
  resV <- .PredictAll(object = result, 
                      newdata = data[elig,], 
                      params = params, 
                      model = mod,
                      txName = txName,
                      txLevels = txLevels)

  result <- new(Class = "DTRSurvStep",
                "txName" = txName,
                "txLevels" = txLevels,
                "model" = mod,
                "survRF" = result,
                "eligibility" = elig,
                "valueAllTx" = resV$predicted,
                "optimal" = resV$optimal)

  return( result )

}
