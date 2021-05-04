#' @include class_DTRSurvStep.R
setClass(Class = "DTRSurv",
         slots = c("call" = "call",
                   "stageResults" = "list",
                   "value" = "ANY",
                   "params" = "Parameters"))

#-------------------------------------------------------------------------------
# method to print key results to screen
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
#' Print Analysis Results
#'
#' Prints the key results of the analysis.
#'
#' @param x A DTRSurv object. The value returned by dtrSurv().
#'
#' @param ... Ignored. 
#'
#' @export
#' @name print
#' @aliases print,DTRSurv-method
#' @returns No return value, called to display key results.
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' # responses must be zero after event
#' evt <- dt[,"D.1"] == 1L
#' dt[evt, "Y.2"] <- 0.0
#'
#' result <- dtrSurv(data = dt, 
#'                   txName = c("A.1", "A.2"), 
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1, 
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' print(x = result)
setMethod(f = "print",
          signature = c(x = "DTRSurv"),
          definition = function(x, ...) {

              cat("\nCall:\n")
              print(x = x@call)

              cat("\n")

              for (i in 1L:length(x@stageResults)) {
                cat("decision point ", i, "\n",
                    "  tx ", x@stageResults[[ i ]]@txName, "\n", 
                    "  tx options: ", x@stageResults[[ i ]]@txLevels, "\n")
              }

              cat("\n")

              if (is(object = x@params, class2 = "Parameters_Mean")) {
                cat("Criterion: Truncated Mean Survival Time\n")
              } else {
                cat("Criterion: Survival Probability at T=", 
                    x@params@survivalTime, 
                    " surv.", x@params@type, "\n")
              }

              cat("\n")

              cat("Estimated Value: ")
              for (i in 1L:length(x = x@value)) {
                cat(" ", names(x@value)[i], ": ", 
                    round(x = x@value[[ i ]], digits = 4), "  ")
              }
              cat("\n")

            })

#-------------------------------------------------------------------------------
# method to show key results to screen
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
#' Show Analysis Results
#'
#' Shows the key results of the analysis.
#'
#' @param object A DTRSurv object. The value returned by dtrSurv().
#'
#' @export
#' @name show
#' @aliases show,DTRSurv-method
#' @returns No return value, called to display key results.
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' # responses must be zero after event
#' evt <- dt[,"D.1"] == 1L
#' dt[evt, "Y.2"] <- 0.0
#'
#' result <- dtrSurv(data = dt, 
#'                   txName = c("A.1", "A.2"), 
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1, 
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' show(object = result)
setMethod(f = "show",
          signature = c(object = "DTRSurv"),
          definition = function(object) {

              cat("\nCall:\n")
              print(x = object@call)

              cat("\n")

              for (i in 1L:length(object@stageResults)) {
                cat("decision point ", i, "\n",
                    "  tx ", object@stageResults[[ i ]]@txName, "\n", 
                    "  tx options: ", object@stageResults[[ i ]]@txLevels, "\n")
              }

              cat("\n")

              if (is(object = object@params, class2 = "Parameters_Mean")) {
                cat("Criterion: Truncated Mean Survival Time\n")
              } else {
                cat("Criterion: Survival Probability at T=", 
                    object@params@survivalTime, 
                    " surv.", object@params@type, "\n")
              }

              cat("\n")

              cat("Estimated Value: ")
              for (i in 1L:length(x = object@value)) {
                cat(" ", names(object@value)[i], ": ", 
                    round(x = object@value[[ i ]], digits = 4), "  ")
              }
              cat("\n")

            })

#-------------------------------------------------------------------------------
# method to return key stage results as a list
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
#' Retrieve Stage Results as a List
#'
#' Returns the key results from all stages or one stage of the Q-learning algorithm.
#'
#' @param object A DTRSurv object. The value returned by dtrSurv().
#'
#' @param ... Ignored. Used to require named inputs.
#'
#' @export
#' @name stage
#' @rdname stage
setGeneric(name = "stage",
           def = function(object, ...) { standardGeneric("stage") })

#' Retrieve Stage Results as a List
#'
#' @rdname dtrSurv-internal-api
#'
setMethod(f = "stage",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

#-------------------------------------------------------------------------------
# method to return key stage results as a list
#-------------------------------------------------------------------------------
# method is exported
#-------------------------------------------------------------------------------
#' Retrieve Stage Results for Decision Point q as a List
#'
#' @param object A DTRSurv object. The value returned by dtrSurv().
#'
#' @param ... Ignored. Used to require named inputs.
#'
#' @param q An integer object. (optional) The stage for which results are 
#'   desired. If q is not provided, all stage results will be returned.
#'
#' @return A list object containing the key results for each requested stage.
#'   If q is not provided, a list of these results will be returned, where the
#'   ith element of that list corresponds to the ith decision point.
#'
#' @rdname stage
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' # responses must be zero after event
#' evt <- dt[,"D.1"] == 1L
#' dt[evt, "Y.2"] <- 0.0
#'
#' result <- dtrSurv(data = dt, 
#'                   txName = c("A.1", "A.2"), 
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1, 
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' tt <- stage(object = result)
setMethod(f = "stage",
          signature = c(object = "DTRSurv"),
          definition = function(object, ..., q) {
              if (missing(x = q)) {
                res <- list()
                for (i in 1L:length(object@stageResults)) {
                  res[[ i ]] <- .stage(object = object@stageResults[[ i ]])
                }
                return( res )
              }
              if (q > length(x = object@stageResults)) {
                stop("q > # of decision points of analysis", call. = FALSE)
              }
              return( .stage(object = object@stageResults[[ q ]]) )
            })


#' Prediction Method
#'
#' Method to estimate the value for new data or to retrieve estimated value for
#'  training data
#'
#' @param object A DTRSurv object. The object returned by a call to dtrSurv().
#'
#' @param ... Ignored. Used to require named inputs.
#'
#' @param newdata NULL or a data.frame object. If NULL, this method retrieves
#'   the estimated value for the training data. If a data.frame, the
#'   value is estimated based on the data provided.
#'
#' @param stage An integer object. The stage for which predictions are desired.
#'
#' @param findOptimal A logical object. If TRUE, the value is estimated for
#'   all treatment options and that leading to the maximum value for each
#'   individual is used to estimate the value.
#'
#' @export
#' @name predict
#' @aliases predict,DTRSurv-method
#' @returns a list object containing a matrix of the predicted survival function
#'   (survFunc), the estimated mean survuval (mean), and the estimated
#'   survival probability (if critical value is surv.mean or surv.prob)
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(101:200,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' # responses must be zero after event
#' evt <- dt[,"D.1"] == 1L
#' dt[evt, "Y.2"] <- 0.0
#'
#' result <- dtrSurv(data = dt, 
#'                   txName = c("A.1", "A.2"), 
#'                   models = list(Surv(Y.1,D.1)~X.1+A.1, 
#'                                 Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' tt <- predict(object = result)
#' tt <- predict(object = result, stage = 1)
#' tt <- predict(object = result, findOptimal = FALSE)
#' tt <- predict(object = result, newdata = dt)
#' tt <- predict(object = result, newdata = dt, stage = 1)
#' tt <- predict(object = result, newdaata = dt, findOptimal = FALSE)

setMethod(f = "predict",
          signature = c(object = "DTRSurv"),
          definition = function(object, 
                                ..., 
                                newdata, 
                                stage = 1,
                                findOptimal = TRUE) {

              if (stage > length(x = object@stageResults)) {
                stop("requested stage not present in analysis", call. = FALSE)
              }

              if (missing(x = newdata)) {
                return( .Predict(object = object@stageResults[[ stage ]],
                                 newdata = NULL,
                                 params = object@params,
                                 findOptimal = findOptimal) )
              } else {

                return( .Predict(object = object@stageResults[[ stage ]],
                                 newdata = newdata,
                                 params = object@params,
                                 findOptimal = findOptimal) )
             }

          })

