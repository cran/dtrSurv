# Class extends CriticalValue to indicate that critical value is survival
#
# Class is not exported and is for internal convenience only
#
# @slot SurvivalTime A numeric object. The time at which the survival 
#   probability is to be estimated
#
# @slot sIndex An integer object. The index of the timePoint vector above
#   which the survival time lies (and it is below the sIndex + 1 element)
#
# @slot sFraction A numeric object. The fractional location of SurvivalTime in
#   t[sIndex] and t[sIndex+1]
#
# @slot type A character object. Indicates of mean is to be used to break ties.
#
# Methods
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {defined}
#  .IsSurvival(object, ...) {new; defined}
#
# Functions
# .criticalValueSurvival(SurvivalTime, timePoints)
#
#' @include class_CriticalValue.R
setClass(Class = "CriticalValueSurvival",
         slots = c(survivalTime = "ANY",
                   sIndex = "integer",
                   sFraction = "numeric",
                   type = "character"),
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------
setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueSurvival"),
          definition = function(object, ...) { 
              if (object@type == "mean") return( "surv.mean" )
              if (object@type == "prob") return( "surv.prob" )
            })

setMethod(f = "initialize",
         signature = c(.Object = "CriticalValueSurvival"),
         def = function(.Object, ..., survivalTime, sIndex, sFraction, type) {

                   obj <- list(...)
                   tst <- sapply(X = obj, 
                                 FUN = function(x){
                                         is(object = x, 
                                            class2 = "CriticalValueSurvival")
                                        })

                   if (any(tst)) {
                     .Object <- obj[[ which(tst) ]]
                   } else if (missing(x = survivalTime) && missing(x = sIndex) &&
                         missing(x = sFraction) && missing(x = type)) {
                     .Object@survivalTime <- Inf
                     .Object@sIndex <- -1L
                     .Object@sFraction <- 0
                     .Object@type <- "none"
                   } else {
                     if (missing(x = survivalTime) || missing(x = sIndex) ||
                         missing(x = sFraction) || missing(x = type)) {
                         gn <- unlist(lapply(list(...),is))
                         if( "CriticalValueSurvival" %in% gn ) return(.Object)
                         stop("insufficient inputs provided")
                     }
                     if (type %in% c("surv.prob", "prob")) {
                       .Object@type <- "prob"
                     } else if (type %in% c("surv.mean", "mean")) {
                       .Object@type <- "mean"
                     }
                     .Object@survivalTime <- survivalTime
                     .Object@sIndex <- sIndex
                     .Object@sFraction <- sFraction
                   }
                   return( .Object )
                 })

#-------------------------------------------------------------------------------
# method to identify if critical value is of survival type
#-------------------------------------------------------------------------------
# method returns a logical
#-------------------------------------------------------------------------------
setGeneric(name = ".IsSurvival",
           def = function(object, ...) { standardGeneric(".IsSurvival") })

setMethod(f = ".IsSurvival",
          signature = c(object = "ANY"),
          definition = function(object, ...) { return( FALSE ) })

setMethod(f = ".IsSurvival",
          signature = c(object = "CriticalValueSurvival"),
          definition = function(object, ...) { return( TRUE ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueSurvival
#-------------------------------------------------------------------------------
# function returns a CriticalValueSurvival object
#-------------------------------------------------------------------------------
.criticalValueSurvival <- function(survivalTime, timePoints, type) {

  nTimes <- length(x = timePoints)

  # index of last time point <= SurvivalTime
  sIndex <- sum(timePoints <= survivalTime)

  if (sIndex < nTimes) {
    # if SurvivalTime is below tau, determine the fraction
    sFraction <- {survivalTime - timePoints[sIndex]} / 
                 {timePoints[sIndex + 1L] - timePoints[sIndex]}
  } else if (sIndex == 0L) {
    # if it is below the minimum, stop -- cannot extrapolate
    stop("survival time is below minimum timepoint")
  } else {
    # if it is above tau, use tau -- cannot extrapolate
    sFraction <- 0.0
    survivalTime <- max(timePoints)
    message("survival time reset to tau")
  }

  return( new("CriticalValueSurvival",
              "survivalTime" = survivalTime,
              "sIndex" = sIndex,
              "sFraction" = sFraction,
              "type" = type) )

}

