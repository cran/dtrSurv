# Class to store information regarding time points
#
# Class is not exported and is for internal convenience only
#
#  @slot timePoints A numeric vector; the timepoints used in the analysis
#
#  @slot timeDiff A numeric vector; the time differences between the timepoints
#
#  @slot tau A numeric object; maximum time
#
# Getters
#   .TimePoints(object, ...) {new; defined}
#   .NTimes(object, ...) {new; defined}
#   .Tau(object, ...) {new; defined}
#   .TimeDiff(object, ...) {new; defined}
#
# Methods
#   .TimeInfoAsList(object, ...) {new; defined}
#
# Functions
#   .timeInfo(timePoints, nTimes, response)
#
setClass(Class = "TimeInfo",
         slots = c(timePoints = "numeric",
                   timeDiff = "numeric",
                   tau = "numeric"))

### Getters

#-------------------------------------------------------------------------------
# Method to retrieve timepoints
#-------------------------------------------------------------------------------
# Method returns the vector of timepoints used for the analysis
#-------------------------------------------------------------------------------
setGeneric(name = ".TimePoints",
           def = function(object, ...) { standardGeneric(".TimePoints") })

setMethod(f = ".TimePoints",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TimePoints",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { return( object@timePoints ) })

#-------------------------------------------------------------------------------
# Method to retrieve the number of timepoints
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------
setGeneric(name = ".NTimes",
           def = function(object, ...) { standardGeneric(".NTimes") })

setMethod(f = ".NTimes",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NTimes",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { 
              return( length(x = object@timePoints) ) 
            })

#-------------------------------------------------------------------------------
# Method to retrieve the maximum timepoint
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------
setGeneric(name = ".Tau",
           def = function(object, ...) { standardGeneric(".Tau") })

setMethod(f = ".Tau",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".Tau",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { return( object@tau ) })

#-------------------------------------------------------------------------------
# Method to retrieve the difference in timepoints
#-------------------------------------------------------------------------------
# Method returns the vector of the differences in timepoints used for the analysis
#-------------------------------------------------------------------------------
setGeneric(name = ".TimeDiff",
           def = function(object, ...) { standardGeneric(".TimeDiff") })

setMethod(f = ".TimeDiff",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TimeDiff",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { return( object@timeDiff ) })


#-------------------------------------------------------------------------------
# Method to retrieve timepoints for printing
#-------------------------------------------------------------------------------
# Method returns a list containing 1 element "timePoints" containing the 
#  vector of timepoints used for the analysis
#-------------------------------------------------------------------------------
setGeneric(name = ".TimeInfoAsList",
           def = function(object, ...) { 
                   standardGeneric(".TimeInfoAsList") 
                 })

setMethod(f = ".TimeInfoAsList",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TimeInfoAsList",
          signature = c(object = "TimeInfo"),
          definition = function(object, ...) { 
              return( list("timePoints" = object@timePoints) )
            })

#-------------------------------------------------------------------------------
# Function to verify inputs and create a TimeInfo object
#-------------------------------------------------------------------------------
# Function returns a TimeInfo object
#-------------------------------------------------------------------------------
#' @include VerifyTimePoints.R
.timeInfo <- function(timePoints, nTimes, tau, response) {

  # ensure that timePoints and nTimes are appropriate. Methods return a vector
  # of unique time points that are sorted in ascending order.
  timePoints <- .VerifyTimePoints(timePoints = timePoints,
                                  tau = tau, 
                                  nTimes = nTimes,
                                  response = response)

  tau <- timePoints$tau
  timePoints <- timePoints$timePoints

  # the total number of times points
  nTimes <- length(x = timePoints)

  # deltaT; T_{i} - T_{i-1} i = 1:nTimes with T_0 = 0
  # timeDiff <- c(timePoints[-1L] - timePoints[-nTimes], 0)
  # this is used for the truncated mean calculation in which it is
  # assumed that all t > tau equals tau to t(nTimes+1)-t(nTimes) = 0

  timeDiff <- timePoints[-1L] - timePoints[-nTimes]
  timeDiff[nTimes] <- 0.0

  return( new(Class = "TimeInfo",
              "tau" = tau,
              "timePoints" = timePoints,
              "timeDiff" = timeDiff) )

}
