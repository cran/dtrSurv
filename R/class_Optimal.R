# Class for storing estimated optimal treatment and value
#
# Class is not exported and is for internal convenience only
#
#  @slot optimalTx A vector object. The index of the estimated optimal tx
#
#  @slot optimalY A vector. The estimated value
#
#  @slot type A character. One of "mean" or "prob"
#
# Getters
#  .OptimalY(object, ...) {new; defined}
#  .OptimalAsList(object, ...) {new; defined}
#
setClass(Class = "Optimal",
         slots = c("optimalTx" = "vector",
                   "optimalY" = "matrix",
                   "type" = "character"))

#-------------------------------------------------------------------------------
# Method to retrieve the estimated optimal value
#-------------------------------------------------------------------------------
# Method returns a vector object
#-------------------------------------------------------------------------------
setGeneric(name = ".OptimalY",
           def = function(object, ...) { standardGeneric(".OptimalY") })

setMethod(f = ".OptimalY",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".OptimalY",
          signature = c(object = "Optimal"),
          definition = function(object, ...) { return( object@optimalY ) })

#-------------------------------------------------------------------------------
# Method to return an Optimal object as a list (used for printing)
#-------------------------------------------------------------------------------
# Method returns a list object
#-------------------------------------------------------------------------------
setGeneric(name = ".OptimalAsList",
           def = function(object, ...) { standardGeneric(".OptimalAsList") })

setMethod(f = ".OptimalAsList",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".OptimalAsList",
          signature = c(object = "Optimal"),
          definition = function(object, ...) { 
              return( list("optimalTx" = object@optimalTx,
                           "optimalY" = object@optimalY,
                           "type" = object@type ) )
            })
