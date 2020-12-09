# Class to store parameters that regulate tree and specify analysis preferences
#
# Class is not exported and is for internal convenience only
#
#  @slot nTree An integer object. The number of trees to be generated in forest.
#
#  @slot nodeSize An integer object. The minimum number of cases allowed in a
#    node
#
#  @slot minEvent An integer object. The minimum number of events allowed in a 
#    node
#
#  @slot pooled A logical object. TRUE indicates that treatment groups are to
#    be considered in combination.
#
#  @slot stratifiedSplit A number object. The coefficient phi for stratified
#    random spliting
#
# Getters
#   .NTree(object, ...) {new; defined}
#   .NodeSize(object, ...) {new; defined}
#   .MinEvent(object, ...) {new; defined}
#   .Pooled(object, ...) {new; defined}
#   .StratifiedSplit(object, ...) {new; defined}
#
# Methods
#   .TreeConditionsAsList(object, ...) {new; defined}
#
# Functions
# .treeConditions(..., nTree, nodeSize, minEvent,
#                 pooled, stratifiedSplit)
#
setClass(Class = "TreeConditions",
         slots = c(nTree = "integer", 
                   nodeSize = "integer", 
                   minEvent = "integer", 
                   pooled = "logical",
                   stratifiedSplit = "numeric"))

## Getters

#-------------------------------------------------------------------------------
# Method to retrieve number of trees in the forest
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".NTree",
           def = function(object, ...) { 
                   standardGeneric(".NTree") 
                 })

setMethod(f = ".NTree",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NTree",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@nTree ) })

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of cases allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".NodeSize",
           def = function(object, ...) { standardGeneric(".NodeSize") })

setMethod(f = ".NodeSize",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".NodeSize",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@nodeSize ) })

#-------------------------------------------------------------------------------
# Method to retrieve minimum number of events allowed in a node
#-------------------------------------------------------------------------------
# Method returns an integer
#-------------------------------------------------------------------------------
setGeneric(name = ".MinEvent",
           def = function(object, ...) { standardGeneric(".MinEvent") })

setMethod(f = ".MinEvent",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".MinEvent",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@minEvent ) })


#-------------------------------------------------------------------------------
# Method to retrieve flag for pooled analysis
#-------------------------------------------------------------------------------
# Method returns a logical
#-------------------------------------------------------------------------------
setGeneric(name = ".Pooled",
           def = function(object, ...) { standardGeneric(".Pooled") })

setMethod(f = ".Pooled",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".Pooled",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@pooled ) })

#-------------------------------------------------------------------------------
# Method to retrieve coefficient for stratified analysis
#-------------------------------------------------------------------------------
# Method returns a numeric
#-------------------------------------------------------------------------------
setGeneric(name = ".StratifiedSplit",
           def = function(object, ...) { standardGeneric(".StratifiedSplit") })

setMethod(f = ".StratifiedSplit",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".StratifiedSplit",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { return( object@stratifiedSplit ) })

#-------------------------------------------------------------------------------
# Method to retrieve primary slots for printing
#-------------------------------------------------------------------------------
# Method returns a list containing 6 elements 
#   "nTree" an integer, the total number of trees in forest
#   "nodeSize" an integer, the minimum number of cases in a node
#   "minEvent" an integer, the minimum number of events in a node
#   "pooled" a logical indicating of treatments were pooled
#   "stratifiedSplit" the coefficient phi for stratified 
#-------------------------------------------------------------------------------
setGeneric(name = ".TreeConditionsAsList",
           def = function(object, ...) { 
                   standardGeneric(".TreeConditionsAsList") 
                 })

setMethod(f = ".TreeConditionsAsList",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

setMethod(f = ".TreeConditionsAsList",
          signature = c(object = "TreeConditions"),
          definition = function(object, ...) { 
              return( list("nTree" = object@nTree, 
                           "nodeSize" = object@nodeSize,
                           "minEvent" = object@minEvent,
                           "pooled" = object@pooled,
                           "stratifiedSplit" = object@stratifiedSplit) )
            })

.treeConditions <- function(...,
                            nTree, 
                            nodeSize, 
                            minEvent, 
                            pooled,
                            stratifiedSplit) {

  # minimum number of events must be integer and > 0
  if (!is.numeric(x = minEvent)) stop("minEvent must be integer", call. = FALSE)
  minEvent <- as.integer(x = minEvent)
  if (minEvent < 1L) stop("minEvent must be non-zero positive", call. = FALSE)

  # minimum number of cases in each node, must be integer and > 0
  if (!is.numeric(x = nodeSize)) stop("nodeSize must be integer", call. = FALSE)
  nodeSize <- as.integer(x = nodeSize)
  if (nodeSize < 1L) stop("nodeSize must be non-zero positive", call. = FALSE)

  # number of trees to grow in forest, must be integer and > 0
  if (!is.numeric(x = nTree)) stop("nTree must be integer", call. = FALSE)
  nTree <- as.integer(x = nTree)
  if (nTree < 1L) stop("nTree must be non-zero positive", call. = FALSE)

  if (!is.logical(x = pooled) || is.na(x = pooled)) {
     stop("pooled must be logical", call. = FALSE)
  }

  if (is.null(x = stratifiedSplit) || stratifiedSplit <= 1e-8) {
    stratifiedSplit <- 0.0
  } else {
    if (stratifiedSplit > 1.0) {
      stop("stratifiedSplit must be [0,1]", call. = FALSE)
    }
  }

  return( new(Class = "TreeConditions",
              "nTree" = nTree, 
              "nodeSize" = nodeSize, 
              "minEvent" = minEvent, 
              "pooled" = pooled,
              "stratifiedSplit" = stratifiedSplit) )

}
