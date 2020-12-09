# Class to store parameters that regulate tree and specify analysis preferences
#
# Class is not exported and is for internal convenience only
#
# Methods
#   .ParametersAsList(object, ...) {new; defined}
#
# Functions
# .parameters(timePoints, nTimes, response, nTree, ERT, uniformSplit, 
#                      randomSplit, splitRule, replace, nodeSize, 
#                      minEvent, tieMethod, criticalValue, 
#                      survivalTime, nSamples, pooled, stratifiedSplit)
#
#' @include class_TimeInfo.R criticalValue.R class_TreeType.R
#' @include class_TreeConditions.R
setClass(Class = "Parameters_Mean",
         contains = c("TimeInfo", "CriticalValueMean", "TreeType", "TreeConditions"))

setClass(Class = "Parameters_Survival",
         contains = c("TimeInfo", "CriticalValueSurvival", "TreeType", "TreeConditions"))

setClassUnion(name = "Parameters",
              members = c("Parameters_Mean", "Parameters_Survival"))

setMethod(f = "initialize",
         signature = c(.Object = "Parameters_Mean"),
         def = function(.Object, ...) {

                   obj <- list(...)

                   for (i in 1L:length(x = obj)) {
                     if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                       as(.Object, "TimeInfo") <- obj[[ i ]]
                     } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                       as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
                     } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                       as(.Object, "TreeType") <- obj[[ i ]]
                     } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                       as(.Object, "TreeConditions") <- obj[[ i ]]
                     } else {
                       stop("unrecognized object sent to Paremeters object")
                     }
                    }
                   return( .Object )

                 })

setMethod(f = "initialize",
         signature = c(.Object = "Parameters_Survival"),
         def = function(.Object, ...) {

                   obj <- list(...)

                   for (i in 1L:length(x = obj)) {
                     if (is(object = obj[[ i ]], class2 = "TimeInfo")) {
                       as(.Object, "TimeInfo") <- obj[[ i ]]
                     } else if (is(object = obj[[ i ]], class2 = "CriticalValueBase")) {
                       as(.Object, is(object = obj[[ i ]])[1L]) <- obj[[ i ]]
                     } else if (is(object = obj[[ i ]], class2 = "TreeType")) {
                       as(.Object, "TreeType") <- obj[[ i ]]
                     } else if (is(object = obj[[ i ]], class2 = "TreeConditions")) {
                       as(.Object, "TreeConditions") <- obj[[ i ]]
                     } else {
                       stop("unrecognized object sent to Paremeters object")
                     }
                    }

                   return( .Object )

                 })

#-------------------------------------------------------------------------------
# Function to verify inputs and create a Parameters object
#-------------------------------------------------------------------------------
# Function returns a Parameters object
#-------------------------------------------------------------------------------
.parameters <- function(timePoints,
                        tau,
                        nTimes,
                        response,
                        nTree, 
                        ERT, 
                        uniformSplit, 
                        randomSplit, 
                        splitRule,
                        replace, 
                        nodeSize, 
                        minEvent, 
                        tieMethod,
                        criticalValue, 
                        survivalTime,
                        nSamples,
                        pooled,
                        stratifiedSplit) {


  # initialize TimeInfo
  # function returns a TimeInfo object
  timeInfo <- .timeInfo(timePoints = timePoints,
                        tau = tau,
                        nTimes = nTimes,
                        response = response)

  cv <- tolower(criticalValue)

  # initialize CriticalValue
  # function returns an object of class CriticalValueMean or
  #   CriticalValueSurvival depending on input survivalTime
  criticalValue <- .criticalValue(criticalValue = criticalValue,
                                  survivalTime = survivalTime,
                                  tau = .Tau(object = timeInfo),
                                  timePoints = .TimePoints(object = timeInfo))

  # initialize tree type info
  # function returns a TreeType object
  treeType <- .treeType(ERT = ERT,
                        nSamples = nSamples,
                        uniformSplit = uniformSplit,
                        replace = replace,
                        randomSplit = randomSplit,
                        splitRule = splitRule,
                        tieMethod = tieMethod,
                        criticalValue = cv)

  # initialize tree conditions info
  # function returns a TreeConditions object
  treeConditions <- .treeConditions(nTree = nTree, 
                                    nodeSize = nodeSize, 
                                    minEvent = minEvent, 
                                    pooled = pooled,
                                    stratifiedSplit = stratifiedSplit)

  if (is(object = criticalValue, class2 = "CriticalValueSurvival")) {
    return( new(Class = "Parameters_Survival",
                timeInfo,
                criticalValue,
                treeType,
                treeConditions) )
  } else {

    return( new(Class = "Parameters_Mean",
                timeInfo,
                criticalValue,
                treeType,
                treeConditions) )
  }

}
