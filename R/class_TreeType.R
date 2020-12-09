# Class to store information regarding tree type ERT or Breiman
#
# Class is not exported and is for internal convenience only
#
#  @slot replace A logical object; TRUE indicates that samples can include
#    duplicate records
#
#  @slot randomSplit A numeric object; The probability of a random split
#
#  @slot ERT A logical object; TRUE indicates that extremely randomized trees
#    methods are to be used
#
#  @slot uniformSplit A logical object; TRUE indicates that cutoffs are to 
#    be selected based on a uniformed random number
#
#  @slot splitRule A character object; must be one of {'logrank', 'mean'}
#
#  @slot tieMethod A character object; must be one of {'first', 'random', 'NA'}
#
# Getters
#
# Methods
#
# Functions
#
#' @include VerifyERT.R VerifyRandomSplit.R 
#' @include VerifyUniformSplit.R VerifyReplace.R VerifySplitRule.R
#' @include VerifyTieMethod.R

setClass(Class = "TreeType",
         slots = c("replace" = "logical",
                   "randomSplit" = "numeric",
                   "ERT" = "logical",
                   "uniformSplit" = "logical",
                   "splitRule" = "character",
                   "tieMethod" = "character"))

## Getters

# initializer
.treeType <- function(ERT, 
                      nSamples,  
                      uniformSplit,  
                      replace,  
                      splitRule,  
                      tieMethod,
                      randomSplit,
                      criticalValue) {

  # ensure that ERT is logical or NULL. Methods return a logical.
  ERT <- .VerifyERT(ERT = ERT)

  # ensure that randomSplit is 0 <= rs < 1. Methods return a numeric.
  randomSplit <- .VerifyRandomSplit(randomSplit = randomSplit)

  # ensure that uniformSplit is logical or NULL. Methods return a logical.
  uniformSplit <- .VerifyUniformSplit(uniformSplit = uniformSplit, ERT = ERT)

  # ensure that replace is logical or NULL. Methods return a logical.
  replace <- .VerifyReplace(replace = replace, ERT = ERT)

  # verify splitRule. methods return the original character object with possible
  # modification to make all lower case
  splitRule <- .VerifySplitRule(splitRule = splitRule, 
                                criticalValue = criticalValue)

  # successful methods return the original character input possibly modified to
  # lower case
  tieMethod <- .VerifyTieMethod(tieMethod = tieMethod)

  return( new(Class = "TreeType", 
              "replace" = replace, 
              "randomSplit" = randomSplit,
              "ERT" = ERT,
              "uniformSplit" = uniformSplit,
              "splitRule" = splitRule,
              "tieMethod" = tieMethod) )

}
