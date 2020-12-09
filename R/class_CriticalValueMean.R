# Class extends CriticalValue to indicate that critical value is non-survival
#  mean
#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueCriterion(object, ...) {defined}
#  .CreateValueObject(object, ...) {defined}
#
# Functions
#  .criticalValueMean(...)
#
#' @include class_CriticalValue.R
setClass(Class = "CriticalValueMean",
         contains = c("CriticalValueBase"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method returns a character (specifically "mean")
#-------------------------------------------------------------------------------
setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "CriticalValueMean"),
          definition = function(object, ...) { return( "mean" ) })

#-------------------------------------------------------------------------------
# internal function to create an object of class CriticalValueMean
#-------------------------------------------------------------------------------
# function returns a CriticalValueMean object
#-------------------------------------------------------------------------------
.criticalValueMean <- function(...) { return( new("CriticalValueMean") ) }
