# Virtual class to store information regarding critical value selection
#
# Class is not exported and is for internal convenience only
#
# Methods
#  .CriticalValueAsList(object, ...) {not allowed}
#  .CriticalValueCriterion(object, ...) {not allowed}
#  .CreateValueObject(object, ...) {not allowed}
#
setClass(Class = "CriticalValueBase",
         contains = c("VIRTUAL"))

#-------------------------------------------------------------------------------
# method to identify if critical value is a mean or a probability
#-------------------------------------------------------------------------------
# method is not defined for a general CriticalValue object
#-------------------------------------------------------------------------------
setGeneric(name = ".CriticalValueCriterion",
           def = function(object, ...) { standardGeneric(".CriticalValueCriterion") })

setMethod(f = ".CriticalValueCriterion",
          signature = c(object = "ANY"),
          definition = function(object, ...) { stop("not allowed") })

resample <- function(x, ...) x[sample.int(n = length(x = x), ...)]

