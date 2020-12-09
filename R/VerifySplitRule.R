# Verify input 'splitRule'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'splitRule' is one of {'logrank', 'mean'}.
#
# successful methods return the original character object with possible
#  modification to make all lower case
#
setGeneric(name = ".VerifySplitRule",
           def = function(splitRule, ...) { 
                   standardGeneric(".VerifySplitRule") 
                 })

# the default method generates an error
setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "ANY"),
          definition = function(splitRule, ...) { 
              stop("splitRule must be one of {'logrank', 'mean'}", call. = FALSE)
            })

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "NULL"),
          definition = function(splitRule, ..., criticalValue) { 
              if (criticalValue == "mean") {
                splitRule = "mean"
              } else {
                splitRule = "logrank"
              }

              return( .VerifySplitRule(splitRule = splitRule, ...) ) 
            })

setMethod(f = ".VerifySplitRule",
          signature = c(splitRule = "character"),
          definition = function(splitRule, ...) { 

              splitRule <- tolower(x = splitRule)
              if (!(splitRule %in% c("logrank", "mean"))) {
                stop("splitRule must be one of {'logrank', 'mean'}")
              }

              return( splitRule ) 
            })
