# Verify input 'tieMethod'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'tieMethod' is in {'random', 'first'}
#
# successful methods return the original input possibly modified to be lower
#   case
#
setGeneric(name = ".VerifyTieMethod",
           def = function(tieMethod, ...) { 
                   standardGeneric(".VerifyTieMethod") 
                 })

# the default method generates an error
setMethod(f = ".VerifyTieMethod",
          signature = c(tieMethod = "ANY"),
          definition = function(tieMethod, ...) { 
              stop("tieMethod must be one of {'random', 'first'}", 
                   call. = FALSE)
            })

setMethod(f = ".VerifyTieMethod",
          signature = c(tieMethod = "character"),
          definition = function(tieMethod, ...) { 

              tieMethod <- tolower(x = tieMethod)
              if (!{tieMethod %in% c("random", "first")}) {
                stop("tieMethod must be one of {'random', 'first'}",
                     call. = FALSE)
              }

              return( tieMethod ) 
            })
