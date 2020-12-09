# Verify input 'ERT'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'ERT' is provided as a logical or NULL object.
#
# successful methods return a logical indicating if the 
# Extremely Randomized Tree method is to be used
#
setGeneric(name = ".VerifyERT",
           def = function(ERT, ...) { standardGeneric(".VerifyERT") })

# the default method generates an error
setMethod(f = ".VerifyERT",
          signature = c(ERT = "ANY"),
          definition = function(ERT, ...) { 
              stop("ERT must be logical or NULL", call. = FALSE)
            })

setMethod(f = ".VerifyERT",
          signature = c(ERT = "logical"),
          definition = function(ERT, ...) { 
              if (is.na(x = ERT)) {
                stop("ERT must be logical or NULL", call. = FALSE)
              }
              return( ERT ) 
            })

setMethod(f = ".VerifyERT",
          signature = c(ERT = "NULL"),
          definition = function(ERT, ...) { return( TRUE ) })
