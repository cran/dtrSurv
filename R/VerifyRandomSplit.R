# Verify input 'randomSplit'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'randomSplit' is numeric satisfying 0 < rs < 1. 
#
# successful methods return the object
#
setGeneric(name = ".VerifyRandomSplit",
           def = function(randomSplit, ...) { 
                   standardGeneric(".VerifyRandomSplit") 
                 })

# the default method generates an error
setMethod(f = ".VerifyRandomSplit",
          signature = c(randomSplit = "ANY"),
          definition = function(randomSplit, ...) { 
              stop("randomSplit must obey 0 < randomSplit < 1", call. = FALSE)
            })

setMethod(f = ".VerifyRandomSplit",
          signature = c(randomSplit = "numeric"),
          definition = function(randomSplit, ...) { 

              if (length(x = randomSplit) > 1L) {
                stop("only 1 value for randomSplit can be given", call. = FALSE)
              }

              if (randomSplit <= 0.0 || randomSplit >= 1.0) {
                stop("randomSplit must obey 0 < randomSplit < 1", call. = FALSE)
              }

              return( randomSplit ) 
            })
