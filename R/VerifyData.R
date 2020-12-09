# Verify input 'data'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'data' is provided as a data.frame or a matrix with named
# columns and that the object does not contain NaN values. 
#
# successful methods return a data.frame object containing the data
#
setGeneric(name = ".VerifyData",
           def = function(data, ...) { standardGeneric(".VerifyData") })

# the default method generates an error
setMethod(f = ".VerifyData",
          signature = c(data = "ANY"),
          definition = function(data, ...) { 
              stop("data must be a data.frame or a matrix with named columns",
                   call. = FALSE)
            })

setMethod(f = ".VerifyData",
          signature = c(data = "data.frame"),
          definition = function(data, ...) { 

              if (any(sapply(X = data, FUN = is.nan))) {
                stop("data cannot include NaN values", call. = FALSE)
              }

              return( data ) 
            })

setMethod(f = ".VerifyData",
          signature = c(data = "matrix"),
          definition = function(data, ...) { 

              if (is.null(x = colnames(x = data))) {
                stop("if a matrix, data must include column headers",
                     call. = FALSE)
              }

              return( .VerifyData(data = as.data.frame(x = data), ...) )
            })
