# Verify input 'replace'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'replace' is logical or NULL. If NULL, replace = ERT.
#
# successful methods return a logical indicating if replacement should be
# used in sampling
#
setGeneric(name = ".VerifyReplace",
           def = function(replace, ...) { standardGeneric(".VerifyReplace") })

# the default method generates an error
setMethod(f = ".VerifyReplace",
          signature = c(replace = "ANY"),
          definition = function(replace, ...) { 
              stop("replace must be logical or NULL", call. = FALSE)
            })

setMethod(f = ".VerifyReplace",
          signature = c(replace = "logical"),
          definition = function(replace, ...) { 
              if (is.na(x = replace)) {
                stop("replace must be logical or NULL", call. = FALSE)
              }
              return( replace ) 
            })

setMethod(f = ".VerifyReplace",
          signature = c(replace = "NULL"),
          definition = function(replace, ..., ERT) { 
              if (is.null(x = ERT)) {
                stop("if replace = NULL, ERT must be set", call. = FALSE)
              }
              return( .VerifyReplace(replace = !ERT, ...) ) 
            })
