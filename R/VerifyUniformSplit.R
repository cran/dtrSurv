# Verify input 'uniformSplit'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'uniformSplit' is provided logical or NULL. If NULL set to value
#  of ERT 
#
# successful methods return a logical object.
#
setGeneric(name = ".VerifyUniformSplit",
           def = function(uniformSplit, ...) { 
                   standardGeneric(".VerifyUniformSplit") 
                 })

# the default method generates an error
setMethod(f = ".VerifyUniformSplit",
          signature = c(uniformSplit = "ANY"),
          definition = function(uniformSplit, ...) { 
              stop("uniformSplit must be logical or NULL", call. = FALSE)
            })

setMethod(f = ".VerifyUniformSplit",
          signature = c(uniformSplit = "logical"),
          definition = function(uniformSplit, ...) { 
              if (is.na(x = uniformSplit)) {
                stop("uniformSplit must be logical or NULL", call. = FALSE)
              }
              return( uniformSplit ) 
            })

setMethod(f = ".VerifyUniformSplit",
          signature = c(uniformSplit = "NULL"),
          definition = function(uniformSplit, ..., ERT) { 
              if (is.null(x = ERT)) {
                stop("if uniformSplit = NULL, ERT must be set", call. = FALSE)
              }
              return( .VerifyUniformSplit(uniformSplit = ERT, ...) ) 
            })
