# Verify input 'sampleSize'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'sampleSize' is numeric, numeric vector, or NULL. If NULL, set 
#   to a default value based on ERT selection.
#
# successful methods return a numeric vector object
#
setGeneric(name = ".VerifySampleSize",
           def = function(sampleSize, ...) { 
                   standardGeneric(".VerifySampleSize") 
                 })

# the default method generates an error
setMethod(f = ".VerifySampleSize",
          signature = c(sampleSize = "ANY"),
          definition = function(sampleSize, ...) { 
              stop("sampleSize must be a numeric or NULL", call. = FALSE)
            })

setMethod(f = ".VerifySampleSize",
          signature = c(sampleSize = "numeric"),
          definition = function(sampleSize, ..., nDP) { 

              # sampleSize must be a fraction
              if (any({sampleSize < 1e-8 || sampleSize > 1.0})) {
                stop("sampleSize must be 0 < sampleSize <= 1", call. = FALSE)
              }

              # sampleSize must be provided for each decision point
              if (length(x = sampleSize) == 1L) {
                sampleSize <- rep(x = sampleSize, times = nDP)
              } else if (length(x = sampleSize) != nDP) {
                stop("if sampleSize provided as vector, ",
                     "must provide values for all decision points", 
                     call. = FALSE)
              }
                
              return( sampleSize ) 
            })

setMethod(f = ".VerifySampleSize",
          signature = c(sampleSize = "NULL"),
          definition = function(sampleSize, ..., ERT, nDP) { 

              if (is.null(x = ERT)) {
                stop("if sampleSize = NULL, ERT must be set", call. = FALSE)
              }

              sampleSize <- ifelse(test = ERT,
                                   yes = 1.0,
                                   no = 0.632)

              return( .VerifySampleSize(sampleSize = sampleSize, ...,
                                        nDP = nDP) ) 

            })
