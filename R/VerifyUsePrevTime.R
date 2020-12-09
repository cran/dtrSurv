# Verify input 'usePrevTime'
#
# method is not exported and is for internal convenience only
#
# ensures that 'usePrevTime' input is logical. If FALSE, user has requested that
# previous times not be included in the model; methods return NULL. If nDP is
# 1L, no previous times are defined; methods return NULL. If TRUE and nDP > 1L,
# method returns a data.frame containing the cumulative response up to each
# decision point, i.e. R1 = 0, R2 = R1, R3 = R1 + R2.

setGeneric(name = ".VerifyUsePrevTimes",
           def = function(usePrevTimes, ...) { 
                   standardGeneric(".VerifyUsePrevTimes") 
                 })

#-------------------------------------------------------------------------------
# method to ensure the input 'usePrevTime' is appropriately defined and to
# obtain previous times for each decision point based on provided response
# matrix
#-------------------------------------------------------------------------------
# the default method generates an error
#-------------------------------------------------------------------------------
setMethod(f = ".VerifyUsePrevTimes",
          signature = c(usePrevTimes = "ANY"),
          definition = function(usePrevTimes, ...) { 
              stop("usePrevTime must be a logical", call. = FALSE)
            })

#-------------------------------------------------------------------------------
# method to ensure the input 'usePrevTime' is appropriately defined and to
# obtain previous times for each decision point based on provided response
# matrix
#-------------------------------------------------------------------------------
# the default method generates an error
#-------------------------------------------------------------------------------
setMethod(f = ".VerifyUsePrevTimes",
          signature = c(usePrevTimes = "logical"),
          definition = function(usePrevTimes, ...) { 

              # NA is a "logical" object -- eliminate it from being a possibilty
              if (is.na(x = usePrevTimes)) {
                stop("usePrevTime must be a logical", call. = FALSE)
              }
             
              return( usePrevTimes )

            })
