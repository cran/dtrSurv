# Verify input 'survivalTime'
#
# methods are not exported and are for internal convenience only
#
# ensures that 'survivalTime' is provided if criticalValue is one of 
# {'surv.prob', 'surv.mean'}. 
#
# successful methods return the numeric survivalTime or NULL.
#
setGeneric(name = ".VerifySurvivalTime",
           def = function(survivalTime, ...) { 
                   standardGeneric(".VerifySurvivalTime") 
                 })

# the default method generates an error
setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "ANY"),
          definition = function(survivalTime, ...) { 
              stop("evalTime must be numeric or NULL", 
                   call. = FALSE)
            })

setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "numeric"),
          definition = function(survivalTime, ..., criticalValue, tau) { 

              if (!{criticalValue %in% c("surv.prob", "surv.mean")}) {
                message("evalTime is ignored if critical value is mean")
                return( NULL )
              }

              if (length(x = survivalTime) > 1L) {
                stop("only 1 value for evalTime can be given",
                     call. = FALSE)
              }

              if (survivalTime <= 0.0 || survivalTime > tau) {
                stop("evalTime must be between 0 and tau", call. = FALSE)
              }

              message("evalTime ", survivalTime)

              return( survivalTime ) 
            })

setMethod(f = ".VerifySurvivalTime",
          signature = c(survivalTime = "NULL"),
          definition = function(survivalTime, ..., tau) { 

              return( .VerifySurvivalTime(survivalTime = tau/2.0,
                                          tau = tau, ... ) )
            })
