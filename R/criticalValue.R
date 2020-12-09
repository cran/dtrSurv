# Function to verify inputs and create appropriate CriticalValue object
#
#  Function is not exported and for internal convenience only
#
# Function returns an object of class CriticalValueMean or CriticalValueSurvival
#
#' @include class_CriticalValue.R class_CriticalValueMean.R
#' @include class_CriticalValueSurvival.R
#' @include VerifyCriticalValue.R VerifySurvivalTime.R
.criticalValue <- function(criticalValue,
                           survivalTime,
                           tau,
                           timePoints) {

  # ensure criticalValue is one of {'mean', 'surv.prob', 'surv.mean'}. 
  # Methods return the original character possibly modified to be lower case.
  criticalValue <- .VerifyCriticalValue(criticalValue = criticalValue)


  # ensure that survivalTime is provided if criticalValue is one of 
  # {'surv.prob', 'surv.man'}. Methods return the numeric survivalTime or NULL.
  survivalTime <- .VerifySurvivalTime(survivalTime = survivalTime,
                                      criticalValue = criticalValue,
                                      tau = tau)

  if (!is.null(x = survivalTime)) {
    # if survivalTime is given as input, verify it and the timePoints input and
    # create a CriticalValueSurvival object
    return( .criticalValueSurvival(survivalTime = survivalTime,
                                   timePoints = timePoints,
                                   type = criticalValue) )
  } else {
    # if survivalTime is not given as input, create a CriticalValueMean object
    return( .criticalValueMean() )
  }

}
