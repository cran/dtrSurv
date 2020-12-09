# @param timePoints A numeric vector. The timepoints used in the analysis
#
# @param survVector A numeric vector. The survival function for a single 
#   individual
#
# @param by A numeric. The amount of time by which the survival function is
#   to be shifted
#
# Function returns the shifted survival function
.shift <- function(timePoints, survVector, by) {

  # the number of timepoints in the analysis
  nTimes <- length(x = timePoints)

  # shift the timepoints by the specified amount
  timePrime <- timePoints - by

  # initialize the shifted survival function
  survPrime <- numeric(length = nTimes)

  for (i in 1L:nTimes) {

    # for each shifted time value, count the number of timepoints <= timeprime_i
    sIndex <- sum(timePoints <= timePrime[i])

    if (sIndex < nTimes && sIndex > 0L) {
      # if the number of timepoints <= timeprime_i is not all or none
      # determine the fraction of dt[i] to which the shift corresponds
      sFraction <- {timePrime[i] - timePoints[sIndex]} / 
                   {timePoints[sIndex + 1L] - timePoints[sIndex]}

      # interpolate the survival function to the shifted time value
      survPrime[i] <- survVector[sIndex] * (1.0-sFraction) + 
                      survVector[sIndex+1L] * sFraction

    } else if (sIndex == 0L) {
      # if no timepoints are <= to the shifted time, survival set to 1
      survPrime[i] <- 1.0

    } else if (sIndex == nTimes) {
      # it should never happen that all timepoints are <= shifted times
      stop("this condition should not happen -- contact maintainer")
    }

    if (i > 1L && survPrime[i] > survPrime[i-1L]) survPrime[i] <- survPrime[i-1L]

  }

  return( survPrime )
}

.shiftMat <- function(timePoints, 
                      survMatrix, 
                      shiftVector,  
                      surv2prob) {

  # number of timepoints
  nTimes <- length(x = timePoints)

  # number of individuals in subset
  nSamples <- length(x = shiftVector)

  # default survival function to zero
  survShifted <- matrix(data = 0.0, nrow = nTimes, ncol = nSamples)

  for (i in 1L:nSamples) {
    # for each individual, shift their survival function down in time by
    # the specified amount given in shiftVector
    survShifted[,i] <- .shift(timePoints = timePoints, 
                              survVector = survMatrix[,i], 
                              by = shiftVector[i])
  }

  if (surv2prob) {
    # if survival function is to be converted to probability mass vector
    # return the change in the survival function at each time step
    survShifted <- survShifted - rbind(survShifted[-1L,], 0.0)
  }

  return( survShifted )
}
