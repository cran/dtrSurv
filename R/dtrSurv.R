#' Dynamic Treatment Regime for Survival Analysis
#'
#' Provides methods for estimating multi-stage optimal dynamic treatment
#'   regimes for survival outcomes with dependent censoring.
#'
#' If using a common formula for all decision points, i.e., 
#'   'models' is a single formula object, your data must follow a specific 
#'   format. Specifically, if 'stageLabel' = ".", covariates must be named as
#'   xxx.1 for the first decision point, 
#'   xxx.2 for the second, xxx.3 for the third, etc. The exact structure of the
#'   'xxx' can be generally defined; however, it cannot contain the stageLabel. 
#'   For example, if the column names are
#'   (Y.1, Y.2, d.1, d.2, A.1, A.2, X.1, X.2) 
#'   'models' = Surv(Y,d) ~ X + A would lead to 
#'   Surv(Y.1,d.1) ~ X.1 + A.1 as the first stage model; and 
#'   Surv(Y.2,d.2) ~ X.2 + A.2 
#'   as the second stage. 
#'   Further, baseline covariates can be used rather than stage dependent. In
#'   this case, the covariates should have no stageLabel.
#'   For example, if the column names are
#'   (Y.1, Y.2, d.1, d.2, A.1, A.2, X1, X2) where X1 and X2 are baseline
#'   'models' = Surv(Y,d) ~ X1 + X2 + A would lead to 
#'   Surv(Y.1,d.1) ~ X1 + X2 + A.1 as the first stage model; and 
#'   Surv(Y.2,d.2) ~ X1 + X2 + A.2 
#'   as the second stage. 
#'
#'   Y.k is the length of Stage k so that (Y.1 + Y.2 + ... + Y.K) is the 
#'   overall observed failure time, d.k is the censoring status at Stage k, 
#'   A.k is the treatment at Stage k, k=1,2,..., K. Note that every quantity 
#'   here is a stage-wide quantity. In other words, Y.2 is the length of Stage 2 
#'   and is not cumulative from the baseline. d.1 is 1 only if a subject 
#'   experiences failure during that stage, and 0 if he/she was censored 
#'   at Stage 1 or moved to Stage 2.
#'   
#'   When one experienced censoring or failure at Stage k, it should be that
#'   Y.j = 0 for all j > k and instantaneous failure (Y.k < 1e-8) is not allowed 
#'   E.g., when d.(k-1) = 0 and Y.k = 0, the person is considered censored at 
#'   Stage k-1, but when d.(k-1) = 0 and Y.k = 2, the person made it to Stage k
#'   and either experienced failure or censoring (depending on d.k) during Stage k.
#'   
#'   Any subject with missing values at Stage k will be ignored.
#' 
#' @param ... Ignored. Present only to require named inputs.
#'
#' @param data A data.frame object. The full dataset including treatments
#'   received, all stage covariates, observed times, and censoring
#'   indicators.
#'   Can be provided as a matrix object if column headers are included.
#'   Can contain missing data coded as NA, but cannot contain NaN.
#'
#' @param txName A character vector object. The treatment variable name for
#'   each decision point. Each element corresponds to the respective decision 
#'   point (element 1 = 1st decision; element 2 = 2nd decision, etc.).
#'
#' @param models A list object or a single formula. The models for each
#'   decision point. For list objects, each 
#'   element corresponds to the respective decision point.
#'   Each element contains a formula defining the 
#'   response as a Surv() object and the covariate structure of the model.
#'   Note that this model should not include any terms of order > 1. If
#'   using a single formula and the number of decision points is > 1, it is
#'   assumed that 'models' is a common formula to be used across all decision
#'   points. See details for further discussion.
#'
#' @param usePrevTime A logical object. If TRUE, previous times are included
#'   in the common formula model given in 'models'. This input is ignored if 
#'   'models' is not specified as a single common formula.
#'
#' @param timePoints A character object or a numeric vector object. If a character
#'   object, must be one of \{"quad", "uni", "exp"\} indicating the distribution
#'   from which the time points are to be calculated. For character input,
#'   input 'nTimes' must also be provided. If a numeric vector, the
#'   time points to be used. If 0 is not the first value, it will be
#'   concatenated by the software.
#'
#' @param nTimes An integer object. The total number of time points to be
#'   generated and considered. Used in conjunction with input 'timePoints'
#'   when 'timePoints' is a character; ignored otherwise.
#'
#' @param tau A numeric object or NULL. The study length. If NULL, the
#'   maximum timePoint is used. 
#'
#' @param criticalValue A character object. Must be one of 
#'   \{"mean", "surv.prob", "surv.mean"\}. The estimator for the value
#'   of a treatment rule. For "mean": the mean survival time; for
#'   "surv.prob": the mean survival probability at time 'evalTime';
#'   for "surv.mean": first the mean survival probability is used, if ties
#'   exist across treatments, the mean survival time is used to identify the
#'   optimal.
#'
#' @param evalTime A numeric object or NULL. If numeric, the time at which
#'   the survival probability is to be estimated to determine
#'   the optimal treatment rule; 'criticalValue' must be one of 
#'   \{"surv.prob", "surv.mean"\}. If NULL, 'criticalValue' must be \{"mean"\}.
#'
#' @param splitRule A character object OR NULL. 
#'   Must be one of \{"logrank", "mean"\}
#'   indicating the test used to determine an optimal split. If NULL and
#'   'criticalValue' = 'mean', takes value 'mean'. If NULL and
#'   'criticalValue' = 'surv.prob' or 'surv.mean', takes value 'logrank'.
#'
#' @param ERT A logical object. If TRUE, the Extremely Randomized Trees 
#'   algorithm is used to select the candidate variable.
#'
#' @param sampleSize A numeric object, numeric vector object, or NULL. 
#'   The fraction (0 < sampleSize <= 1) of the data to be used for each 
#'   tree in the forest. If only
#'   one value is given, it is assumed to be the fraction for all decision
#'   points. If a vector is given, the length must be equal to the total
#'   number of decision points and each element corresponds to its respective
#'   decision point. If NULL and 'ERT' is TRUE,
#'   sampleSize defaults to 1.0. If NULL and 'ERT'
#'   is FALSE, sampleSize defaults to 0.632.
#'
#' @param uniformSplit A logical object. If 'ERT' and 'uniformSplit' are TRUE,
#'   the random cutoff is sampled from a uniform distribution over the range
#'   of available covariate values. If 'ERT' is TRUE and 'uniformSplit' is 
#'   FALSE, a case is randomly selected and the cutoff is taken to be the mean
#'   cutoff between it and the next largest covariate value. If 'ERT' is FALSE,
#'   input is ignored.
#'
#' @param replace A logical object or NULL. If TRUE, the sample drawn for each 
#'   of the nTree trees may have duplicate records. If FALSE, no individual is 
#'   present in the sample for than once. If NULL, 'replace' = !'ERT'.
#'
#' @param randomSplit A numeric object. The probability that a random split
#'   will occur. Must be 0 < randomSplit < 1.
#'
#' @param tieMethod A character object. Must be one of 
#'   \{"first", "random"\}. If multiple splits lead to the same
#'   value, the method by which the tie is broken.
#'
#' @param minEvent An integer object. The minimum number of events that must be
#'   present in a node.
#'
#' @param nodeSize An integer object. The minimum number of individuals that
#'   must be present in a node.
#'
#' @param nTree An integer object. The number of trees to grow.
#'
#' @param mTry An integer or integer vector object. The maximum number of 
#'   covariates to sample for each split. If a vector, each element
#'   corresponds to its respective decision point.
#'
#' @param pooled A logical object. If TRUE, data are pooled for the analysis.
#'    If FALSE, data is separated into groups based on treatment
#'    received and a tree is grown for each treatment group.
#'
#' @param stratifiedSplit A numeric object. The stratified random split
#'    coefficient. Covariates for which the number of splits (s_i) is less
#'    than s*stratifiedSplit/d are explored preferentially 
#     (s is the total number of splits, d is the
#'    total number of covariates under consideration).
#'
#' @param stageLabel A character object. If using a common formula, the 
#'    character used to separate the covariate from the decision point label.
#'    See details.
#'
#' @references Cho, H., Holloway, S.T., and Kosorok, M.R.
#'   Multi-stage optimal dynamic treatment regimes for survival outcomes
#'   with dependent censoring. Submitted.
#'
#' @include VerifyData.R VerifyTxName.R VerifyModels.R VerifySampleSize.R
#' @include VerifyUsePrevTime.R class_Parameters.R
#' @include class_DTRSurvStep.R class_DTRSurv.R 
#' @import methods
#' @export
#' @useDynLib dtrSurv
#' @import survival
#'
#' @returns An S4 object of class DTRSurv containing the key results and
#'   input parameters of the analysis. The information contained therein
#'   should be accessed through convenience functions stage(), show(), print(),
#'   and predict().
#'
#' @examples
#'
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(1:100,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X.1" = rnorm(100), "X.2" = rnorm(100))
#'
#' # responses must be zero after event
#' evt <- dt[,"D.1"] == 1L
#' dt[evt, "Y.2"] <- 0.0
#'
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = list(Surv(Y.1,D.1)~X.1+A.1, Surv(Y.2,D.2)~X.2+A.2+Y.1))
#'
#' # common formula
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = Surv(Y,D)~X+A,
#'         usePrevTime = TRUE,
#'         stageLabel = ".")
#'
#' # common formula and pooled analysis
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = Surv(Y,D)~X+A,
#'         usePrevTime = TRUE,
#'         stageLabel = ".",
#'         pooled = TRUE)
#'
#' dt <- data.frame("Y.1" = sample(1:100,100,TRUE), "Y.2" = sample(1:100,100,TRUE),
#'                  "D.1" = rbinom(100, 1, 0.9), "D.2" = rbinom(100,1,0.9),
#'                  "A.1" = rbinom(100, 1, 0.5), "A.2" = rbinom(100,1,0.5),
#'                  "X1" = rnorm(100), "X2" = rnorm(100))
#'
#' # responses must be zero after event
#' evt <- dt[,"D.1"] == 1L
#' dt[evt, "Y.2"] <- 0.0
#'
#' # common formula with only baseline covariates
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = Surv(Y,D)~X1+X2+A)
#'
#' # common formula with only baseline covariates
#' # cutoff selected from indices
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = Surv(Y,D)~X1+X2+A,
#'         ERT = TRUE, uniformSplit = FALSE)
#'
#' # common formula with only baseline covariates
#' # not extremely random trees
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = Surv(Y,D)~X1+X2+A,
#'         ERT = FALSE)
#'
#' # common formula with only baseline covariates
#' # survival probability
#' dtrSurv(data = dt, 
#'         txName = c("A.1", "A.2"), 
#'         models = Surv(Y,D)~X1+X2+A,
#'         criticalValue = 'surv.prob')
#'

dtrSurv <- function(data,
                    txName,
                    models, 
                    ...,
                    usePrevTime = TRUE,
                    timePoints = "quad",
                    nTimes = 100L,
                    tau = NULL,
                    criticalValue = "mean",
                    evalTime = NULL,
                    splitRule = NULL,
                    ERT = TRUE,
                    uniformSplit = NULL,
                    sampleSize = NULL,
                    replace = NULL,
                    randomSplit = 0.2,
                    tieMethod = "random",
                    minEvent = 3L,
                    nodeSize = 6L, 
                    nTree = 10L,
                    mTry = NULL,
                    pooled = FALSE,
                    stratifiedSplit = NULL,
                    stageLabel = ".") {

  # ensure that 'data' is provided as a data.frame or a matrix and does not
  # contain NaN values. If 'data' is appropriate, method returns a data.frame
  # object
  data <- .VerifyData(data = data)

  # total number of individuals in dataset
  nSamples <- nrow(x = data)

  # ensure that 'txName' is provided as a character or character vector and
  # that the provided names are present in 'data'. This input defines the
  # number of decision points for the analysis. If 'txName' is appropriate,
  # the object returned is the original input without modification.
  txName <- .VerifyTxName(txName = txName, data = data)

  # number of decision points in the analysis
  nDP <- length(x = txName)

  # ensures that the usePrevTime input is of appropriate type
  usePrevTime <- .VerifyUsePrevTimes(usePrevTimes = usePrevTime)

  # ensure that 'models' is provided as a formula or a list of formula and
  # that the provided models can be generated by the data. If the input is 
  # appropriate, the object returned is list containing
  #   "models" - the original input.
  #   "response" - matrix of the survival response variables
  models <- .VerifyModels(models = models, 
                          nDP = nDP, 
                          data = data, 
                          txName = txName,
                          stageLabel = stageLabel,
                          usePrevTime = usePrevTime)

  response <- models$response
  del <- models$delta
  models <- models$models

  # Ensure that responses are zero after event
  if (ncol(x = del) > 1L) {
    for (i in 1L:ncol(x = del)) {
      evt <- del[,i] > 0.5
      j <- i + 1L
      while (j <= ncol(x = response)) {
        if (any(abs(response[evt,j]) > 1e-8)) {
          stop("for ", sum(abs(response[evt,j]) > 1e-8), 
               " participant(s), delta = 1 at stage ", i, 
               "; non-zero stage response found in stage ", j,
               call. = FALSE)
        }
        j <- j + 1L
      }
    }
  }

  # convert to list for single decision analyses for convenience
  if (!is.list(x = models)) models <- list(models)

  # combine all inputs that regulate tree and specify analysis preferences
  # function returns a Parameters object

  params <- .parameters(timePoints = timePoints,
                        tau = tau,
                        nTimes = nTimes,
                        response = response,
                        nTree = nTree, 
                        ERT = ERT, 
                        uniformSplit = uniformSplit, 
                        randomSplit = randomSplit, 
                        splitRule = splitRule,
                        replace = replace, 
                        nodeSize = nodeSize, 
                        minEvent = minEvent, 
                        tieMethod = tieMethod,
                        criticalValue = criticalValue, 
                        survivalTime = evalTime,
                        nSamples = nSamples,
                        pooled = pooled,
                        stratifiedSplit = stratifiedSplit)

  # store basic information on fortran side

  # retrieve index and fraction if survival type value estimator
  # set as 0's if mean estimator

  crit <- .CriticalValueCriterion(params)
  if (crit == "mean") {
    ind = 0L
    frac = 0.0
  } else if (crit == "surv.mean" || crit == "surv.prob") {
    ind = params@sIndex
    frac = params@sFraction
  }

  # set basic parameter values in Fortran
  res = .Fortran("setUpBasics",
                 t_nt = as.integer(x = .NTimes(object = params)),
                 t_dt = as.double(x = .TimeDiff(object = params)),
                 t_rs = as.double(x = params@randomSplit),
                 t_ERT = as.integer(x = params@ERT),
                 t_uniformSplit = as.integer(x = params@uniformSplit),
                 t_nodeSize = as.integer(x = .NodeSize(object = params)),
                 t_minEvent = as.integer(x = .MinEvent(object = params)),
                 t_rule = as.integer(x = params@splitRule == 'logrank'),
                 t_sIndex = as.integer(x = ind),
                 t_sFraction = as.double(x = frac),
                 t_stratifiedSplit = as.double(x = params@stratifiedSplit),
                 t_replace = as.integer(params@replace),
                 PACKAGE = "dtrSurv")

  # ensure that if given, sampleSize is 0 < sampleSize <= 1 and that
  # a value is provided for each decision point. If only 1 value is given,
  # it is assumed to be used for all decision points
  sampleSize <- .VerifySampleSize(sampleSize = sampleSize, 
                                  ERT = params@ERT,  
                                  nDP = nDP)

  # ensure that mTry is provided as a vector. At this point, there is
  # no verification of an appropriate value
  if (length(x = mTry) == 1L) {
    mTry <- rep(x = mTry, times = nDP)
  } else if (!is.null(x = mTry) && {length(x = mTry) != nDP}) {
    stop("if provided as vector, mTry must be provided for each dp", 
         call. = FALSE)
  }

  # Q-learning algorithm

  stageResults <- list()

  # final stage analysis

  message("Stage ", nDP)

  stageResults[[ nDP ]] <- .dtrSurvStep(model = models[[ nDP ]], 
                                        data = data, 
                                        priorStep = NULL,
                                        params = params,
                                        txName = txName[nDP],
                                        mTry = mTry[nDP],
                                        sampleSize = sampleSize[nDP])

  q <- nDP - 1L
  
  # backward recursion

  while( q >= 1L ) {

    message("Stage ", q)

    stageResults[[ q ]] <- .dtrSurvStep(model = models[[ q ]], 
                                        data = data, 
                                        priorStep = stageResults[[ q + 1L ]], 
                                        params = params,
                                        txName = txName[q],
                                        mTry = mTry[q],
                                        sampleSize = sampleSize[q])

    q <- q - 1L
  }

  # value obtained from the first stage analysis
  valueTrain <- .meanValue(object = stageResults[[ 1L ]])

  message("Estimated Value:", appendLF = FALSE)
  for (i in 1L:length(valueTrain)) {
    message(" ", names(valueTrain)[i], ": ", valueTrain[[ i ]], appendLF = FALSE)
  }
  message()


  # store values in call structure for returned object
  cl <- match.call()
  cl[[ 1L ]] <- as.name("dtrSurv")

  res <- new(Class = "DTRSurv", 
             "stageResults" = stageResults,
             "value" = valueTrain,
             "call" = cl,
             "params" = params)

  return( res )

}


#' Hidden methods
#'
#' @name dtrSurv-internal-api
#' @keywords internal
#' @import methods
NULL

