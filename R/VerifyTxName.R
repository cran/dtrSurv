# Verify input 'txName'
#
# method is not exported and is for internal convenience only
#
# ensures that 'txName' is provided as a character or character vector and
# that the provided names are present in data. This input defines the
# number of decision points for the analysis.
#
# successful methods return the original input without modification.
#
setGeneric(name = ".VerifyTxName",
           def = function(txName, ...) { standardGeneric(".VerifyTxName") })

# the default method generates an error
setMethod(f = ".VerifyTxName",
          signature = c(txName = "ANY"),
          definition = function(txName, ...) { 
              stop("txName must be a vector of character objects", 
                   call. = FALSE)
            })

setMethod(f = ".VerifyTxName",
          signature = c(txName = "character"),
          definition = function(txName, ..., data) { 

              if (length(x = txName) == 0L) {
                stop("txName must be provided", call. = FALSE)
              }

              # if treatment name is in data, it is *the* treatment name
              # if it is not, it is a single item indicating the treatment
              # variable name in the common formula
              # We assume a dot between name and decision point

              test <- tryCatch(expr = data[,txName,drop = FALSE],
                               error = function(e) { return( NULL ) })

              if (is.null(x = test) && length(x = txName) == 1L) {

                dataNames <- colnames(x = data)

                # split the data.frame names on dots
                cov <- strsplit(x = dataNames, split = ".", fixed = TRUE)

                # assume that first element is the common name
                areAs <- lapply(X = cov, FUN = function(x){x[[ 1L ]] == txName})

                if (sum(areAs) > 0L) {

                  nDP <- sum(areAs)

                  message("detected ", nDP, "decision points")

                  txName <- dataNames[areAs]

                  return( .VerifyTxName(txName = txName, data = data) )
                }
              }

              test <- tryCatch(expr = data[,txName,drop = FALSE],
                               error = function(e) { 
                                         stop("unable to retrieve 'txName' from data",
                                              e$message, call. = FALSE)
                                       })

              if (any(sapply(X = test, FUN = is.nan))) {
                stop("txName cannot include NaN values", call. = FALSE)
              }

              # ensure tx is factor or integer-like
              for (i in 1L:ncol(x = test)) {
                if (!is.factor(x = test[,i])) {
                  if (is.numeric(x = test[,i])) {
                    if (!isTRUE(all.equal(target = test[,i], 
                                          current = round(x = test[,i], digits = 0L)))) {

                      stop("treatment variable must be integer or factor", 
                           call. = FALSE)
                    }
                  } else {
                    stop("treatment variable must be integer or factor", 
                         call. = FALSE)
                  }
                }
              }

              return( txName )
            })
