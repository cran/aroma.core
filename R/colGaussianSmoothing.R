# M_i' = w*M = w*(T-R) = w*T - w*R = T_i' - R'
#
# Before smoothing, the reference R_i == median(T_i). 
# Keep this property for R' too.
#
# R' = median(T_i')
# T_i' = M_i' - R'
#
# => w*T = w*M + w*R = M' + w*R
setMethodS3("colGaussianSmoothing", "matrix", function(Y, x=seq_len(nrow(Y)), w=NULL, xOut=x, sd=1, censorSd=3, na.rm=FALSE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'Y'
  n <- nrow(Y)
  k <- ncol(Y)
  
  # Argument 'x'
  if (length(x) != n) {
    throw("Argument 'x' has different number of values than rows in 'Y': ", 
                                                     length(x), " != ", n)
  }

  # Argument 'w'
  if (is.null(w)) {
  } else {
    if (length(w) != n) {
      throw("Argument 'w' has different number of values that rows in 'Y': ", 
                                                       length(w), " != ", n)
    }
  }

  # Argument 'xOut'
  if (is.null(xOut)) {
    xOut <- x
  } else {
    xOut <- Arguments$getNumerics(xOut)
  }
  nOut <- length(xOut)

  # Argument 'censorSd':
  censorSd <- Arguments$getNumeric(censorSd, range=c(0,Inf))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Smoothing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Allocate vector of smoothed signals
  naValue <- NA_real_
  Ys <- matrix(naValue, nrow=nOut, ncol=k)
  colnames(Ys) <- colnames(Y)

#  wKernelMax <- dnorm(x, sd=sd)

  verbose && enter(verbose, "Estimating signals at given locations")

  verbose && cat(verbose, "Output locations:")
  verbose && str(verbose, xOut)

  censorThreshold <- censorSd * sd
  isCensored <- (censorThreshold < Inf)
  if (isCensored) {
    # Default weights - all zeros
    wZeroKernel <- rep(0, times=length(x))
  }

  # At each position in 'xOut', calculate the weighed average 
  # using a Gaussian kernel.
  for (kk in seq_len(nOut)) {
    if (kk %% 100 == 0)
      verbose && cat(verbose, kk)

    # Weights centered around x[kk]
    xDiff <- (x-xOut[kk])
    if (isCensored) {
      keep <- which(abs(xDiff) <= censorThreshold)

      # Nothing to do?
      if (length(keep) == 0) {
        next
      }

      wKernel <- dnorm(xDiff[keep], mean=0, sd=sd)
      Y2 <- Y[keep,,drop=FALSE]
    } else {
      wKernel <- dnorm(xDiff, mean=0, sd=sd)
      Y2 <- Y
    }

    wKernel <- wKernel / sum(wKernel)
#    verbose && str(verbose, wKernel)

    # Exclude NAs
    if (na.rm) {
      value <- colSums(wKernel*Y2, na.rm=TRUE)
    } else {
      value <- colSums(wKernel*Y2)
    }

    # Fix: Smoothing over a window with all missing values give zeros, not NA.
    idxs <- which(value == 0)
    if (length(idxs) > 0) {
      # Are these real zeros or missing values?
      Y2 <- Y2[idxs,,drop=FALSE]
      Y2 <- !is.na(Y2)
      idxsNA <- idxs[colSums(Y2) == 0]
      value[idxsNA] <- NA
    }

#    verbose && str(verbose, value)
    Ys[kk,] <- value
  } # for (kk ...)

  verbose && exit(verbose)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Weighting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  if (!is.null(w)) {
    w <- w / sum(w, na.rm=TRUE)
    Ys <- w*Ys
  }

  Ys
}) # colGaussianSmoothing()


setMethodS3("gaussianSmoothing", "matrix", function(Y, ...) {
  colGaussianSmoothing(Y, ...)
})

setMethodS3("gaussianSmoothing", "numeric", function(y, ...) {
  y <- colGaussianSmoothing(as.matrix(y), ...)
  dim(y) <- NULL
  y
})
