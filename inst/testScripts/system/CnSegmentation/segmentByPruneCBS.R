setMethodS3("findAtomicRegions", "CopyNumberRegions", function(cnr, rcn, alpha=0.02, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rcn':
  rcn <- Arguments$getInstanceOf(rcn, "RawCopyNumbers");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  nbrOfRegions <- nbrOfRegions(cnr);

  # Nothing to do?
  if (nbrOfRegions < 3) {
    res <- list(
      atomicRegions=integer(0),
      atomicIslands=integer(0)
    );
    return(res);
  }

  verbose && enter(verbose, "Call equivalent copy-number states by pruning");

  # Initial set of atomic regions
  atomicRegions <- NULL;

  start <- cnr$start;
  stop <- cnr$stop;
  for (rr in 2:(nbrOfRegions-1)) {
    verbose && enter(verbose, sprintf("Region #%d of %d", rr, nbrOfRegions-1L));

    xRangeL <- c(start[rr-1], stop[rr-1]);
    xRangeR <- c(start[rr+1], stop[rr+1]);
    rcnL <- extractRegion(rcn, region=xRangeL);
    rcnR <- extractRegion(rcn, region=xRangeR);
    yL <- getSignals(rcnL);
    yR <- getSignals(rcnR);
    yL <- yL[is.finite(yL)];
    yR <- yR[is.finite(yR)];

    # Test H0: muL = muR against H1: muL != muR.
    fit <- t.test(yL, yR, paired=FALSE, var.equal=TRUE, alternative="two.sided");
    t <- fit$statistic;
    p <- fit$p.value;
    isSignificant <- (p < alpha);

    isSame <- (!isSignificant);
    verbose && printf(verbose, "t=%.3f (p=%g), (L==R)=%s\n", t, p, isSame);

    if (isSame) {
      atomicRegions <- c(atomicRegions, rr);
      verbose && print(verbose, atomicRegions);
    }

    verbose && exit(verbose);
  } # for (rr ...)

  # Atomic islands = atomic regions that are not next 
  # to another atomic region
  dups <- which(diff(atomicRegions) == 1);
  if (length(dups) > 0) {
    dups <- c(dups, dups+1L);
    atomicIslands <- atomicRegions[-dups];
  } else {
    atomicIslands <- atomicRegions;
  }

  res <- list(
    atomicRegions=atomicRegions,
    atomicIslands=atomicIslands,
    ambigousRegions=setdiff(atomicRegions, atomicIslands)
  );

  verbose && exit(verbose);

  res;
}, protected=TRUE) # findAtomicRegions()



setMethodS3("segmentByPruneCBS", "RawGenomicSignals", function(this, ...,       strict=FALSE, normalMean=0.0, debugPlot=FALSE, ylim=c(0,6), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'normalMean':
  normalMean <- Arguments$getDouble(normalMean);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Segment using PruneCBS");

  if (debugPlot) {
    cols <- c("atomic island"="red", 
              "ambigous atomic region"="orange", 
              "extreme region"="purple");
  }

  res <- list();

  depth <- 1;
  cnr <- NULL;
  cn <- this;
  while(TRUE) {
    verbose && enter(verbose, sprintf("Segmentation depth %d", depth));
    cnrPrev <- cnr;
  
    fit <- segmentByCBS(cn);
    cnr <- extractCopyNumberRegions(fit);

    # Converged?
    if (equals(cnr, cnrPrev)) {
      cnr$type <- "constant";
      res[[depth]] <- cnr;
      verbose && cat(verbose, "Quitting, because nothing has changed since last round.");
      verbose && exit(verbose);
      break;
    }
  
  
    # Decrease the weights close to change points.
    # This will lower the risk for false change points
    # in succeeding segmentation iterations.
    # The weighting should be on the same "x scale" as the
    # segmentation operates on, that is, if it segments with
    # physical distances, then the weighting function should
    # be on that too, and if it segments with index distances,
    # then the weighting function should to that too.
    verbose && enter(verbose, "Updating locus-specific weights");
    w <- cn$w;
    if (is.null(w)) {
      w <- rep(1, times=nbrOfLoci(cn));
    }
    dx <- 100e3;
    regions <- as.data.frame(cnr)[,c("start","stop")];
    regions <- as.matrix(regions);
    regions2 <- regions;
    regions2[,"start"] <- regions2[,"start"] + dx;
    regions2[,"stop"] <- regions2[,"stop"] - dx;
    nok <- (regions2[,"start"] > regions[,"stop"]);
    regions2[nok,"start"] <- regions[nok,"stop"];
    nok <- (regions2[,"stop"] < regions[,"start"]);
    regions2[nok,"stop"] <- regions[nok,"start"];
    regions3 <- regions;
    regions3[,"stop"] <- regions2[,"start"];
    regions4 <- regions;
    regions4[,"start"] <- regions2[,"stop"];
    regions <- rbind(regions3, regions4);
    x <- getPositions(cn);
    for (kk in seq(length=nrow(regions))) {
      region <- regions[kk,,drop=TRUE];
      keep <- (region[1] <= x & x <= region[2]);
      w[keep] <- 0.2*w[keep];
    }
    cn$w <- w;
    rm(w);
    verbose && exit(verbose);
  
    if (debugPlot) {
      if (depth > 1 && (depth-1) %% np == 0) {
        readline("wait...");
      }
      plot(cn, col="#aaaaaa", ylim=ylim);
      drawLevels(cnr, col="black");
      verbose && print(verbose, cnr);
    }
  
    # Done?
    if (nbrOfRegions(cnr) == 1) {
      cnr$type <- "constant";
      res[[depth]] <- cnr;
      verbose && exit(verbose);
      break;
    }
  

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Find minimal region(s) to be pruned
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fitAR <- findAtomicRegions(cnr, rcn=cn, ..., verbose=less(verbose,10));
    verbose && print(verbose, fitAR);

    # (a) Any atomic islands?
    minimalRegions <- fitAR$atomicIslands;
    n <- length(minimalRegions);
    if (n > 0) {
      type <- "atomic island";
    } else {
      # (b) Any ambigous atomic islands?
      minimalRegions <- fitAR$ambigousRegions;
      n <- length(minimalRegions);
      if (n > 0) {
        type <- "ambigous atomic region";
        cnrT <- subset(cnr, fitAR$ambigousRegions);
        verbose && cat(verbose, "atomicSiblings:");
        verbose && print(verbose, cnrT);
  
        # Find the furthest away from the normal state
        data <- as.data.frame(cnrT);
        normalMean <- 0.0;
        dMu <- (data$mean - normalMean);
        idx <- which.max(abs(dMu));
        minimalRegions <- fitAR$atomicRegions[idx];
      } else {
        type <- "extreme region";
        # (c) Extreme regions...
        # Find the furthest away from the normal state
        data <- as.data.frame(cnr);
        dMu <- (data$mean - normalMean);
        minimalRegions <- which.max(abs(dMu));
      }
    }
  
    n <- length(minimalRegions);
    verbose && printf(verbose, "Selected %d minimal regions that are %ss:\n", n, type);
    verbose && print(verbose, minimalRegions);

    if (n > 1 && strict) {
      verbose && cat(verbose, "Keeping only one minimal region (argument strict=TRUE)");
      # TO DO: For now, just the first one. /HB 2010-07-24.
      minimalRegions <- minimalRegions[1];
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Record the identified minimal regions
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cnrX <- subset(cnr, minimalRegions);
    cnrX$type <- type;
    res[[depth]] <- cnrX;

    if (debugPlot) {
      drawLevels(cnrX, col=cols[type]);
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Drop minimal region(s)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    cnrT <- subset(cnr, -minimalRegions);
    cn <- extractRegions(cn, regions=cnrT);

    depth <- depth + 1;

    verbose && exit(verbose);
  } # while(...)

  types <- sapply(res, FUN=function(cnr) cnr$type);
  names(res) <- types;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Reverse history
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- rev(res);

  verbose && exit(verbose);

  res;
}) # segmentByPruneCBS()



############################################################################
# HISTORY:
# 2010-07-24
# o CLEAN UP: Now the notation of the code better repflect the algorithm.
# o Now findAtomicRegions() returns ambigous atomic regions too.
# o Added argument 'ylim'.
# 2010-07-20
# o Added argument 'debugPlot'.
# 2010-07-19
# o Added trial version of segmentByPruneCBS().
# o TO DO: Down-weight loci that were close to earlier 
#   change points in the succeeding segmentations.
# o Added prototype version of findAtomicRegions().
# o Added prototype version of callByPruning().
# o Created.
############################################################################
