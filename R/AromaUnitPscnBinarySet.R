###########################################################################/**
# @RdocClass AromaUnitPscnBinarySet
#
# @title "The AromaUnitPscnBinarySet class"
#
# \description{
#  @classhierarchy
#
#  An AromaUnitPscnBinarySet object represents a set of
#  @see "AromaUnitPscnBinaryFile"s with \emph{identical} chip types.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AromaUnitSignalBinarySet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author
#*/###########################################################################
setConstructorS3("AromaUnitPscnBinarySet", function(...) {
  extend(AromaUnitSignalBinarySet(...), c("AromaUnitPscnBinarySet", uses("CopyNumberDataSet")))
})


setMethodS3("byName", "AromaUnitPscnBinarySet", function(static, name, tags=NULL, ..., chipType=NULL, paths=c("totalAndFracBData"), pattern=".*,pscn[.]asb$") {
  suppressWarnings({
    path <- findByName(static, name=name, tags=tags, chipType=chipType,
                                           ..., paths=paths, mustExist=TRUE)
  })

  suppressWarnings({
    byPath(static, path=path, ..., pattern=pattern)
  })
}, static=TRUE)


setMethodS3("getAverageFile", "AromaUnitPscnBinarySet", function(this, name=NULL, prefix="average", indices="remaining", mean=c("median", "mean"), sd=c("mad", "sd"), na.rm=TRUE, g=NULL, h=NULL, ..., unitsPerChunk=ram*10^7/length(this), ram=1, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'this':
  if (length(this) == 0L) {
    throw("Cannot get average file. ", class(this)[1L], " is empty: ", getFullName(this))
  }

  # Argument 'mean':
  if (is.character(mean)) {
    mean <- match.arg(mean)
    meanName <- mean
    if (mean == "mean") {
      mean <- base::rowMeans
    } else if (mean == "median") {
      mean <- rowMedians
    }
  } else if (is.function(mean)) {
    meanName <- "customMean"
  } else {
    throw("Argument 'mean' must be either a character or a function: ", mode(mean))
  }

  # Argument 'sd':
  if (is.character(sd)) {
    sd <- match.arg(sd)
    sdName <- sd
    if (sd == "sd") {
      sd <- rowSds
    } else if (sd == "mad") {
      sd <- rowMads
    }
  } else if (is.function(sd)) {
    sdName <- "customSd"
  } else {
    throw("Argument 'sd' must be either a character or a function: ",
                                                           mode(sd))
  }

  # Argument 'name':
  if (is.null(name)) {
    key <- list(method="getAverageFile", class=class(this)[1],
                arrays=sort(getNames(this)), mean=meanName, sd=sdName)
    # assign mean and sd to an empty environment so that digest() doesn't
    # pick up any "promised" objects from the original environment.
    # A bit ad hoc, but it works for now. /2007-01-03
    key <- lapply(key, FUN=function(x) {
      if (is.function(x))
        environment(x) <- emptyenv()
      x
    })
    id <- getChecksum(key)
    field <- "signals"
    name <- sprintf("%s-%s-%s-%s,%s", prefix, field, meanName, sdName, id)
  }

  # Argument 'indices':
  df <- getOneFile(this)
  nbrOfUnits <- nbrOfUnits(df)
  if (force) {
    if (identical(indices, "remaining")) {
      indices <- NULL
    }
  }

  if (is.null(indices)) {
    indices <- 1:nbrOfUnits
  } else if (identical(indices, "remaining")) {
  } else {
    indices <- Arguments$getIndices(indices, max=nbrOfUnits)
  }

  # Argument 'unitsPerChunk':
  unitsPerChunk <- Arguments$getInteger(unitsPerChunk, range=c(1,Inf))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Retrieving average unit signals across ", length(this), " arrays")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create CEL file to store the average array
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a private filename (with a dot prefix) to make sure it is not
  # identified as a regular data file when the directory is scanned for files.
  ext <- getFilenameExtension(df)
  filename <- sprintf(".%s.%s", name, ext)
  if (is.null(this$.averageFiles))
    this$.averageFiles <- list()
  res <- this$.averageFiles[[filename]]

  # Has file been deleted since last time?
  if (!is.null(res) && !isFile(res)) {
    warning("Will recalculate average file, because it seems to have been deleted since last time: ", getPathname(res))
    res <- NULL
  }

  if (is.null(res)) {
    verbose && enter(verbose, "Searching for an existing file")

    # Searching for the output file in multiple directories
    path <- getPath(this)
    paths <- c(path)

    # Drop tags from root path?
    if (getOption(aromaSettings, "devel/dropRootPathTags", TRUE)) {
      path <- dropRootPathTags(path, depth=2, verbose=less(verbose, 5))
      paths <- c(paths, path)
      paths <- unique(paths)
    }

    verbose && cat(verbose, "Paths:")
    verbose && print(verbose, paths)
    verbose && cat(verbose, "Filename: ", filename)

    pathname <- NULL
    for (kk in seq_along(paths)) {
      path <- paths[kk]
      verbose && enter(verbose, sprintf("Searching path #%d of %d", kk, length(paths)))

      verbose && cat(verbose, "Path: ", path)
      pathnameT <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE)
      verbose && cat(verbose, "Pathname: ", pathnameT)
      if (isFile(pathnameT)) {
        pathname <- pathnameT
        verbose && cat(verbose, "Found an existing file.")
        verbose && exit(verbose)
        break
      }

      verbose && exit(verbose)
    } # for (kk ...)
    verbose && cat(verbose, "Located pathname: ", pathname)

    verbose && exit(verbose)


    if (isFile(pathname)) {
      verbose && enter(verbose, "Loading existing data file")
      res <- newInstance(df, pathname)
      verbose && exit(verbose)
    } else {
      verbose && enter(verbose, "Allocating empty data file to store average signals")
      ugp <- getAromaUgpFile(df)

      path <- paths[length(paths)]
      verbose && cat(verbose, "Path: ", path)
      verbose && cat(verbose, "Filename: ", filename)
      pathname <- Arguments$getWritablePathname(filename, path=path, mustNotExist=TRUE)

      pathnameT <- pushTemporaryFile(pathname, verbose=verbose)

      res <- df$allocateFromUnitAnnotationDataFile(udf=ugp, filename=pathnameT, verbose=less(verbose))
      naValue <- NA_real_
      res[,1L] <- naValue
      res[,2L] <- naValue

      pathname <- popTemporaryFile(pathnameT, verbose=verbose)

      # Don't forget to update 'res' accordingly
      # BTW, should there be a protected setPathname() or a popTemporaryFile().
      res$.pathname <- pathname

      verbose && exit(verbose)
    }

    this$.averageFiles[[filename]] <- res
  }

  verbose && print(verbose, res)

  pathname <- getPathname(res)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify which indices to use
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (identical(indices, "remaining")) {
    values <- res[,nbrOfColumns(res),drop=TRUE]
    indices <- which(is.na(values) | (values == 0))
    # Not needed anymore
    values <- NULL
  }

  nbrOfIndices <- length(indices)
  # Nothing more to do?
  if (nbrOfIndices == 0L) {
    verbose && exit(verbose)
    return(res)
  }

  verbose && cat(verbose, "Number of units to be updated: ", nbrOfIndices)

  # Get the pathnames of all data files to average
  pathnames <- getPathnames(this)
  nbrOfArrays <- length(pathnames)

  # Update foot with number of arrays averaged over
  footer <- readFooter(res)
  params <- footer$params
  if (length(params) == 0L) {
    srcFiles <- lapply(this, FUN=function(df) {
      list(
        fullname = getFullName(df),
        fileSize = getFileSize(df),
        checkSum = getChecksum(df)
      )
    })
    params <- list(
      meanName = meanName,
      sdName = sdName
    )
    srcDetails <- list(
      nbrOfFiles = length(srcFiles),
      checkSum = getChecksum(srcFiles)
    )
    footer$srcDetails <- srcDetails
    footer$params <- params
    writeFooter(res, footer)
  }

  # Since we might want to do this robustly, but also because we want to
  # estimate the standard deviation, for each unit we need all data across
  # arrays at once.  In order to this efficiently, we do this in chunks

  arrays <- seq_len(nbrOfArrays)
  naValue <- NA_real_
  lapplyInChunks(indices, function(idxs, ...) {
    verbose && enter(verbose, "Processing chunk")
    verbose && str(verbose, "Indices in chunk:")
    verbose && str(verbose, idxs)

    for (cc in seq_len(nbrOfColumns)) {
      verbose && enter(verbose, sprintf("Column #%d of %d", cc, nbrOfColumns))

      verbose && enter(verbose, "Reading data")
      X <- matrix(naValue, nrow=length(idxs), ncol=nbrOfArrays)
      for (kk in arrays) {
        df <- this[[kk]]
        X[,kk] <- df[idxs,cc, drop=TRUE]
      }
      verbose && exit(verbose)

      if (!is.null(g)) {
        verbose && enter(verbose, "Transforming data using y = g(x)")
        X <- g(X)
        verbose && exit(verbose)
      }

      verbose && enter(verbose, "Estimating averages and standard deviations")
      if (na.rm) {
        n <- base::apply(X, MARGIN=1, FUN=function(x) { sum(!is.na(x)) })
      }
      # Calculate the mean signal
      mu <- mean(X, na.rm=na.rm) # Special mean()!
      # Calculate the standard deviation of the signals
      sigma <- sd(X, mean=mu, na.rm=na.rm) # Special sd()!
      verbose && exit(verbose)

      if (!is.null(h)) {
        verbose && enter(verbose, "Back-transforming estimates using x = h(y)")
        mu <- h(mu)
        sigma <- h(sigma)
        verbose && exit(verbose)
      }

      # Write estimates to result file
      verbose && enter(verbose, "Writing estimates")
      res[idxs,cc] <- mu
      ## Only mu is supported. \HB 2009-11-19 ##, stdvs=sigma, pixels=n)
      verbose && exit(verbose)

      verbose && exit(verbose)
    } # for (cc ...)

    verbose && exit(verbose)

    NA
  }, chunkSize=unitsPerChunk, useNames=FALSE, verbose=less(verbose, 10)) # lapplyInChunks()

  verbose && exit(verbose)

  res
}) # getAverageFile()


setMethodS3("writeDataFrame", "AromaUnitPscnBinarySet", function(this, filename=sprintf("%s,total.txt", getFullName(this)), ...) {
  NextMethod("writeDataFrame", filename=filename)
})
