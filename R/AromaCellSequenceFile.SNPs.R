setMethodS3("getSnpPositions", "AromaCellSequenceFile", function(this, cells, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (!is.matrix(cells)) {
    throw("Argument 'cells' must be a matrix: ", class(cells)[1])
  }
  dim <- dim(cells)
  if (!any(dim == 2)) {
    throw("Argument 'cells' must have either two rows or two columns: ", paste(dim, collapse="x"))
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Identifying SNP positions of cell allele pairs")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this)
  key <- list(method="getSnpPositions", class=class(this)[1],
              version="20110830",
              chipType=chipType, tags=getTags(this), cells=cells, ...)
  if (getOption(aromaSettings, "devel/useCacheKeyInterface", FALSE)) {
    key <- getCacheKey(this, method="getSnpPositions", chipType=chipType,
                                     tags=getTags(this), cells=cells, ...)
  }
  dirs <- c("aroma.affymetrix", chipType)
  if (!force) {
    verbose && enter(verbose, "Checking for cached results")
    res <- loadCache(key=key, dirs=dirs)
    if (!is.null(res)) {
      verbose && cat(verbose, "Found cached results")
      verbose && exit(verbose)
      verbose && exit(verbose)
      return(res)
    }
    verbose && exit(verbose)
  }

  byRow <- (dim[2] == 2)
  if (byRow) {
    cells <- t(cells)
    dim <- dim(cells)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read probe sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Reading sequences matrix")
  seqs <- readSequenceMatrix(this, cells=cells, what="raw")
  map <- attr(seqs, "map")
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify SNP positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Searching for mismatch position in pairs")
  probeLength <- dim(seqs)[2]
  dim(seqs) <- c(dim, probeLength)
  seqsA <- seqs[1,,,drop=FALSE]
  seqsB <- seqs[2,,,drop=FALSE]
  dim(seqsA) <- dim(seqsB) <- c(dim[2], probeLength)
  # Not needed anymore
  seqs <- NULL

  # Identify the *last* difference
  naValue <- NA_integer_
  pos <- rep(naValue, times=ncol(cells))
  for (pp in seq_len(probeLength)) {
    idxs <- which(seqsA[,pp] != seqsB[,pp])
    pos[idxs] <- pp
  }
  # Not needed anymore
  idxs <- pp <- NULL

  # Sanity check
  .stop_if_not(length(pos) == ncol(cells))

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Caching result")
  saveCache(pos, key=key, dirs=dirs)
  verbose && exit(verbose)

  verbose && exit(verbose)

  pos
}) # getSnpPositions()


setMethodS3("getSnpShifts", "AromaCellSequenceFile", function(this, ...) {
  pos <- getSnpPositions(this, ...)
  pos <- pos - ((getProbeLength(this) %/% 2) + 1L)
  pos
}) # getSnpShifts()



setMethodS3("getSnpNucleotides", "AromaCellSequenceFile", function(this, cells, shifts=0, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'cells':
  if (!is.matrix(cells)) {
    throw("Argument 'cells' must be a matrix: ", class(cells)[1])
  }
  dim <- dim(cells)
  if (!any(dim == 2)) {
    throw("Argument 'cells' must have either two rows or two columns: ", paste(dim, collapse="x"))
  }

  # Argument 'shifts':
  shifts <- Arguments$getIntegers(shifts, range=c(-5,5))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Retrieving SNP nucleotides")
  verbose && cat(verbose, "Probe shifts:")
  verbose && str(verbose, shifts)

  verbose && cat(verbose, "Argument 'cells':")
  verbose && str(verbose, cells)

  byRow <- (dim[2] == 2)
  if (byRow) {
    cells <- t(cells)
    dim <- dim(cells)
  }

  verbose && cat(verbose, "Number of pairs: ", ncol(cells))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify SNP positions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pos <- getSnpPositions(this, cells=cells, verbose=less(verbose, 1))
  verbose && cat(verbose, "SNP positions:")
  verbose && str(verbose, pos)

  # Sanity check
  .stop_if_not(length(pos) == ncol(cells))

  verbose && cat(verbose, "Tabulated SNP positions:")
  verbose && print(verbose, table(pos))

  uPos <- sort(unique(pos))
  verbose && cat(verbose, "Unique SNP positions:")
  verbose && print(verbose, uPos)

  # Sanity check
  .stop_if_not(length(uPos) > 0)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Read probe sequences
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  naValue <- NA_character_
  seqs <- rep(naValue, times=length(cells))
  dim(seqs) <- dim(cells)
  for (pp in seq_along(uPos)) {
    snpPosition <- uPos[pp]

    # Cells to read
    idxs <- which(pos == snpPosition)
    cellsPP <- cells[,idxs,drop=FALSE]

    # Sequence positions to read
    positions <- snpPosition + shifts

    seqsPP <- readSequences(this, cells=cellsPP, positions=positions)
    seqs[,idxs] <- seqsPP
  }
  # Not needed anymore
  idxs <- seqsPP <- positions <- cellsPP <- snpPosition <- cells <- pos <- NULL

  if (byRow) {
    seqs <- t(seqs)
  }

  verbose && cat(verbose, "Probe sequences:")
  verbose && str(verbose, seqs)
  verbose && exit(verbose)

  seqs
}) # getSnpNucleotides()



setMethodS3("groupBySnpNucleotides", "AromaCellSequenceFile", function(this, cells, ignoreOrder=TRUE, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ignoreOrder':
  ignoreOrder <- Arguments$getLogical(ignoreOrder)

  # Argument 'force':
  force <- Arguments$getLogical(force)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Identifying groups of SNP nucleotide sequence pairs")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  chipType <- getChipType(this)
  key <- list(method="groupBySnpNucleotides", class=class(this)[1],
              chipType=chipType, tags=getTags(this),
              cells=cells, ignoreOrder=ignoreOrder,
              version="2008-12-04", ...)
  dirs <- c("aroma.affymetrix", chipType)
  if (!force) {
    verbose && enter(verbose, "Checking for cached results")
    res <- loadCache(key=key, dirs=dirs)
    if (!is.null(res)) {
      verbose && cat(verbose, "Found cached results")
      verbose && exit(verbose)
      verbose && exit(verbose)
      return(res)
    }
    verbose && exit(verbose)
  }


  verbose && enter(verbose, "Get SNP nucleotide sequence pairs")
  pairs <- getSnpNucleotides(this, cells=cells, ..., verbose=verbose)
  verbose && str(verbose, pairs)
  verbose && exit(verbose)

  dim <- dim(pairs)

  byRow <- (dim[2] == 2)
  if (byRow) {
    verbose && enter(verbose, "Transposing matrices")
    cells <- t(cells)
    pairs <- t(pairs)
    dim <- dim(cells)
    verbose && exit(verbose)
  }

  verbose && enter(verbose, "Generating names of pairs")
  pairNames <- paste(pairs[1,], pairs[2,], sep="/")
  pairNames[is.na(pairs[1,])] <- NA
  uniquePairs <- sort(unique(pairNames))
  verbose && cat(verbose, "uniquePairs:")
  verbose && str(verbose, uniquePairs)
  verbose && exit(verbose)

  verbose && enter(verbose, "Identifying names of pairs to map to")
  pairsToBuild <- uniquePairs
  if (ignoreOrder) {
    pairsToBuild <- strsplit(pairsToBuild, split="/", fixed=TRUE)
    pairsToBuild <- lapply(pairsToBuild, FUN=sort)
    pairsToBuild <- sapply(pairsToBuild, FUN=paste, collapse="/")
    pairsToBuild <- unique(pairsToBuild)
    pairsToBuild <- sort(pairsToBuild)
  }
  verbose && cat(verbose, "pairsToBuild:")
  verbose && str(verbose, pairsToBuild)
  verbose && exit(verbose)

  res <- vector("list", length(pairsToBuild)+1)
  names(res) <- c(pairsToBuild, "missing")
  for (kk in seq_along(uniquePairs)) {
    pair <- uniquePairs[kk]
    verbose && enter(verbose, sprintf("Pair %d ('%s') of %d",
                                          kk, pair, length(uniquePairs)))

    idxs <- which(pairNames == pair)
    verbose && cat(verbose, "Matching pairs:")
    verbose && str(verbose, idxs)
    cellsKK <- cells[,idxs,drop=FALSE]

    # Swap?
    if (ignoreOrder) {
      verbose && enter(verbose, "Order pair")
      pair <- strsplit(pair, split="/", fixed=TRUE)[[1]]
      if (pair[1] > pair[2]) {
        cellsKK <- cellsKK[2:1,,drop=FALSE]
        pair <- rev(pair)
      }
      pair <- paste(pair, collapse="/")
      verbose && exit(verbose)
    }

    res[[pair]] <- cbind(res[[pair]], cellsKK)

    verbose && exit(verbose)
  } # for (kk ...)

  idxs <- which(is.na(pairNames))
  cellsKK <- cells[,idxs,drop=FALSE]
  if (length(cellsKK) > 0)
    res[["missing"]] <- cellsKK
  # Not needed anymore
  idxs <- cellsKK <- NULL

  for (kk in seq_along(res)) {
    cells <- res[[kk]]
    pair <- names(res)[kk]
    if (pair == "missing") {
      dimnames(cells) <- NULL
    } else {
      pair <- strsplit(pair, split="/", fixed=TRUE)[[1]]
      rownames(cells) <- pair
    }
    if (byRow)
      cells <- t(cells)
    res[[kk]] <- cells
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Save to cache
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Caching result")
  saveCache(res, key=key, dirs=dirs)
  verbose && exit(verbose)

  verbose && exit(verbose)

  res
}) # groupBySnpNucleotides()
