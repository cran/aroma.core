## setMethodS3("str", "raster", function(object, ...) {
##   str(unclass(as.array(object), ...));
## })

setConstructorS3("RasterImage", function(z=NULL, dim=NULL, ...) {
  if (!is.null(z)) {
    if (!is.matrix(z) && !is.array(z)) {
      z <- array(z, dim=dim, ...);
    }
    dim <- dim(z);
    ndim <- length(dim);
    if (ndim < 1 || ndim > 4) {
      throw("The number of dimension of argument 'z' is out of range [1,4]: ", ndim);
    }
  }

  extend(BasicObject(z), "RasterImage");
})

setMethodS3("as.character", "RasterImage", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);
  s <- c(s, sprintf("Dimensions: %s", paste(dim(this), collapse="x")));
  s;
})


setMethodS3("writeImage", "RasterImage", function(x, ...) {
  write(x, ...);
}, protected=TRUE)


setMethodS3("save", "RasterImage", function(this, file=NULL, ...) {
  if (!is.null(file)) {
    ext <- gsub("(.*)[.]([^.]*)$", "\\2", file);
    ext <- tolower(ext);
    if (is.element(ext, c("png", "jpg", "jpeg", "gif"))) {
      throw("Cannot save ", class(this)[1], " as an ", toupper(ext), " image file. Did you mean to use write()?");
    }
  }

  NextMethod("save", this, file=file, ...);
})

setMethodS3("write", "RasterImage", function(x, file, path=".", overwrite=FALSE, format=c("auto", "png"), ...) {
  # To please R CMD check
  this <- x;

  # Argument 'format':
  format <- match.arg(format);

  # Infer file format from filename extension?
  if (format == "auto") {
    format <- gsub("(.*)[.]([^.]*)$", "\\2", file);
    format <- tolower(format);

    # Is format supported?
    knownFormats <- as.character(formals(write.RasterImage)$format[-1]);
    format <- match.arg(format, choices=knownFormats);
  }

  # Argument 'file' and 'path':
  pathname <- Arguments$getWritablePathname(file, path=path, mustNotExist=!overwrite, ...);

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose);

  if (format == "png") {
    if (!isPackageInstalled("png")) {
      throw("Cannot write PNG to file. Package not installed: png");
    }
    png::writePNG(this, target=pathnameT);
  }

  # Assert that image file was created
  pathnameT <- Arguments$getReadablePathname(pathnameT, mustExist=TRUE);

  # Rename temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose);

  invisible(pathname);
})

setMethodS3("read", "RasterImage", function(static, ..., format=c("auto", "png")) {
  # Argument 'file' and 'path':
  pathname <- Arguments$getReadablePathname(..., mustExist=TRUE);

  # Argument 'format':
  format <- match.arg(format);

  # Infer file format from filename extension?
  if (format == "auto") {
    format <- gsub("(.*)[.]([^.]*)$", "\\2", pathname);
    format <- tolower(format);

    # Is format supported?
    knownFormats <- formals(read.RasterImage)$format[-1];
    format <- match.arg(format, choices=knownFormats);
  }

  output <- capture.output({
    if (!isPackageInstalled("png")) {
      throw("Cannot read PNG from file. Package not installed: png");
    }
    img <- png::readPNG(pathname, native=FALSE);
  });

  RasterImage(img);
})

setMethodS3("[", "RasterImage", function(...) {
  img <- NextMethod("[");
  RasterImage(img);
}, protected=TRUE)


setMethodS3("dropAlpha", "RasterImage", function(this, ...) {
  dim <- dim(this);
  dim3 <- 1:min(3,dim[3]);
  this[,,dim3, drop=FALSE];
})

setMethodS3("getImageData", "RasterImage", function(this, ...) {
  unclass(this);
}, protected=TRUE)

setMethodS3("setImageData", "RasterImage", function(this, data, ...) {
  RasterImage(data, ...);
}, protected=TRUE)


setMethodS3("rescale", "RasterImage", function(this, scale=1, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'scale':
  scale <- Arguments$getDouble(scale, range=c(0,Inf));
  if (scale == 0) {
    throw("Argument 'scale'is out of range (0,Inf): ", scale);
  }

  # Nothing to do?
  if (scale == 1) {
    return(this);
  }



  z <- getImageData(this);

  # Source image dimensions
  dim0 <- dim <- dim(z);

  # Force a 4-d array while rescaling
  ndim <- length(dim);
  if (ndim < 4) {
    dim <- c(dim, 1L, 1L)[1:4];
    dim(z) <- dim;
  }

  # Rescaled image dimensions
  dimN <- dim;
  dimN[1:2] <- round(scale * dimN[1:2]);

  # Auxillary
  nrow <- dim[1];
  ncol <- dim[2];
  nrow2 <- dimN[1];
  ncol2 <- dimN[2];

  if (nrow2 < 2 || ncol2 < 2) {
    throw("Scale factor is too small: ", scale);
  }


  # Shrink or enlarge?
  if (scale < 1) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # (1) Scale by rows
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Map rows
    destRows <- seq(from=1, to=nrow2, length.out=nrow);
    destRows <- round(destRows);

    # Sanity check
    stopifnot(min(destRows) == 1);
    stopifnot(max(destRows) == nrow2);
    stopifnot(length(unique(destRows)) == nrow2);

    # Allocate nrow2-by-ncol image 
    dim2 <- dim;
    dim2[1] <- nrow2;
    z2 <- array(0, dim=dim2);

    # Sum up the signals
    for (rr in seq(length=nrow)) {
      rr2 <- destRows[rr];
      z2[rr2,,,] <- z2[rr2,,,] + z[rr,,,];
    }

    # Not needed anymore
    rm(z);

    # Rescale intensities
    counts <- as.integer(table(destRows));
    # Sanity check
    stopifnot(length(counts) == nrow2);
    z2 <- z2 / counts;


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # (2) Scale by columns
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Map rows
    destCols <- seq(from=1, to=ncol2, length.out=ncol);
    destCols <- round(destCols);

    # Sanity check
    stopifnot(min(destCols) == 1);
    stopifnot(max(destCols) == ncol2);
    stopifnot(length(unique(destCols)) == ncol2);

    # Allocate nrow2-by-ncol2 image 
    dim3 <- dim2;
    dim3[2] <- ncol2;
    z3 <- array(0, dim=dim3);

    # Sum up the signals
    for (cc2 in seq(length=ncol)) {
      cc3 <- destCols[cc2];
      z3[,cc3,,] <- z3[,cc3,,] + z2[,cc2,,];
    }

    # Not needed anymore
    rm(z2);

    # Rescale intensities
    counts <- as.integer(table(destCols));
    # Sanity check
    stopifnot(length(counts) == ncol2);
    for (cc in seq(length=ncol2)) {
      z3[,cc,,] <- z3[,cc,,] / counts[cc];
    }

    zN <- z3;

    # Not needed anymore
    rm(z3);
  } else if (scale > 1) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # (1) Scale by rows
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Map rows
    srcRows <- seq(from=1, to=nrow, length.out=nrow2);
    srcRows <- round(srcRows);

    # Sanity check
    stopifnot(min(srcRows) == 1);
    stopifnot(max(srcRows) == nrow);
    stopifnot(length(unique(srcRows)) == nrow);

    # Allocate nrow2-by-ncol image 
    dim2 <- dim;
    dim2[1] <- nrow2;
    z2 <- array(0, dim=dim2);

    # Sum up the signals
    for (rr2 in seq(length=nrow2)) {
      rr <- srcRows[rr2];
      z2[rr2,,,] <- z2[rr2,,,] + z[rr,,,];
    }

    # Not needed anymore
    rm(z);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # (2) Scale by columns
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Map columns
    srcCols <- seq(from=1, to=ncol, length.out=ncol2);
    srcCols <- round(srcCols);
    # Sanity check
    stopifnot(min(srcCols) == 1);
    stopifnot(max(srcCols) == ncol);
    stopifnot(length(unique(srcCols)) == ncol);

    # Allocate nrow2-by-ncol2 image 
    dim3 <- dim2;
    dim3[2] <- ncol2;
    z3 <- array(0, dim=dim3);

    # Sum up the signals
    for (cc3 in seq(length=ncol2)) {
      cc2 <- srcCols[cc3];
      z3[,cc3,,] <- z3[,cc3,,] + z2[,cc2,,];
    }

    zN <- z3;

    # Not needed anymore
    rm(z2, z3);
  }

  # Sanity check
  stopifnot(all(dim(zN) == dimN));
  stopifnot(min(zN) >= 0);
  stopifnot(max(zN) <= 1);

  # Preserve the original number of dimensions
  dimN <- dimN[1:ndim];
  dim(zN) <- dimN;

  imgN <- RasterImage(zN, ...);  

  imgN;
})


setMethodS3("str", "RasterImage", function(object, ...) {
  str(unclass(as.array(object), ...));
})


setMethodS3("display", "RasterImage", function(this, dropAlpha=TRUE, interpolate=FALSE, ...) {
  if (!exists("rasterImage", mode="function")) {
    throw("Cannot display ", class(this)[1], ". Only possible in R v2.11.0 or newer: ", as.character(getRversion()));
  }

  dim <- dim(this);
  ndim <- length(dim);

  img <- this;

  if (ndim < 3) {
    img <- colorize(img);
  }

  if (dropAlpha) {
    img <- dropAlpha(img);
  }

  img <- as.raster(img);

  par(mar=c(0,0,0,0), xaxs="i", yaxs="i"); ##, mgp=c(0,0,0), oma=c(0,0,0,0));
  plot(NA, type="n", xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE);
  rasterImage(img, xleft=0, ybottom=0, xright=1, ytop=1, interpolate=interpolate, ...);
})



setMethodS3("colorize", "RasterImage", function(this, palette=gray.colors(256), lim=c(-Inf,Inf), outlierCol="white", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'this':

  # Argument 'lim':

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # Nothing todo?
  if (!is.na(dim(this)[3])) {
    return(this);
  }

  verbose && enter(verbose, "Colorizing image");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Bin the signals
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Binning image signals");
  # Data
  x <- getImageData(this);
  dim <- dim(x);

  # Outliers
  if (lim[1] > -Inf)
    x[x < lim[1]] <- NA;
  if (lim[2] < +Inf)
    x[x > lim[2]] <- NA;
  x[!is.finite(x)] <- NA;

  # Standardize to [0,1]
  r <- range(x, na.rm=TRUE);
  verbose && cat(verbose, "Before:");
  verbose && print(verbose, summary(as.vector(x)));
  x <- (x - r[1])/(r[2]-r[1]);
  verbose && cat(verbose, "After:");
  verbose && print(verbose, summary(as.vector(x)));

  # Standardize to [0,n] where n is number of bins
  n <- length(palette);
  x <- x*n;

  # Bin to [1,n] by moving any values at the left extreme to bin one.
  x <- as.integer(x);
  x[x == 0L] <- 1L;

  verbose && str(verbose, x);
  verbose && summary(verbose, x);
  verbose && exit(verbose);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Assign colors to each bin
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  verbose && enter(verbose, "Binning image data to color map");
  # Translate 'palette' of length n into a 4xn (R,G,B,A) matrix
  paletteRGBA <- col2rgb(palette, alpha=TRUE);
  colnames(paletteRGBA) <- palette;
  verbose && cat(verbose, "Color map:");
  verbose && str(verbose, paletteRGBA);

  naRGBA <- col2rgb(outlierCol, alpha=TRUE);
  colnames(naRGBA) <- outlierCol;
  verbose && cat(verbose, "Missing-value color:");
  verbose && str(verbose, naRGBA);

  # Map image to colors
  z <- paletteRGBA[,x,drop=FALSE];
  z[,is.na(x)] <- naRGBA;
  dim(z) <- c(4L, dim);

  # Rescale to [0,1]
  z <- z / 255;

  verbose && str(verbose, z);
  verbose && summary(verbose, z);

  # Sanity checks
  stopifnot(min(z, na.rm=TRUE) >= 0);
  stopifnot(max(z, na.rm=TRUE) <= 1);

  verbose && exit(verbose);


  verbose && enter(verbose, "Restructuring");

  z <- aperm(z, perm=c(2,3,1));
  verbose && str(verbose, z);
  verbose && summary(verbose, z);

  # Is the alpha channel needed?
  if (all(z[,,4,drop=TRUE] == 1)) {
    verbose && cat(verbose, "Dropped the alpha channel, because it is not needed.");
    z <- z[,,1:3,drop=FALSE];
    verbose && str(verbose, z);
    verbose && summary(verbose, z);
  }

  verbose && exit(verbose);


  # 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Create a new Image object
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  img <- RasterImage(z, ...);

  verbose && exit(verbose);

  img;
}, protected=TRUE)



############################################################################
# HISTORY:
# 2011-02-24
# o Now write() and read() for RasterImage throws an informative error
#   message explaining that the 'png' package is needed.
# 2011-01-31
# o Added interleave() for RasterImage.
# o Added colorize() for RasterImage.
# o Now write() is atomic by writing to a temporary file.
# o Added "smart" save() for RasterImage that checks if user really
#   mean to call write() instead.
# o Added "smart" read() and write() for RasterImage.
# o Added rescale() for RasterImage.  Works both for shrinking and 
#   enlarging images.
# o Added RasterImage class for supporting 'raster' images.
# o Created.
############################################################################