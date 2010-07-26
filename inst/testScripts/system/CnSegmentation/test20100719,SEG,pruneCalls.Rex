if (interactive()) savehistory();
library("aroma.core");

verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Simulating copy-number data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cnr <- function(start, stop, mean=1, chromosome=1, xScale=1e6, ...) {
  CopyNumberRegions(chromosome=chromosome, 
                    start=xScale*start, stop=xScale*stop, mean=mean);
} 

# Build up CN profile generation by generation
t <- cnr(1,2000, 2) +
     cnr(1000,1500) + 
     cnr(1000,1250) +
     cnr(1650,1800) + 
     cnr(200,300) - cnr(650,800);
print(t);

# Simulate data from the track
cn <- simulateRawCopyNumbers(t, n=2000, sd=1/2);
print(cn);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Draw
np <- 8;
subplots(np, ncol=1);
par(mar=c(2,4,1,1)+0.1);

plot(cn, col="#aaaaaa", ylim=c(0,5));
drawLevels(t, col="orange", lwd=4); 

res <- segmentByPruneCBS(cn, normalMean=2, strict=TRUE, debugPlot=TRUE, ylim=c(0,5), verbose=-10);
print(res);
