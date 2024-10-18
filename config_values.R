debugpref <- "PRU"
ppi <- 100
options( warn = -1 )
datadir <- "data/"
dataprocessed <- "dataprocessed/"
lweightrf <- c("none")#,"sqrt","ln")
dbaseb <- "smodels/"
odirb <- "plots/"
rdirb <- "results/"
dirnullsb <- "nullmatrix/"
eigenplots <- "eigenplots/"
matrixplots <- "matrixplots/"
# Configuration parameters
MIN_LINKS <- 20  # Smaller networks are discarded
seed <- 122
num_experiments <- 10
ignore_GC_results <- FALSE # Ignore the results of Giant Component files
plottofile <- TRUE # Save individual network distributions plot
plotzigs <- FALSE  # Plotting ziggurats of all models is rather slow. So when TRUE magnitudes are
# not saved. 
cold_start <- TRUE # If TRUE removes the output directory tree
NetworkMagsFile <- "NetworkMagnitudes.csv" # Stores network magnitudes and average spectral measures
nmagnitudes <- list("spect_rad","adj_energy","lpl_spect_rad","lpl_energy")#,"algebraic_connectivity")
# Null models for binary/binarized and weighted networks
mnamesbin <- c("RND","MGEN","SHUFFLE","VAZ","SYTR")
mnamesweighted <- c("SWAP","WRND","BVAZ","BSHUFFLE","PATEFIELD")

create_dirs <- function(weightrf){
  dbase <<- paste0(debugpref,dbaseb)
  dir.create(dbase, showWarnings = FALSE)
  dbase <<- paste0(dbase,weightrf,"/")
  dir.create(dbase, showWarnings = FALSE)
  odir <<- paste0(dbase,odirb)
  rdir <<- paste0(dbase,rdirb)
  eiplotsdir <<- paste0(dbase,eigenplots)
  matrixplotsdir <<- paste0(dbase,matrixplots)
  dirnulls <<- paste0(dbase,dirnullsb)
  dir.create(odir, showWarnings = FALSE)
  dir.create(rdir, showWarnings = FALSE)
  dir.create(dirnulls, showWarnings = FALSE)
  dir.create(eiplotsdir, showWarnings = FALSE)
  dir.create(matrixplotsdir, showWarnings = FALSE)
}

minmaxnorm <- function(val,minv,maxv){
  return((val-minv)/(maxv-minv))
}
