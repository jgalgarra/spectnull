debugpref <- "DEB"
ppi <- 100
options( warn = -1 )
datadir <- "data/"
dataprocessed <- "dataprocessed/"
lweightrf <- c("none")#,"sqrt","ln")
dbaseb <- "smodels/"
odirb <- "plots/"
rdirb <- "results/"
dirnullsb <- "nullmatrix/"
# Configuration parameters
seed <- 122
num_experiments <- 5
ignore_GC_results <- FALSE # Ignore the results of Giant Component files
plottofile <- TRUE # Save individual network distributions plot
plotzigs <- FALSE  # Plotting ziggurats of all models is rather slow. So when TRUE magnitudes are
# not saved. Run the script with a big number of experiments (~1000) to compute
# magnitudes and plotzigs FALSE. Run it again selecting just the networks you need
# to plot and a small number of experiments (~10)
NetworkMagsFile <- "NetworkMagnitudes.csv" # Stores network magnitudes and average spectral measures
nmagnitudes <- list("spect_rad","adj_energy","lpl_spect_rad","lpl_energy")#,"algebraic_connectivity")
# Null models for binary/binarized and weighted networks
mnamesbin <- c("RND","MGEN","SHUFFLE","VAZ","SYTR")
mnamesweighted <- c("SWAP","WRND","BVAZ","BSHUFFLE","PATEFIELD")
MIN_LINKS_SIZE <- 20  # Smaller networks are discarded


create_dirs <- function(weightrf){
  dbase <<- paste0(debugpref,dbaseb)
  dir.create(dbase, showWarnings = FALSE)
  dbase <<- paste0(dbase,weightrf,"/")
  dir.create(dbase, showWarnings = FALSE)
  odir <<- paste0(dbase,odirb)
  rdir <<- paste0(dbase,rdirb)
  dirnulls <<- paste0(dbase,dirnullsb)
  dir.create(odir, showWarnings = FALSE)
  dir.create(rdir, showWarnings = FALSE)
  dir.create(dirnulls, showWarnings = FALSE)
}
