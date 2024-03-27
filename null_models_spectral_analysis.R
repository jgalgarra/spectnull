rm(list = ls())
library(kcorebip)
library(ggplot2)
library(patchwork)
library(forcats)
source("synthrade_nullmodel.R")
source("colormodels.R")
source("config_values.R")
source("plot_null_distributions.R")
source("manage_null_models.R")

INCOMPLETE_MODEL_VALUE <<- -1

network_null_spectral_distances <- function(netw,weightrf,numexperiments,mnamesbin,mnamesweighted,datadir,rdir,odir,dirnulls,
                                            guild_a_label = "Plant", guild_b_label = "Pollinator",min_links=20,
                                            plottofile=TRUE,plotzigs=FALSE,networkmagsfile="NetworkMagnitudes.csv")
{
  create_dirs(weightrf)
  print(paste("Network:",netw,"transform",weightrf))
  result_analysis <- analyze_network(directory = datadir, netw, guild_a = guild_a_label, guild_b = guild_b_label, only_NODF = TRUE)
  if (result_analysis$links <= min_links)      # Discard tiny networks
    return()
  # Full network process
  nnm <- process_network_null_models(netw,result_analysis,num_experiments,mnamesbin,
                                     mnamesweighted,weightrf)
  # Stack null model results: value, magnitude, model
  datamod <- stack_models_results(nnm)
  # Compute mean, quantile and raw spectral distance
  compute_statistic_indexes(netw,nnm$weighted_network,datamod,nnm$networkspect,
                            nnm$modresults,nnm$mnames)
  # Save model magnitudes and statistic files
  write.csv(datamod,paste0(rdir,"MODS_",netw),row.names = FALSE)  # Model results
  write.csv(nnm$dfnested,paste0(rdir,"NESTED_",netw),row.names = FALSE)  # Model results
  # Compute spectral magnitudes of null HNESTED and WNESTED models
  store_network_magnitudes_global_file(netw,nnm)
  if (plottofile)
    plot_all_distr(gsub(".csv","",netw),plotzigs,nnm)
  if (plotzigs)
    plot_ziggurats(nnm,netw,dirnulls,datadir,odir)
}

######################################################################################################
# 
# PROGRAM BODY
#
#
# Configuration parameters
seed <- 122
num_experiments <- 5
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

# Here, the list of data files to process
filenames <- Sys.glob(paste0(datadir,"*PL*01?.csv"))
# Network names
lnetw <- gsub(datadir,"",filenames)
for (netw in (lnetw)){                 # Each network
  for (weightrf in lweightrf)
    network_null_spectral_distances(netw,weightrf,numexperiments,mnamesbin,mnamesweighted,datadir,rdir,odir,dirnulls,
                                              guild_a_label = "Plant", guild_b_label = "Pollinator",min_links=20,
                                              plottofile=TRUE,plotzigs=plotzigs,networkmagsfile="NetworkMagnitudes.csv")
}