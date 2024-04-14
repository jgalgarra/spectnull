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
                                            guild_a_label = "Plant", guild_b_label = "Pollinator",min_links=MIN_LINKS,
                                            plottofile=TRUE,plotzigs=FALSE,networkmagsfile=NetworkMagsFile)
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

dbase <- paste0(debugpref,dbaseb)
if (cold_start){
  inp <- readline(paste("The script is configured to remove all results stored at",dbase,". Are you sure (Y/N)? "))
  inp <- as.character(inp)
  if (inp=="Y")
    unlink(dbase, recursive=TRUE)
  else
    stop("OK, cancelled mission!")
}
  # Here, the list of data files to process
filenames <- Sys.glob(paste0(datadir,"*_*.csv"))
# Network names
lnetw <- gsub(datadir,"",filenames)
for (netw in (lnetw)){                 # Each nesdtwork
  for (weightrf in lweightrf)
    network_null_spectral_distances(netw,weightrf,numexperiments,mnamesbin,mnamesweighted,datadir,rdir,odir,dirnulls,
                                              guild_a_label = "Plant", guild_b_label = "Pollinator",min_links=MIN_LINKS,
                                              plottofile=TRUE,plotzigs=plotzigs,networkmagsfile=NetworkMagsFile)
}