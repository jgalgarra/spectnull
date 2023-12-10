# Script to create a unique file of normalized distances ALLNORMALIZED.cand
# an auxiliary file with the types of networks networtypes.csv
# Run it after bipartite_graph_spectrum_NULL_random.R

minmaxnorm <- function(val,minv,maxv){
  (val-minv)/(maxv-minv)
}

resultsdir <- "results/"
datadir <- "data/"
normresults <- "results/normalized/"
dir.create(normresults, showWarnings = FALSE)
filenames <- Sys.glob(paste0(resultsdir,"MODS_*.csv"))
lnetw <- gsub(".csv","",gsub(resultsdir,"",filenames))
#lnetw <- c("MODS_M_PL_065")
networktypes <- data.frame("network"=c(),"type"=c())
for (netw in lnetw)
{
  nname <- gsub("MODS_","",netw)
  nmatrix <- read.csv(paste0(datadir,nname,".csv"))
  if (sum(nmatrix[2:nrow(nmatrix),2:ncol(nmatrix)]>1)>0)
    matrixtype = "WEIGHTED"
  else
    matrixtype = "BINARY"
  networktypes <- rbind(networktypes,data.frame("network"=nname,"type"=matrixtype))
  print(netw)
  values_network <- read.csv(paste0(resultsdir,netw,".csv"))
  values_energy <- values_network[values_network$ind %in% c("adj_energy","lpl_energy","spect_rad"),]
  models <- unique(values_energy$MODEL)
  adj_energy_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="adj_energy"),]$values
  lpl_energy_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="lpl_energy"),]$values
  spect_rad_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="spect_rad"),]$values
  
  # min_adj_energy <- min(values_energy[(values_energy$ind=="adj_energy"),]$values)
  # max_adj_energy <- max(values_energy[(values_energy$ind=="adj_energy"),]$values)
  # min_lpl_energy <- min(values_energy[(values_energy$ind=="lpl_energy"),]$values)
  # max_lpl_energy <- max(values_energy[(values_energy$ind=="lpl_energy"),]$values)
  # min_spect_rad <- min(values_energy[(values_energy$ind=="spect_rad"),]$values)
  # max_spect_rad <- max(values_energy[(values_energy$ind=="spect_rad"),]$values)
  quantmin = 0.05
  quantmax = 0.95
  min_adj_energy <- quantile(values_energy[(values_energy$ind=="adj_energy"),]$values,c(quantmin))
  max_adj_energy <- quantile(values_energy[(values_energy$ind=="adj_energy"),]$values,c(quantmax))
  min_lpl_energy <- quantile(values_energy[(values_energy$ind=="lpl_energy"),]$values,c(quantmin))
  max_lpl_energy <- quantile(values_energy[(values_energy$ind=="lpl_energy"),]$values,c(quantmax))
  min_spect_rad <- quantile(values_energy[(values_energy$ind=="spect_rad"),]$values,c(quantmin))
  max_spect_rad <- quantile(values_energy[(values_energy$ind=="spect_rad"),]$values,c(quantmax))
  
  adj_energy_empirical_norm <- minmaxnorm(adj_energy_empirical,min_adj_energy,max_adj_energy)
  lpl_energy_empirical_norm <- minmaxnorm(lpl_energy_empirical,min_lpl_energy,max_lpl_energy)
  spect_rad_empirical_norm <- minmaxnorm(spect_rad_empirical,min_spect_rad,max_spect_rad)
  
  
  norm_values_energy <- values_energy
  norm_values_energy$normdist <- 0
  for (i in 1:nrow(norm_values_energy)){
    if (norm_values_energy$ind[i]=="adj_energy"){
      norm_values_energy$values[i] = minmaxnorm(norm_values_energy$values[i],min_adj_energy,max_adj_energy)
      norm_values_energy$normdist[i] <- norm_values_energy$values[i] - adj_energy_empirical_norm
    }
    if (norm_values_energy$ind[i]=="lpl_energy"){
      norm_values_energy$values[i] = minmaxnorm(norm_values_energy$values[i],min_lpl_energy,max_lpl_energy)
      norm_values_energy$normdist[i] <- norm_values_energy$values[i] - lpl_energy_empirical_norm
    }
    if (norm_values_energy$ind[i]=="spect_rad"){
      norm_values_energy$values[i] = minmaxnorm(norm_values_energy$values[i],min_spect_rad,max_spect_rad)
      norm_values_energy$normdist[i] <- norm_values_energy$values[i] - spect_rad_empirical_norm
    }
  }
  norm_values_energy[norm_values_energy$ind == "adj_energy",]$ind = "adj_norm_energy"
  norm_values_energy[norm_values_energy$ind == "lpl_energy",]$ind = "lpl_norm_energy"
  norm_values_energy[norm_values_energy$ind == "spect_rad",]$ind = "spect_rad_norm"
  write.csv(norm_values_energy,paste0(normresults,"MINMAX_",netw,".csv"),row.names=FALSE)
}
write.csv(networktypes,paste0(normresults,"networktypes.csv"),row.names=FALSE)

print("Merging all normalized distance files. It may take a while, be patient please.")
filenames <- Sys.glob(paste0(normresults,"MINMAX_MODS_*.csv"))
alldistances <- read.csv(filenames[1])
nname <- gsub(".csv","",gsub(normresults,"",filenames[1]))
alldistances$network <- nname
for (i in 2:length(filenames)){
  ndistances <- read.csv(filenames[i])
  nname <- gsub(".csv","",gsub(normresults,"",filenames[i]))
  ndistances$network <- nname
  alldistances <- rbind(alldistances,ndistances)
  alldistances <- alldistances[!is.na(alldistances$network),]
}
write.csv(alldistances,paste0(normresults,"ALLNORMALIZED.csv"),row.names=FALSE)