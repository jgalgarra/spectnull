# Script to create a unique file of normalized distances ALLNORMALIZED.cand
# an auxiliary file with the types of networks networtypes.csv
# Run it after bipartite_graph_spectrum_NULL_random.R

library("data.table")
source("config_values.R")

minmaxnorm <- function(val,minv,maxv){
  return((val-minv)/(maxv-minv))
  }


find_reference <- function(values_energy,mag,model,stat="median"){
  d <- values_energy[(values_energy$ind==mag) & (values_energy$MODEL==model),]$values
  if (stat=="median")
    return(median(d))
  else
    return(mean(d))
}

computemeandist <- function(alldistances,magnitude,netw){
  magnitude_data <- alldistances[(alldistances$ind==magnitude) & (alldistances$network==netw),]
  models <- unique(magnitude_data$MODEL)
  for (m in models){
    #if (m!="NETWORK")
    meandist <- rbind(meandist,data.frame("network"=netw,
                                          "meannormdist"=mean(magnitude_data[magnitude_data$MODEL==m,]$normdist),
                                          "MODEL"=m,
                                          "disttype"=magnitude,
                                          "networktype"=magnitude_data$type[1]))
  }
  return(meandist)
}
for (weightrf in lweightrf)
{
  print(paste("Transformation",weightrf))
  rdir <- paste0(debugpref,dbaseb,weightrf,"/",rdirb)
  resultsdir <- rdir
  normresults <- paste0(rdir,"normalized/")
  dir.create(normresults, showWarnings = FALSE)
  filenames <- Sys.glob(paste0(resultsdir,"MODS_*.csv"))
  lnetw <- gsub(".csv","",gsub(resultsdir,"",filenames))
  networktypes <- data.frame("network"=c(),"type"=c())
  for (netw in lnetw)
  {
    nname <- gsub("MODS_","",netw)
    nmatrix <- fread(paste0(datadir,nname,".csv"))
  
    if (sum(nmatrix[2:nrow(nmatrix),2:ncol(nmatrix)]>1)>0)
      matrixtype = "WEIGHTED"
    else
      matrixtype = "BINARY"
    networktypes <- rbind(networktypes,data.frame("network"=nname,"type"=matrixtype))
    print(netw)
    values_network <- fread(paste0(resultsdir,netw,".csv"))
    if (matrixtype=="WEIGHTED")
      values_energy <- values_network[values_network$ind %in% c("adj_energy","lpl_energy","spect_rad","lpl_weighted_energy","adj_weighted_energy","spect_rad_weighted","algebraic_connectivity"),]
    else
      values_energy <- values_network[values_network$ind %in% c("adj_energy","lpl_energy","spect_rad"),]
    models <- unique(values_energy$MODEL)
    adj_energy_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="adj_energy"),]$values
    lpl_energy_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="lpl_energy"),]$values
    spect_rad_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="spect_rad"),]$values
    if (matrixtype=="WEIGHTED"){
      adj_weighted_energy_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="adj_weighted_energy"),]$values
      lpl_weighted_energy_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="lpl_weighted_energy"),]$values
      spect_rad_weighted_empirical <- values_energy[(values_energy$MODEL=="NETWORK") & (values_energy$ind=="spect_rad_weighted"),]$values
    }
  
    min_adj_energy <- find_reference(values_energy,"adj_energy","HNESTED")
    max_adj_energy <- find_reference(values_energy,"adj_energy","RND")
    min_spect_rad <- find_reference(values_energy,"spect_rad","RND")
    max_spect_rad <- find_reference(values_energy,"spect_rad","HNESTED")
    min_lpl_energy <- find_reference(values_energy,"lpl_energy","RND")
    max_lpl_energy <- find_reference(values_energy,"lpl_energy","HNESTED")
  
    if (matrixtype=="WEIGHTED"){
      min_adj_weighted_energy <- find_reference(values_energy,"adj_weighted_energy","WRND")
      max_adj_weighted_energy <- find_reference(values_energy,"adj_weighted_energy","WNESTED")
      min_spect_rad_weighted <- find_reference(values_energy,"spect_rad_weighted","WRND")
      max_spect_rad_weighted <- find_reference(values_energy,"spect_rad_weighted","WNESTED")
      min_lpl_weighted_energy <- find_reference(values_energy,"lpl_weighted_energy","WRND")
      max_lpl_weighted_energy <- find_reference(values_energy,"lpl_weighted_energy","WNESTED")
    }
  
  
    adj_energy_empirical_norm <- minmaxnorm(adj_energy_empirical,min_adj_energy,max_adj_energy)
    lpl_energy_empirical_norm <- minmaxnorm(lpl_energy_empirical,min_lpl_energy,max_lpl_energy)
    spect_rad_empirical_norm <- minmaxnorm(spect_rad_empirical,min_spect_rad,max_spect_rad)
    if (matrixtype=="WEIGHTED"){
      adj_weighted_energy_empirical_norm <- minmaxnorm(adj_weighted_energy_empirical,min_adj_weighted_energy,max_adj_weighted_energy)
      lpl_weighted_energy_empirical_norm <- minmaxnorm(lpl_weighted_energy_empirical,min_lpl_weighted_energy,max_lpl_weighted_energy)
      spect_rad_weighted_empirical_norm <- minmaxnorm(spect_rad_weighted_empirical,min_spect_rad_weighted,max_spect_rad_weighted)
    }
  
    norm_values_energy <- values_energy
    norm_values_energy$normdist <- 0
    for (i in 1:nrow(norm_values_energy)){
      if (norm_values_energy$ind[i]=="algebraic_connectivity")
        norm_values_energy$normdist[i] = norm_values_energy$values[i]
      if (norm_values_energy$ind[i]=="adj_energy")
        norm_values_energy$normdist[i] = minmaxnorm(norm_values_energy$values[i],min_adj_energy,max_adj_energy)
      if (norm_values_energy$ind[i]=="lpl_energy")
        norm_values_energy$normdist[i] = minmaxnorm(norm_values_energy$values[i],min_lpl_energy,max_lpl_energy)
      if (norm_values_energy$ind[i]=="spect_rad")
        norm_values_energy$normdist[i] = minmaxnorm(norm_values_energy$values[i],min_spect_rad,max_spect_rad)
      if (matrixtype=="WEIGHTED"){
        if (norm_values_energy$ind[i]=="adj_weighted_energy")
          norm_values_energy$normdist[i] = minmaxnorm(norm_values_energy$values[i],min_adj_weighted_energy,max_adj_weighted_energy)
        if (norm_values_energy$ind[i]=="lpl_weighted_energy")
          norm_values_energy$normdist[i] = minmaxnorm(norm_values_energy$values[i],min_lpl_weighted_energy,max_lpl_weighted_energy)
        if (norm_values_energy$ind[i]=="spect_rad_weighted")
          norm_values_energy$normdist[i] = minmaxnorm(norm_values_energy$values[i],min_spect_rad_weighted,max_spect_rad_weighted)
      }
    }
    fwrite(norm_values_energy,paste0(normresults,"MINMAX_",netw,".csv"),row.names=FALSE)
  }
  fwrite(networktypes,paste0(normresults,"networktypes.csv"),row.names=FALSE)
  
  print("Merging all normalized distance files. It may take a while, be patient please.")
  filenames <- Sys.glob(paste0(normresults,"MINMAX_MODS_*.csv"))
  alldistances <- fread(filenames[1])
  nname <- gsub(".csv","",gsub(normresults,"",filenames[1]))
  alldistances$network <- nname
  if (length(filenames)>1)
    for (i in 2:length(filenames)){
      ndistances <- fread(filenames[i])
      nname <- gsub(".csv","",gsub(normresults,"",filenames[i]))
      ndistances$network <- nname
      alldistances <- rbind(alldistances,ndistances)
      alldistances <- alldistances[!is.na(alldistances$network),]
    }
  fwrite(alldistances,paste0(normresults,"ALLNORMALIZED.csv"),row.names=FALSE)
  
  # Merging nestedness measures
  
  print("Merging nestedness measures.")
  filenames <- Sys.glob(paste0(resultsdir,"NESTED_*.csv"))
  nestedmeasures <- fread(filenames[1])
  nname <- gsub("NESTED_","",gsub(".csv","",gsub(resultsdir,"",filenames[1])))
  nname <- gsub(".csv","",nname)
  nestedmeasures$network <- nname
  if (length(filenames)>1)
    for (i in 2:length(filenames)){
      ndistances <- fread(filenames[i])
      nname <- gsub("NESTED_","",gsub(".csv","",gsub(resultsdir,"",filenames[i])))
      ndistances$network <- nname
      nestedmeasures <- rbind(nestedmeasures,ndistances)
      nestedmeasures <- nestedmeasures[!is.na(nestedmeasures$network),]
    }
  nestedmeasures <-nestedmeasures[,-c("links","totalweight")]
  fwrite(nestedmeasures,paste0(normresults,"ALLNESTEDMEASURES.csv"),row.names=FALSE)
}
