# This script computes the average normalized distances.
# It must be run after normalize_distances.R

library("data.table")
source("config_values.R")


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
  print("Reading all normalized distances from file")
  alldistances <- fread(paste0(normresults,"ALLNORMALIZED.csv"))
  alldistances$network <- gsub("MINMAX_MODS_","",alldistances$network)
  networktypes <- read.csv(paste0(normresults,"networktypes.csv"))
  
  print("Computing averages")
  lnetw <- unique(alldistances$network)
  
  #lnetw <- lnetw[!grepl("HP",lnetw)]
  
  alldistances$type <- "BINARY"
  meandist <- data.frame("network"=c(),"meannormdist"=c(),"MODEL"=c(), "disttype"=c(), "networktype"=c())
  for (netw in lnetw){
    ntype <- networktypes[networktypes$network==netw,]
    nestedmags <- read.csv(paste0(resultsdir,"NESTED_",netw,".csv"))
    alldistances[alldistances$network==netw,]$type <- ntype$type
    print(netw)
    mdistalgcon <- computemeandist(alldistances,"algebraic_connectivity",netw)
    mdistadj <- computemeandist(alldistances,"adj_energy",netw)
    mdistlpl <- computemeandist(alldistances,"lpl_energy",netw)
    mdistspect <- computemeandist(alldistances,"spect_rad",netw)
    mdistlplspect <- computemeandist(alldistances,"lpl_spect_rad",netw)
    meandist <- rbind(meandist,mdistalgcon,mdistadj,mdistlpl,mdistspect,mdistlplspect)
    if (ntype$type=="WEIGHTED"){
      mdistadj <- computemeandist(alldistances,"adj_weighted_energy",netw)
      mdistlpl <- computemeandist(alldistances,"lpl_weighted_energy",netw)
      mdistspect <- computemeandist(alldistances,"spect_rad_weighted",netw)
      mdistlplspect <- computemeandist(alldistances,"lpl_spect_rad_weighted",netw)
      meandist <- rbind(meandist,mdistadj,mdistlpl,mdistspect,mdistlplspect)
    }
    meandist <- meandist[!duplicated(meandist),]
  }
  print("Saving normalized mean distances file")
  meandist <- meandist[!(grepl("weighted",meandist$disttype) & (meandist$networktype=="BINARY")),]
  meandist[abs(meandist$meannormdist)<0.000001,]$meannormdist=0
  fwrite(meandist,paste0(normresults,"meannormdistances.csv"),row.names=FALSE)
  print("Computing correlations...")
  networkmags <- read.csv(paste0(resultsdir,"networkmagnitudes.csv"))
  normdistavgs <- read.csv(paste0(normresults,"meannormdistances.csv"))
  normdistavgs$interaction <- ""
  if (sum(grepl("_PL_",normdistavgs$network))>1)
    normdistavgs[grepl("_PL_",normdistavgs$network),]$interaction = "PL"
  if (sum(grepl("_HP_",normdistavgs$network))>1)
    normdistavgs[grepl("_HP_",normdistavgs$network),]$interaction = "HP"
  if (sum(grepl("_SD_",normdistavgs$network))>1)
    normdistavgs[grepl("_SD_",normdistavgs$network),]$interaction = "SD"
  normdistavgs$links <- 0
  normdistavgs$nodes <- 0
  normdistavgs$weight <- 0
  normdistavgs$corlinks <- 0
  normdistavgs$cornodes <- 0
  normdistavgs$corweight <- 0
  lntype <- unique(normdistavgs$networktype)
  ldisttype <- unique(normdistavgs$disttype)
  lmodel <- unique(normdistavgs$MODEL)
  
  
  datacorr <- data.frame("disttype"=c(),"model"=c(),"corrlinks"=c(),"corrnodes"=c())
  nets <- unique(normdistavgs$network)
  for (i in 1:nrow(normdistavgs)){
     dtype <- normdistavgs$disttype[i]
     dfn <- networkmags[networkmags$Network==normdistavgs$network[i],]
     normdistavgs$links[i] <- dfn$Links[1]
     normdistavgs$nodes[i] <- dfn$NodesA[1]+dfn$NodesB[1]
     normdistavgs$weight[i] <- dfn$Weight[1]
     
  }
  
  for (disttype in ldisttype)
    for (model in lmodel){
      for (nclass in c("ALL")){# c(unique(normdistavgs$interaction),"ALL")){
      print(paste(disttype,model,nclass))
      corrdata <- meandist[meandist$MODEL==model & meandist$disttype==disttype,]
      if (nclass!="ALL")
        corrdata <- corrdata[grepl(nclass,corrdata$network),]
      if (nrow(corrdata)>0 )
        if(sd(corrdata$meannormdist)>0){
        corrdata$links <- 0
        corrdata$nodes <- 0
        corrdata$weight <- 0
        corrdata$connectance <- 0
        for (i in 1:nrow(corrdata)){
         dfn <- networkmags[networkmags$Network==corrdata$network[i],]
         corrdata[i,]$links <- dfn$Links[1]
         corrdata[i,]$nodes <- dfn$NodesA[1]+dfn$NodesB[1]
         corrdata[i,]$weight <- dfn$Weight[1]
         corrdata[i,]$connectance <- dfn$Connectance[1]
        }
        mycorrlinks <- cor(corrdata$meannormdist,corrdata$links,method = "pearson")
        mycorrnodes <- cor(corrdata$meannormdist,corrdata$nodes,method = "pearson")
        mycorrweight <- cor(corrdata$meannormdist,corrdata$weight,method = "pearson")
        mycorrconnectance <- cor(corrdata$meannormdist,corrdata$connectance,method = "pearson")
        
        if (!is.na(mycorrlinks))
          datacorr <- rbind(datacorr,data.frame("disttype"=disttype,"model"=model,
                                                "nclass"=nclass,
                                                "corrlinks"=mycorrlinks,
                                                "corrnodes"=mycorrnodes,
                                                "corrweight"=mycorrweight,
                                                "corrconnectance"=mycorrconnectance))
      }
    }
  }
  datacorr <- datacorr[!duplicated(datacorr), ]
  print("Saving correlations file")
  fwrite(datacorr,paste0(normresults,"corrnormdistances.csv"),row.names=FALSE)
}  
