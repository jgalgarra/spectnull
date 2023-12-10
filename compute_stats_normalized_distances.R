# This script computes the average normalized distances.
# It must be run after normalize_distances.R

normresults <- "results/normalized/"
results <- "results/"
print("Reading all normalized distances from file")
alldistances <- read.csv(paste0(normresults,"ALLNORMALIZED.csv"))
alldistances$network <- gsub("MINMAX_MODS_","",alldistances$network)
networktypes <- read.csv(paste0(normresults,"networktypes.csv"))

print("Computing averages")
lnetw <- unique(alldistances$network)
alldistances$type <- "BINARY"
for (netw in lnetw){
  ntype <- networktypes[networktypes$network==netw,]
  alldistances[alldistances$network==netw,]$type <- ntype$type
#   adj_norm_net <- alldistances[(alldistances$network==netw) & (alldistances$MODEL =="NETWORK") & (alldistances$ind=="adj_norm_energy"),]$values
#   lpl_norm_net <- alldistances[(alldistances$network==netw) & (alldistances$MODEL =="NETWORK") & (alldistances$ind=="lpl_norm_energy"),]$values
#   spectrad_norm_net <- alldistances[(alldistances$network==netw) & (alldistances$MODEL =="NETWORK") & (alldistances$ind=="spect_rad_norm"),]$values
#   alldistances[(alldistances$network==netw) & (alldistances$ind=="adj_norm_energy"),]$normdist = alldistances[(alldistances$network==netw) & (alldistances$ind=="adj_norm_energy"),]$values - adj_norm_net
#   alldistances[(alldistances$network==netw) & (alldistances$ind=="lpl_norm_energy"),]$normdist = alldistances[(alldistances$network==netw) & (alldistances$ind=="lpl_norm_energy"),]$values - lpl_norm_net
#   alldistances[(alldistances$network==netw) & (alldistances$ind=="spect_rad_norm"),]$normdist = alldistances[(alldistances$network==netw) & (alldistances$ind=="spect_rad_norm"),]$values - spectrad_norm_net
#   
#   print(paste(netw,adj_norm_net))
}

meandist <- data.frame("network"=c(),"meannormdist"=c(),"MODEL"=c(), "disttype"=c(), "networktype"=c())
for (netw in lnetw){
  adj_data <- alldistances[(alldistances$ind=="adj_norm_energy") & (alldistances$network==netw),]
  models <- unique(adj_data$MODEL)
  for (m in models){
    if (m!="NETWORK")
      meandist <- rbind(meandist,data.frame("network"=netw,
                                            "meannormdist"=mean(adj_data[adj_data$MODEL==m,]$normdist),
                                            "MODEL"=m,
                                            "disttype"="adj_norm_energy",
                                            "networktype"=adj_data$type[1]))
  }
  lpl_data <- alldistances[(alldistances$ind=="lpl_norm_energy") & (alldistances$network==netw),]
  models <- unique(lpl_data$MODEL)
  for (m in models){
    if (m!="NETWORK")
      meandist <- rbind(meandist,data.frame("network"=netw,
                                            "meannormdist"=mean(lpl_data[lpl_data$MODEL==m,]$normdist),
                                            "MODEL"=m,
                                            "disttype"="lpl_norm_energy",
                                            "networktype"=lpl_data$type[1]))
    spect_data <- alldistances[(alldistances$ind=="spect_rad_norm") & (alldistances$network==netw),]
    models <- unique(spect_data$MODEL)
    for (m in models){
      if (m!="NETWORK")
        meandist <- rbind(meandist,data.frame("network"=netw,
                                              "meannormdist"=mean(spect_data[spect_data$MODEL==m,]$normdist),
                                              "MODEL"=m,
                                              "disttype"="spect_rad_norm",
                                              "networktype"=spect_data$type[1]))
    }
  }
}
print("Saving normalized mean distances file")
write.csv(meandist,paste0(normresults,"meannormdistances.csv"),row.names=FALSE)

print("Computing correlations...")
networkmags <- read.csv(paste0(results,"networkmagnitudes.csv"))
normdistavgs <- read.csv(paste0(normresults,"meannormdistances.csv"))
normdistavgs$interaction <- ""
if (sum(grepl("_PL_",normdistavgs$network))>1)
  normdistavgs[grepl("_PL_",normdistavgs$network),]$interaction = "HP"
if (sum(grepl("_HP_",normdistavgs$network))>1)
  normdistavgs[grepl("_HP_",normdistavgs$network),]$interaction = "HP"
if (sum(grepl("_SD_",normdistavgs$network))>1)
  normdistavgs[grepl("_SD_",normdistavgs$network),]$interaction = "SD"
normdistavgs$links <- 0
normdistavgs$nodes <- 0
normdistavgs$corlinks <- 0
normdistavgs$cornodes <- 0
lntype <- unique(normdistavgs$networktype)
ldisttype <- unique(normdistavgs$disttype)
lmodel <- unique(normdistavgs$MODEL)

datacorr <- data.frame("nettype"=c(),"disttype"=c(),"model"=c(),"corrlinks"=c(),"corrnodes"=c())
nets <- unique(normdistavgs$network)
for (i in 1:nrow(normdistavgs)){
  dtype <- normdistavgs$disttype[i]
  dfn <- networkmags[networkmags$Network==normdistavgs$network[i],]
  normdistavgs$links[i] <- dfn$Links[1]
  normdistavgs$nodes[i] <- dfn$NodesA[1]+dfn$NodesB[1]
}
for (i in 1:nrow(normdistavgs[normdistavgs$network==nets[1],])){
  nmdata <- normdistavgs[i,]
  for (ntype in lntype){
    for (nmodel in lmodel){
      distcorrdata <- normdistavgs[normdistavgs$MODEL==nmodel &  
                                     normdistavgs$disttype==nmdata$disttype & 
                                     normdistavgs$networktype == ntype,]
      mycorrlinks <- cor(distcorrdata$meannormdist,distcorrdata$links,method = "spearman")
      mycorrnodes <- cor(distcorrdata$meannormdist,distcorrdata$nodes,method = "spearman")
      if (!is.na(mycorrlinks)){
        datacorr <- rbind(datacorr,data.frame("nettype"=ntype, "disttype"=nmdata$disttype,
                                              "model"=nmodel,"corrlinks"=abs(mycorrlinks),"corrnodes"=abs(mycorrnodes)))
      }
    }
  }
}
datacorr <- datacorr[!duplicated(datacorr), ]
print("Saving correlations file")
write.csv(datacorr,paste0(normresults,"corrnormdistances.csv"),row.names=FALSE)

