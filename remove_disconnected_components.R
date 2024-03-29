rm(list = ls())
library(kcorebip)
source("config_values.R")
source("manage_null_models.R")

dir.create(dataprocessed, showWarnings = FALSE)
filenames <- Sys.glob(paste0(datadir,"*_*.csv"))
# Network names
lnetw <- gsub(datadir,"",filenames)
datanets <- data.frame("Network"=c(),"FullyConnected"=c())
for (netw in lnetw){
  print(netw)
  result_analysis <- analyze_network(directory = datadir, netw, only_NODF = TRUE)
  madj <- sq_adjacency(result_analysis$matrix,result_analysis$num_guild_a, result_analysis$num_guild_b)[[1]]
  mlapl <- create_laplacian_matrix(madj)
  ev <- sort(eigen(mlapl,only.values = TRUE)$values)
  if (ev[2]>0.000001){
    print("connected")
    fullyconn <- 1
  } else {
    print("disconnected")
    fullyconn <- 0
    rg <- V(result_analysis$graph)
    g <-  rg[rg$kradius != Inf]
    outsider <- rg[rg$kradius == Inf]
    pls <- as.numeric(gsub("pl","",outsider$name[grepl("pl",outsider$name)]))
    pols <- as.numeric(gsub("pol","",outsider$name[grepl("pol",outsider$name)]))
    trfmatrix <- result_analysis$matrix[-c(pols),-c(pls)]
    write.csv(trfmatrix,paste0(dataprocessed,gsub(".csv","",netw),"_GC.csv"))
  }
  datanets <- rbind(datanets,data.frame("Network"=gsub(".csv","",netw),"FullyConnected"=fullyconn))
}
write.csv(datanets,paste0(dataprocessed,"dataconnection.csv"),row.names=FALSE)