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

LaplweightedEnergy <- function(eigenvalues,num_links,num_nodes){
  return(sum(eigenvalues^2))
}

LaplEnergy <- function(eigenvalues,num_links,num_nodes){
  return(sum(abs(eigenvalues-2*num_links/num_nodes)))
}

AdjEnergy <- function(eigenvalues){
  return(sum(abs(eigenvalues)))
}

AdjweightedEnergy <- function(eigenvalues){
  return(sum(eigenvalues^2))
  #return(sum(abs(eigenvalues)))
}

store_spect_values <- function(val,typem){
  return(data.frame("val"=val,"ind"=seq(1:length(val)),"type" = typem))
}

create_networkspect_dataframe <- function(weighted_network,bm,wm){
  if(weighted_network)
    networkspect <- list("spect_rad"=bm$adj_spect$values[1],
                         "adj_energy"=bm$adj_energy,
                         "lpl_energy"=bm$lpl_energy,
                         "algebraic_connectivity"=bm$lpl_spect$values[length(bm$lpl_spect$values)-1],
                         "spect_rad_weighted"=wm$adj_weighted_spect$values[1],
                         "adj_weighted_energy"=wm$adj_weighted_energy,
                         "lpl_weighted_energy"=wm$lpl_weighted_energy)
  else
    networkspect <- list("spect_rad"=bm$adj_spect$values[1],
                         "adj_energy"=bm$adj_energy,
                         "lpl_energy"=bm$lpl_energy,
                         "algebraic_connectivity"=bm$lpl_spect$values[length(bm$lpl_spect$values)-1],
                         "spect_rad_weighted"=-1,
                         "adj_weighted_energy"=-1,
                         "lpl_weighted_energy"=-1)
    return(networkspect)
}

# Updates Networ magnitudes file
store_network_magnitudes_global_file <- function(netw,nnm){
  # Write Network magnitudes to global results file
  network_data <- data.frame("Network"=netw,"Weighted"=nnm$weighted_network,
                             "NodesA"=nnm$nodes_a,"NodesB"=nnm$nodes_b,"Links"=nnm$num_links,
                             "Weight"=nnm$network_total_weight,"Model"="NETWORK")
  network_values <- cbind(network_data,as.data.frame(nnm$networkspect))
  
  network_values <- rbind(network_values,network_values)
  network_values[2,]$Model = "NESTED"
  for (i in which(names(network_values)=="spect_rad"):length(names(network_values)))
    network_values[2,i] = INCOMPLETE_MODEL_VALUE
  network_values <- rbind(network_values,network_values[2,])
  nested_df <- nnm$datamod[nnm$datamod$MODEL=="NESTED",]
  for (i in nested_df$ind)
    network_values[2,][i] <- nested_df[nested_df$ind==i,]$values
  network_values[3,]$Model = "WNESTED"
  nested_df <- nnm$datamod[nnm$datamod$MODEL=="WNESTED",]
  for (i in nested_df$ind)
    network_values[3,][i] <- nested_df[nested_df$ind==i,]$values
  for (modelname in nnm$mnames){
    network_null <- cbind(network_data,as.data.frame(lapply(nnm$modresults[[modelname]],mean)) )
    network_null$Model <- modelname
    network_values <- rbind(network_values,network_null)
  }
  NFile <- paste0(rdir,NetworkMagsFile)
  if (file.exists(NFile)){
    NMags <- read.csv(NFile)
    NMags <- NMags[NMags$Network!=netw,]
    NMags <- rbind(NMags,network_values)
  } 
  else       
    NMags <- network_values
  if (!plotzigs)
    write.csv(NMags,NFile,row.names = FALSE)
}

compute_bin_magnitudes <- function(result_analysis,matrix,network_nested_values,nodes_a,nodes_b,num_links,num_nodes,dfnested){
  dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="NETWORK",
                                        "NODF"=network_nested_values["NODF"],
                                        "wine"=network_nested_values["wine"],
                                        "links"=result_analysis$links,"totalweight"=sum(matrix)))

  mtx <- sq_adjacency(matrix, nodes_a, nodes_b)
  adj_sq_matrix <- mtx[[1]]  # binarized matrix
  adj_sq_weighted_matrix <- mtx[[2]] # original matrix, both are equal if the network is binary
  adj_spect <- eigen(adj_sq_matrix)
  adj_energy <- AdjEnergy(adj_spect$values)
  print(sprintf("Sum adjacency spectrum^2 %.2f Energy %.2f",sum((adj_spect$values)^2),adj_energy))
  lapl_matrix <- 0-adj_sq_matrix
  degrees <- rowSums(adj_sq_matrix)
  for (i in 1:nrow(lapl_matrix))
    lapl_matrix[i,i] <- degrees[i]
  lpl_spect = eigen(lapl_matrix)
  lpl_energy <- LaplEnergy(lpl_spect$values,num_links,num_nodes)
  print(sprintf("Sum Laplacian spectrum %.2f",sum(lpl_spect$values)))
  calc_values <- list("adj_spect"=adj_spect,"adj_energy"=adj_energy,"lpl_spect"=lpl_spect,
                      "lpl_energy"=lpl_energy,"adj_sq_matrix"=adj_sq_matrix,
                      "adj_sq_weighted_matrix"=adj_sq_weighted_matrix,
                      "lapl_matrix"=lapl_matrix,"dfnested"=dfnested)
  return(calc_values)
}

compute_weighted_magnitudes <- function(adj_sq_weighted_matrix)
{
  adj_weighted_spect <- eigen(adj_sq_weighted_matrix)
  spect_rad_weighted <- adj_weighted_spect$values[1]
  adj_weighted_energy <- AdjweightedEnergy(adj_weighted_spect$values)
  lapl_weighted_matrix <- 0-adj_sq_weighted_matrix
  sumweights <- rowSums(adj_sq_weighted_matrix)
  for (i in 1:nrow(lapl_weighted_matrix))
    lapl_weighted_matrix[i,i] <- sumweights[i]
  lpl_weighted_spect = eigen(lapl_weighted_matrix)
  lpl_weighted_energy <- LaplweightedEnergy(lpl_weighted_spect$values,num_links,num_nodes)
  print(sprintf("Sum weighted Laplacian spectrum %.2f",sum(lpl_weighted_spect$values)))
  calc_values <- list("adj_weighted_spect"=adj_weighted_spect,
                      "adj_weighted_energy"=adj_weighted_energy,
                      "spect_rad_weighted" = spect_rad_weighted,
                      "lpl_weighted_spect"=lpl_weighted_spect,
                      "lpl_weighted_energy"=lpl_weighted_energy)
  return(calc_values)
}

compute_statistic_indexes <- function(netw,weighted_network,datamod,networkspect,modresults,mnames)
{
  lheader <- length(networkspect)
  numvals <- (nrow(datamod)-lheader)/num_experiments
  ntry <- c()
  for (i in 1:lheader)
    ntry <- c(ntry,1:num_experiments)
  dtry <- replicate(lheader,0)
  for (i in 1:length(names(modresults)))
    dtry <- c(dtry,ntry)

  testvalues = data.frame("MODEL"=c(),"magnitude"=c(), "networkvalue"=c(),"modelmean"=c(),
                          "quantile"=c(),"distance"=c())
  testvalues = data.frame("MODEL"=c(),"magnitude"=c(), "networkvalue"=c(),"modelmean"=c(),
                          "quantile"=c(),"distance"=c())
  if (weighted_network)
    tmagnitudes <- list("spect_rad","adj_energy","lpl_energy",
                        "spect_rad_weighted", "adj_weighted_energy","lpl_weighted_energy")
  else
    tmagnitudes <- list("spect_rad","adj_energy","lpl_energy")
  for (modeln in mnames)
    for (magnitude in tmagnitudes)
    {
      datos <- datamod[(datamod$MODEL==modeln) & (datamod$ind==magnitude),]$values
      mvalue <- mean(datos)
      datored <- datamod[(datamod$MODEL=="NETWORK") & (datamod$ind==magnitude),]$values
      percentile <- ecdf(datos)
      testvalues <- rbind(testvalues,data.frame("MODEL"=modeln,
                                                "magnitude"=magnitude, "networkvalue"=datored,
                                                "modelmean"=mvalue,
                                                "quantile"=100*percentile(datored),
                                                "distance"=datored-mvalue))
    }
  if (!plotzigs)
    write.csv(testvalues,paste0(rdir,"TESTVALUES_",netw),row.names = FALSE)
}

process_network_null_models <- function(netw,result_analysis,num_experiments,mnamesbin,mnamesweighted,weightrf)
{
  orimatrix <- result_analysis$matrix
  dfnested <- create_empty_dfnested()
  weighted_network <- sum(orimatrix > 1)>0
  if ((weighted_network) && (weightrf != "none"))
  {
    if (weightrf=="sqrt")
      trfmatrix <- sqrt(orimatrix)
    if (weightrf=="ln")
      trfmatrix <- log(1+orimatrix)
    trfmatrix <- ceiling(trfmatrix)
  } else {
    trfmatrix <- orimatrix
  }
  network_total_weight <- sum(trfmatrix)
  network_nested_values <- nested(trfmatrix, c("NODF","wine"))
  nodes_a <- result_analysis$num_guild_a
  nodes_b <- result_analysis$num_guild_b
  num_links <- result_analysis$links
  num_nodes <- nodes_a+nodes_b
  bmag <- compute_bin_magnitudes(result_analysis,trfmatrix,network_nested_values,nodes_a,nodes_b,num_links,num_nodes,dfnested)
  #Nested model
  pnm <- process_nested_model_bin(nodes_a, nodes_b, num_links)
  dfnested <- rbind(dfnested,bmag$dfnested,data.frame("network"=netw,"model"="NESTED",
                                        "NODF"=pnm$nnst["NODF"],
                                        "wine"=pnm$nnst["wine"],
                                        "links"= sum(pnm$nstmodel>0),
                                        "totalweight"= sum(pnm$nstmodel)))
  
  # Weighted magnitudes, only for weighted networks
  if (weighted_network){
    wmag <- compute_weighted_magnitudes(bmag$adj_sq_weighted_matrix)
    pnw <- process_nested_model_weighted(nodes_a, nodes_b, num_links, trfmatrix)
    dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="WNESTED",
                                          "NODF"=pnw$weightednnst["NODF"],
                                          "wine"=pnw$weightednnst["wine"],
                                          "links"= sum(pnw$weightednstmodel>0),"totalweight"= sum(pnw$weightednstmodel)))
  
  } else {
    lpl_weighted_energy = INCOMPLETE_MODEL_VALUE
    wmag <- NULL
    pnw <- NULL
  }

  save_null_model(netw,dirnulls,pnm$nstmodel,"WNESTED",pnm$nstadj_spect,pnm$nstlpl_spect,
                    wadj_spect_val=pnw$weightednstadj_spect,wlpl_spect_val=pnw$weightednstlpl_spect)
  networkspect <- create_networkspect_dataframe(weighted_network,bmag,wmag)
  mnames <- mnamesbin
  if (weighted_network)
    mnames <- c(mnames,mnamesweighted)
  nmodels <- length(mnames)
  model_full <- create_models_list(mnames,NULL)
  nname <- gsub(".csv","",netw)
  nullsinfo <- create_nullsinfo(num_experiments)
  modresults <- lapply(1:nmodels, function(x) nullsinfo)
  names(modresults) <- mnames
  specresults <- lapply(1:nmodels, function(x) c())
  names(specresults) <- mnames
  datamod <- create_datamodels(weighted_network,networkspect,pnm,pnw)
  nodes_guild_a <- seq(1,result_analysis$num_guild_a)  # guild a top rows/left columns
  nodes_guild_b <- seq(1,result_analysis$num_guild_b)  # guild b bottom rows/right columns
  for (k in 1:num_experiments) {
    if (!k%%20)
      print(paste("experiment",k))
    for (tmodel in mnames)
    {
      for (i in 1:20){     # This section is protectec against some null model failures, specially SWAP
        result <- tryCatch(
          {
                p <- gen_null_model(tmodel,datamod,result_analysis,trfmatrix,nodes_guild_a,nodes_guild_b,
                            num_links,num_nodes,nmagnitudes,weighted_network)
          },
          error = function(cond) {
            message(paste("Null model generation error:", tmodel))
            message("Here's the original error message:")
            message(conditionMessage(cond))
            # Choose a return value in case of error
            NA
          }
        )
        if (class(result) != "try-error") {   # Null model succesfully generated
          break
        }
      }
      model_full[[tmodel]] <- p$incidmatrix
      # Save null matrices of first experiment
      if (k==1)
        save_null_model(netw,dirnulls,p$incidmatrix,tmodel,p$resp$null_spect,
                        p$resp$lpl_spect_nulls,wadj_spect_val=p$resp$null_weighted_spect,
                        wlpl_spect_val=p$resp$lpl_weighted_spect_nulls)  
      modresults[[tmodel]][k,] <- p$mres
      specresults[[tmodel]] <- p$resp
      dfnested <- rbind(dfnested,data.frame("network"=netw,"model"=tmodel,"NODF"=p$resp$NODF,"wine"=p$resp$wine,
                                            "links"= sum(p$incidmatrix>0),"totalweight"= sum(p$incidmatrix)))
    }
  }
  calc_values <- list("matrix"=trfmatrix,"modresults"=modresults,"specresults"=specresults,"datamod"=datamod,
                      "dfnested" = dfnested,"networkspect"=networkspect,
                      "mnames"=mnames,"weighted_network"=weighted_network,
                      "nodes_a"=nodes_a,"nodes_b"=nodes_b,"num_links"=num_links,
                      "network_total_weight"=network_total_weight,
                      "bmag"=bmag,"wmag"=wmag,"pnm"=pnm,"pnw"=pnw)
  return(calc_values)
}

create_empty_dfnested <- function()
{
  p <- data.frame("network"=c(),"model"=c(),"NODF"=c(),"wine"=c(), 
                  "links"=c(),"totalweight"=c())
  return(p)
}

store_nested_values <- function(dfnested){
  dfnested <- rbind(dfnested,data.frame("network"=netw,"model"=tmodel,"NODF"=p$resp$NODF,"wine"=p$resp$wine,
                                        "links"= sum(p$incidmatrix>0),"totalweight"= sum(p$incidmatrix)))
  return(dfnested)
}

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
  # Compute spectral magnitudes of null NESTED and WNESTED models
  store_network_magnitudes_global_file(netw,nnm)
  if (plottofile)
    plot_all_distr(netw,plotzigs,nnm)
  if (plotzigs)
    plot_ziggurats(nnm,dirnulls,odir)
}
# Configuration parameters
seed <- 122
num_experiments <- 50
plottofile <- TRUE # Save individual network distributions plot
plotzigs <- FALSE  # Plotting ziggurats of all models is rather slow. So when TRUE magnitudes are
                   # not saved. Run the script with a big number of experiments (~1000) to compute
                   # magnitudes and plotzigs FALSE. Run it again selecting just the networks you need
                   # to plot and a small number of experiments (~10)
NetworkMagsFile <- "NetworkMagnitudes.csv" # Stores network magnitudes and average spectral measures
nmagnitudes <- list("spect_rad","adj_energy","lpl_energy")#,"algebraic_connectivity")
# Null models for binary/binarized and weighted networks
mnamesbin <- c("RND","MGEN","SHUFFLE","VAZ","SYTR")
mnamesweighted <- c("SWAP","WRND","BVAZ","BSHUFFLE","PATEFIELD")
MIN_LINKS_SIZE <- 20  # Smaller networks are discarded

# Here, the list of data files to process
filenames <- Sys.glob(paste0(datadir,"*PL*055*.csv"))
# Network names
lnetw <- gsub(datadir,"",filenames)
for (netw in (lnetw)){                 # Each network
  for (weightrf in lweightrf)
    network_null_spectral_distances(netw,weightrf,numexperiments,mnamesbin,mnamesweighted,datadir,rdir,odir,dirnulls,
                                              guild_a_label = "Plant", guild_b_label = "Pollinator",min_links=20,
                                              plottofile=TRUE,plotzigs=FALSE,networkmagsfile="NetworkMagnitudes.csv")
}