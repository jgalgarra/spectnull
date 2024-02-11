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

compute_nested_magnitudes <- function(network_values){
  network_values <- rbind(network_values,network_values)
  network_values[2,]$Model = "NESTED"
  for (i in which(names(network_values)=="spect_rad"):length(names(network_values)))
    network_values[2,i] = INCOMPLETE_MODEL_VALUE
  network_values <- rbind(network_values,network_values[2,])
  nested_df <- datamod[datamod$MODEL=="NESTED",]
  for (i in nested_df$ind)
    network_values[2,][i] <- nested_df[nested_df$ind==i,]$values
  network_values[3,]$Model = "WNESTED"
  nested_df <- datamod[datamod$MODEL=="WNESTED",]
  for (i in nested_df$ind)
    network_values[3,][i] <- nested_df[nested_df$ind==i,]$values
  for (modelname in mnames){
    network_null <- cbind(network_data,as.data.frame(lapply(modresults[[modelname]],mean)) )
    network_null$Model <- modelname
    network_values <- rbind(network_values,network_null)
  }
  NFile <- paste0(rdir,NetworkMagsFile)
  if (file.exists(NFile)){
    NMags <- read.csv(NFile)
    NMags <- NMags[NMags$Network!=nname,]
    NMags <- rbind(NMags,network_values)
  } 
  else       
    NMags <- network_values
  if (!plotzigs)
    write.csv(NMags,NFile,row.names = FALSE)
}

compute_bin_magnitudes <- function(){
  dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="NETWORK",
                                        "NODF"=network_nested_values["NODF"],
                                        "wine"=network_nested_values["wine"],
                                        "links"=result_analysis$links,"totalweight"=sum(result_analysis$matrix)))

  mtx <- sq_adjacency(result_analysis$matrix, nodes_a, nodes_b)
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
                      "adj_sq_weighted_matrix"=adj_sq_weighted_matrix,"lapl_matrix"=lapl_matrix)
  return(calc_values)
}

compute_weighted_magnitudes <- function(adj_sq_weighted_matrix,adj_weighted_spect)
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

compute_statistic_indexes <- function(networkspect,netw)
{
  lheader <- length(networkspect)
  numvals <- (nrow(datamod)-lheader)/num_experiments
  ntry <- c()
  for (i in 1:lheader)
    ntry <- c(ntry,1:num_experiments)
  dtry <- replicate(lheader,0)
  for (i in 1:length(names(modresults)))
    dtry <- c(dtry,ntry)
  if (!plotzigs){  
    write.csv(datamod,paste0(rdir,"MODS_",netw),row.names = FALSE)  # Model results
    write.csv(dfnested,paste0(rdir,"NESTED_",netw),row.names = FALSE)  # Model results
  }
  testvalues = data.frame("MODEL"=c(),"magnitude"=c(), "networkvalue"=c(),"modelmean"=c(),
                          "quantile"=c(),"distance"=c(),"reldist"=c())
  testvalues = data.frame("MODEL"=c(),"magnitude"=c(), "networkvalue"=c(),"modelmean"=c(),
                          "quantile"=c(),"distance"=c(),"reldist"=c())
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
                                                "distance"=datored-mvalue,"reldist"=(datored-mvalue)/datored))
    }
  if (!plotzigs)
    write.csv(testvalues,paste0(rdir,"TESTVALUES_",netw),row.names = FALSE)
}

filenames <- Sys.glob(paste0(datadir,"*PL*004*.csv"))
seed <- 122
num_experiments <- 1
plottofile <- TRUE
plotzigs <- FALSE  # Plotting ziggurats of al models is rather slow. So when TRUE magnitudes are
                   # not saved. Run the script with a big number of experiments (~1000) to compute
                   # magnitudes and plotzigs FALSE. Run it again selecting just the networks you need
                   # to plot and a small number of experiments (~10)
NetworkMagsFile <- "NetworkMagnitudes.csv"

fract_dummy <- c(5)
nmagnitudes <- list("spect_rad","adj_energy","lpl_energy")#,"algebraic_connectivity")

lnetw <- gsub(datadir,"",filenames)

for (weightrf in lweightrf){
  create_dirs(weightrf)
  for (netw in (lnetw)){
    print(paste("Network:",netw,"transform",weightrf))
    dfnested <- data.frame("network"=c(),"model"=c(),"NODF"=c(),"wine"=c(), "links"=c(),"totalweight"=c())
    result_analysis <- analyze_network(directory = datadir, netw, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
    if (result_analysis$links < 21)
      next
    weighted_network <- sum(result_analysis$matrix > 1)>0
    if ((weighted_network) && (weightrf != "none"))
    {
      if (weightrf=="sqrt")
        result_analysis$matrix <- sqrt(result_analysis$matrix)
      if (weightrf=="ln")
        result_analysis$matrix <- log(1+result_analysis$matrix)
      result_analysis$matrix <- ceiling(result_analysis$matrix)
    }
    
    network_total_weight <- sum(result_analysis$matrix)  
    network_nested_values <- nested(result_analysis$matrix, c("NODF","wine"))
    nodes_a <- result_analysis$num_guild_a
    nodes_b <- result_analysis$num_guild_b
    num_links <- result_analysis$links
    num_nodes <- nodes_a+nodes_b
    bmag <- compute_bin_magnitudes()
    #Nested model
    pnm <- process_nested_model_bin(nodes_a, nodes_b, num_links)
    dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="NESTED",
                                          "NODF"=pnm$nnst["NODF"],
                                          "wine"=pnm$nnst["wine"],
                                          "links"= sum(pnm$nstmodel>0),"totalweight"= sum(pnm$nstmodel)))
    
    # Weighted magnitudes, only for weighted networks
    if (weighted_network){
      wmag <- compute_weighted_magnitudes(bmag$adj_sq_weighted_matrix,bmag$adj_weighted_spect)
      pnw <- process_nested_model_weighted(nodes_a, nodes_b, num_links)
      dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="WNESTED",
                                            "NODF"=pnw$weightednnst["NODF"],
                                            "wine"=pnw$weightednnst["wine"],
                                            "links"= sum(pnw$weightednstmodel>0),"totalweight"= sum(pnw$weightednstmodel)))
      
    } else {
      lpl_weighted_energy = INCOMPLETE_MODEL_VALUE
    }
    networkspect <- create_networkspect_dataframe(weighted_network,bmag,wmag)
    mnames <- c("RND","MGEN","SHUFFLE","VAZ","SYTR")
    if (weighted_network)
      mnames <- c(mnames,c("SWAP","WRND","BVAZ","BSHUFFLE","PATEFIELD"))#"WSYTR"
    nmodels <- length(mnames)
    model_full <- create_models_list(mnames,NULL)
    nname <- gsub(".csv","",netw)
    nullsinfo <- create_nullsinfo(num_experiments)
    modresults <- lapply(1:nmodels, function(x) nullsinfo)
    names(modresults) <- mnames
    specresults <- lapply(1:nmodels, function(x) c())
    names(specresults) <- mnames
    datamod <- create_datamodels(networkspect,pnm,pnw)
    nodes_guild_a <- seq(1,result_analysis$num_guild_a)  # guild a top rows/left columns
    nodes_guild_b <- seq(1,result_analysis$num_guild_b)  # guild b bottom rows/right columns
    for (k in 1:num_experiments) {
      if (!k%%100)
        print(paste("experiment",k))
      for (tmodel in mnames)
      {
        p <- gen_null_model(tmodel,result_analysis,nodes_guild_a,nodes_guild_b,
                            num_links,num_nodes,nmagnitudes,weighted_network)
        model_full[[tmodel]] <- p$incidmatrix
        # Save null matrices of first experiment
        if (k==1)
          save_null_model(p$incidmatrix,netw,dirnulls,tmodel)
        modresults[[tmodel]][k,] <- p$mres
        specresults[[tmodel]] <- p$resp
        dfnested <- rbind(dfnested,data.frame("network"=netw,"model"=tmodel,"NODF"=p$resp$NODF,"wine"=p$resp$wine,
                                              "links"= sum(p$incidmatrix>0),"totalweight"= sum(p$incidmatrix)))
      }
    }
  
    # Save null model results
    for (m in names(modresults)){
      sm <- stack(modresults[[m]])
      sm$MODEL <- m
      datamod <- rbind(datamod,sm)
      if (!weighted_network)
        datamod <- datamod[!grepl("weighted",datamod$ind),]
    }
    
    compute_statistic_indexes(networkspect,netw)
    # Write Network magnitudes to global results file
    network_data <- data.frame("Network"=nname,"Weighted"=weighted_network,"NodesA"=nodes_a,"NodesB"=nodes_b,"Links"=num_links,
                   "Weight"=network_total_weight,"Model"="NETWORK")
    network_values <- cbind(network_data,as.data.frame(networkspect))
    # NESTED values
    compute_nested_magnitudes(network_values)
    if (plottofile)
      plot_all_distr(mnames,bmag,wmag)
    
  }
}
