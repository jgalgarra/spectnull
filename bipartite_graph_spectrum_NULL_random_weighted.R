rm(list = ls())
library(kcorebip)
library(ggplot2)
library(patchwork)
library(forcats)
source("synthrade_nullmodel.R")

INCOMPLETE_MODEL_VALUE <<- -1

# In the original files guild a is stored by columns and guild_b by rows

get_colors <- function(models,colorlist){
  lnames <- names(colorlist)
  lcolors <- colorlist
  for (i in lnames){
    if (!(i %in% models)){
      lnames <- lnames[-which(lnames==i)]
      lcolors <- lcolors[-which(names(lcolors)==i)]
    }
  }
  calc_values <- list("lnames"= lnames, "lcolors" = lcolors)
  return(calc_values)
}

remove_model_data <- function(m,mod){
  if (sum(names(m)==mod)==0)
    return(m)
  else
    return(m[-which(names(m)==mod)])
}

plot_distr_null <- function(df,nvalue,networkname="",title="",nestedvalue=""){
  lmodels <- unique(df$type)
  clist <- get_colors(lmodels,colormodels)
  plot <- ggplot(data=df)+geom_histogram(aes(x=vals,fill=type),bins=30,position="identity",alpha=0.6)+
      xlab("")+ scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
    theme_bw()+
    theme(legend.position = "bottom",legend.text = element_text(size = 8),
          panel.grid = element_blank(),
          legend.key.size = unit(0.3, 'cm'))+
          guides(fill=guide_legend(nrow=2, byrow=TRUE))
  plot <- plot+geom_vline(data = data.frame("val"=nvalue), aes(xintercept = val), 
                          color = "blue", size=0.5,alpha=0.5)+ggtitle(title)
  if (nestedvalue!="")
    plot <- plot+geom_vline(data = data.frame("val"=nestedvalue), aes(xintercept = val), 
                            color = "red", size=0.75,alpha=0.5, linetype="dotted")+ggtitle(title)
  
  return(plot)
}

plot_spectr <- function(df,title=""){
  df$type <- fct_relevel(df$type, "RND", "VAZ", "NETWORK")
  if (nrow(df)<20)
    psize <- 3
  else
    psize <- 1.5
  lmodels <- unique(df$type)
  clist <- get_colors(lmodels,colormodels)
  plot <- ggplot(data=df)+geom_point(aes(x=ind,y=val,color=type),size=psize,
                                     position="identity",alpha=0.5)+
                          xlab("")+ylab("Eigen values")+theme_bw()+ggtitle(title)+
                          scale_y_continuous(expand=c(0,0))+
                          scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
                          theme_bw()+ theme(legend.position = "bottom",legend.text = element_text(size = 8),
                                            legend.key.size = unit(0.3, 'cm'),
                                            axis.title.x = element_text(size = 1))+
                          guides(color=guide_legend(title= "Null model"))
  return(plot)
}

plot_bipartite <- function(bg, aspect_ratio = 9/35, vframecolor = "grey70", vlabelcex = 4,
                           vsize = 4, vcolor = c("lightblue","pink2"), labelcolor = c("blue","red"),
                           framedisp = FALSE,  color_link = "grey50", vertical = FALSE)
{
  l <- layout.bipartite(bg)
  if (vertical)
    la <- l[, c(2,1)]
  else
    la <- l
  plot.igraph(bg, layout= la,asp=aspect_ratio,vertex.frame.color=vframecolor,
              vertex.label.cex=vlabelcex,vertex.label.color="black",
              vertex.size=vsize, edge.color= color_link, frame = framedisp,
              vertex.color=vcolor[V(bg)$type+1])
}

# Plots the ziggurat diagram of the model
plot_null_model <- function(modeldata,modeltype,netname,dirn,plotdir,ptfile){
  filenull <- paste0(gsub(".csv","",netname),"_NULL_",modeltype,".csv")
  print(paste("FILENULL",filenull))
  write.csv(modeldata,paste0(dirnulls,filenull))
  result_analysis <- analyze_network(directory = dirn, filenull, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  pgr_null <- ziggurat_graph(dirn,filenull,plotsdir=plotdir,print_to_file = ptfile,show_title = FALSE,weighted_links = "log10")
}

save_nested_model <- function(nstm,netname,dirn,nestedname){
  filenull <- paste0(gsub(".csv","",netname),"_NULL_",nestedname,".csv")
  print(paste("FILENULL",filenull))
  write.csv(nstm,paste0(dirnulls,filenull))
}

# Plots the probability distributions of all models for a particular magnitude and network
plot_distributions <- function(mresults,magnitude,stitle,vtitle,nname,fdummy,nestedv=""){
  mdls <- names(mresults)
  col <- which(names(mresults[[1]])==magnitude)
  df_nulls <- data.frame("vals"=mresults[[1]][,col],"type"=mdls[[1]])
  for (i in 2:length(mdls))
    df_nulls <- rbind(df_nulls,data.frame("vals"= mresults[[mdls[i]]][,col],"type"=mdls[i]))
  pimage <- plot_distr_null(df_nulls,vtitle,networkname=nname,title=sprintf("%s %.2f",stitle,vtitle),nestedvalue = nestedv)
  calc_values <- list("df_nulls"= df_nulls, "pimage" = pimage)
  return(calc_values)
}

LaplweightedEnergy <- function(eigenvalues,num_links,num_nodes){
  return(sum(eigenvalues^2))
}

LaplEnergy <- function(eigenvalues,num_links,num_nodes){
  return(sum(abs(eigenvalues-2*num_links/num_nodes)))
}

AdjEnergy <- function(eigenvalues){
  return(sum(abs(eigenvalues)))
}

erdosrmodel <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links,weighted=FALSE){
  model <- sample_bipartite(
    length(nodes_guild_a),
    length(nodes_guild_b),
    type = "gnm",
    m=num_links,
    directed = FALSE,
    mode = "all"
  )
  mmodel <- as.matrix(as_adjacency_matrix(model))
  mo <- mmodel[(1+length(nodes_guild_a)):nrow(mmodel),(1:length(nodes_guild_a))]
  if (weighted){
    empmatrix <- unname(adjmatrix[unname(adjmatrix)>1])
    ws <-  shuffle(unname(adjmatrix[unname(adjmatrix)>1]))
    mo[mo==1] <- ws
  }
  return(mo)
}

rndnullmodel_build <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links,allconnected = FALSE)
{
  mmatrix <- adjmatrix * 0
  # Throw links at random, some nodes may remain disconnected
  if (!allconnected)
    mmatrix[sample(1:length(adjmatrix),num_links,replace = FALSE)] <- 1
  else{
    for (i in 1:nrow(mmatrix))
      mmatrix[i,sample(1:ncol(mmatrix),1)] <- 1   # One link per row
    for (i in 1:ncol(mmatrix))
      mmatrix[sample(1:nrow(mmatrix),1),i] <- 1   # One link per column
    mmatrix[sample(1:length(adjmatrix),num_links-sum(mmatrix),replace = FALSE)] <- 1
  }
  return(mmatrix)

}

parsebinoutcome <- function(nrand){
  nums <- c()
  count <- 0
  for (i in 1:length(nrand)){
    if (nrand[i]==0){
      if (count > 0)
        nums <- c(nums,count)
      count = 0
      nums <- c(nums,0)
      
    }
    else
      count=count+1
  }
  return(unlist(nums))
}


# Weighted Random Model following Garlaschelli's algorithm

WRG_build <- function(empmatrix,nodes_guild_a,nodes_guild_b,num_links)
{
  binmatrix <- erdosrmodel(empmatrix,nodes_guild_a,nodes_guild_b,num_links)
  meanemp <- mean(empmatrix[empmatrix>0])
  #sdemp <- sd(empmatrix[empmatrix>0])
  wmatrix <- binmatrix
  rtokens <- rnorm(num_links,meanemp,meanemp/2)
  rtokens[rtokens<0] <- 0
  rtokens <- floor(rtokens)
  # W <- sum(empmatrix)-num_links
  # na <- length(nodes_guild_a)
  # nb <- length(nodes_guild_b)
  # N <- na*nb
  # p <-  W/(N+W)
  # mmatrix <- unname(empmatrix * 0)
  # #rtokens <- rbinom(1.2*W,W/N,p)
  # rtokens <- rbinom(num_links,ceiling(W/num_links),p)
  
  k=1
  for (i in nodes_guild_b)
    for (j in nodes_guild_a)
      if (wmatrix[i,j]>0){
        wmatrix[i,j] <- 1+rtokens[k]
        k <- k+1
      }
  return(wmatrix)
}


# WRG_build <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links)
# {
#   W <- sum(adjmatrix)
#   na <- length(nodes_guild_a)
#   nb <- length(nodes_guild_b)
#   N <- na*nb
#   p <-  W/(N+W)
#   mmatrix <- unname(adjmatrix * 0)
#   ntrials <- 1.2*W/p
#   rtokens <- rbinom(ntrials,1,p)
#   sweights <- parsebinoutcome(rtokens)
#   sweights <- sweights[1:(nrow(mmatrix)*ncol(mmatrix))]
#   mmatrix <- matrix(sweights,nrow=nrow(mmatrix))
#   return(mmatrix)
# }
# 



# Intercambia un porcentaje de enlaces

dummynullmodel_build <- function(nodes_guild_a,nodes_guild_b,num_links,re,swap_perc)
{
  mmatrix <- re$matrix
  swap_links <- max(1,round(0.01 * swap_perc * num_links))
  ones <- which(result_analysis$matrix > 0)
  zeroes <- which(result_analysis$matrix == 0)
  j = 1
  found <- FALSE
  while ((!found) & (j<10)){
    if (j>1)
      print(paste("dummy",j))
    pmatrix <- mmatrix
    onestozeroes <- sample(ones,swap_links,replace=FALSE)
    zeroestoones <- sample(zeroes,swap_links,replace=FALSE)
    pmatrix[onestozeroes] <- 0
    pmatrix[zeroestoones] <- 1
    j <- j +1
    found <- ((sum(rowSums(pmatrix))>0) & (sum(colSums(pmatrix))>0))
  } 
  return(pmatrix)
}

null_model_process <- function(admatrix,re,weighted_network)
{
  mtx <- sq_adjacency(admatrix,re$num_guild_a,re$num_guild_b)
  model_matrix <- mtx[[1]]  # binarized matrix
  model_weighted_matrix <- mtx[[2]] # original matrix
  
  nullNMP_spect <- eigen(model_matrix)
  lapl_nullNMP <- 0-model_matrix
  degreesnullNMP <- rowSums(model_matrix)
  for (i in 1:nrow(lapl_nullNMP))
    lapl_nullNMP[i,i] <- degreesnullNMP[i]
  lpl_spect_nullNMP = eigen(lapl_nullNMP)

  if (weighted_network){
    nullNMP_weighted_spect <- eigen(model_weighted_matrix)
    lapl_weighted_nullNMP <- 0-model_weighted_matrix
    nullweightedNMP <- rowSums(model_weighted_matrix)
    for (i in 1:nrow(lapl_weighted_nullNMP))
      lapl_weighted_nullNMP[i,i] <- nullweightedNMP[i]
    lpl_weighted_spect_nullNMP = eigen(lapl_weighted_nullNMP)
  } else {
    lpl_weighted_spect_nullNMP = lpl_spect_nullNMP
  }
  
  spect_rad_nulls_NMP <- nullNMP_spect$values[1]
  adj_energy_nulls_NMP <- AdjEnergy(nullNMP_spect$values)
  lpl_energy_nulls_NMP <- LaplEnergy(lpl_spect_nullNMP$values,num_links,num_nodes)
  algebraic_connectivity_nulls_NMP <- lpl_spect_nullNMP$values[length(lpl_spect_nullNMP$values)-1]
  if (weighted_network){
    spect_rad_weighted_nulls_NMP <- nullNMP_weighted_spect$values[1]
    lpl_weighted_energy_nulls_NMP <- LaplweightedEnergy(lpl_spect_nullNMP$values,num_links,num_nodes)
    adj_weighted_energy_nulls_NMP <- AdjEnergy(nullNMP_weighted_spect$values)
  } else {
    spect_rad_weighted_nulls_NMP <- INCOMPLETE_MODEL_VALUE
    lpl_weighted_energy_nulls_NMP <- INCOMPLETE_MODEL_VALUE
    adj_weighted_energy_nulls_NMP <- INCOMPLETE_MODEL_VALUE
    nullNMP_weighted_spect <- NULL
  }
    
  
  nested_values <- nested(as.matrix(admatrix), c("NODF","wine"))
  calc_values <- list("model_matrix" = model_matrix,
                      "null_spect" = nullNMP_spect,
                      "lpl_spect_nulls" = lpl_spect_nullNMP,
                      "spect_rad_nulls" = spect_rad_nulls_NMP,
                      "adj_energy_nulls" = adj_energy_nulls_NMP,
                      "lpl_energy_nulls" = lpl_energy_nulls_NMP,
                      "algebraic_connectivity_nulls" = algebraic_connectivity_nulls_NMP,
                      "null_weighted_spect" = nullNMP_weighted_spect,
                      "lpl_weighted_spect_nulls" = lpl_weighted_spect_nullNMP,
                      "spect_rad_weighted_nulls" = spect_rad_weighted_nulls_NMP,
                      "adj_weighted_energy_nulls" = adj_weighted_energy_nulls_NMP,
                      "lpl_weighted_energy_nulls" = lpl_weighted_energy_nulls_NMP,
                      "NODF"=nested_values["NODF"],"wine"=nested_values["wine"])
  return(calc_values)
}

binarize_matrix <- function(m){
  m[m>1]=1
  return(m)
}

sq_adjacency <- function(rmatrix, num_guild_a, num_guild_b)
{
  z_matrix <- matrix(0, ncol = num_guild_b, nrow = num_guild_b)
  incidmatrix <- cbind(rmatrix,z_matrix)
  colnames(incidmatrix)<-NULL
  rownames(incidmatrix)<-NULL
  z_matrix <- matrix(0, ncol = num_guild_a, nrow = num_guild_a)
  incidmatrix <- rbind(cbind(z_matrix,t(incidmatrix[,1:num_guild_a])),incidmatrix)
  binarized_incidmatrix <- incidmatrix
  if (max(incidmatrix)>1)
    binarized_incidmatrix[incidmatrix>1]=1
  return(list(binarized_incidmatrix, incidmatrix))  # They are equal for binary networks
}

store_model_results <- function(dp,num_links,num_nodes,magnitudes,weighted){
  mr <- replicate(length(magnitudes),c())
  names(mr) <- magnitudes
  mr$spect_rad <- dp$null_spect$values[1]
  mr$adj_energy <- AdjEnergy(dp$null_spect$values)
  mr$lpl_energy <- LaplEnergy(dp$lpl_spect_nulls$values,num_links,num_nodes)
  mr$algebraic_connectivity <- dp$lpl_spect_nulls$values[length(dp$lpl_spect_nulls$values)-1]
  if (weighted){
    mr$spect_rad_weighted <- dp$null_weighted_spect$values[1]
    mr$adj_weighted_energy <- AdjEnergy(dp$null_weighted_spect$values)
    mr$lpl_weighted_energy <- LaplweightedEnergy(dp$lpl_weighted_spect_nulls$values,num_links,num_nodes)
  }
  else{
    mr$spect_rad_weighted <- INCOMPLETE_MODEL_VALUE
    mr$adj_weighted_energy <- INCOMPLETE_MODEL_VALUE
    mr$lpl_weighted_energy <- INCOMPLETE_MODEL_VALUE
  }
  return(as.data.frame(mr))
}

gen_null_model <- function(typemodel,resanalysis,nodesga,nodesgb,nlinks,nnodes,magnitudes,weighted)
{
  if (typemodel == "WRND"){
    incidmatrix <- WRG_build(resanalysis$matrix,nodesga,nodesgb,nlinks)
   }
  if (typemodel == "RND"){
    incidmatrix <- erdosrmodel(binarize_matrix(resanalysis$matrix),nodesga,nodesgb,nlinks)
  }
  if (typemodel == "VAZ")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="vaz")[[1]]
  if (typemodel == "BVAZ")
    incidmatrix <- nullmodel(binarize_matrix(resanalysis$matrix),N=1,method="vaz")[[1]]
  if (typemodel == "SWAP")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="swap.web")[[1]]
  # if (typemodel == "BSWAP")
  #   incidmatrix <- nullmodel(binarize_matrix(resanalysis$matrix),N=1,method="swap.web")[[1]]
  if (typemodel == "MGEN")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="mgen")[[1]]
  if (typemodel == "SHUFFLE")
    incidmatrix <- nullmodel(resanalysis$matrix,N=1,method="shuffle.web")[[1]]
  if (typemodel == "BSHUFFLE")
    incidmatrix <- nullmodel(binarize_matrix(resanalysis$matrix),N=1,method="shuffle.web")[[1]]
  if(typemodel == "SYTR"){
    repeat{
    incidmatrix <- SynthTradeNull(resanalysis)$matrix_synth
    if(sum(incidmatrix>0)==resanalysis$links)
      break
    }
  }
  if(typemodel == "WSYTR"){
    # totweight <- sum(resanalysis$matrix)-nlinks
    # repeat{
    #   incidmatrix <- SynthTradeNull(resanalysis)$matrix_synth
    #   if(sum(incidmatrix>0)==resanalysis$links)
    #     break
    # }
    # wmatrix <- unname(incidmatrix)
    # totwmatrix <- sum(wmatrix)
    # for (i in 1:nrow(wmatrix))
    #   for (j in 1:ncol(wmatrix)){
    #     if (wmatrix[i,j]>0)
    #       wmatrix[i,j] <- 1+sum(rbinom(totweight,1,wmatrix[i,j]/totwmatrix))
    #   }
    # incidmatrix <- wmatrix


    repeat{
      incidmatrix <- SynthTradeNull(resanalysis)$matrix_synth
      if(sum(incidmatrix>0)==resanalysis$links)
        break
    }
    totweight <- sum(resanalysis$matrix)
    wmatrix <- unname(incidmatrix)
    totwmatrix <- sum(wmatrix)
    for (i in 1:nrow(wmatrix))
      for (j in 1:ncol(wmatrix)){
        if (wmatrix[i,j]>0)
          wmatrix[i,j] <- max(1,sum(rbinom(totweight,1,wmatrix[i,j]/totwmatrix)))
      }
    incidmatrix <- wmatrix
  }
  
  if (grepl('^DU_', typemodel)){
    fdummy <- as.integer(gsub("DU_","",typemodel))
    incidmatrix <- dummynullmodel_build(nodesga,nodesgb,nlinks,resanalysis,fdummy)
  }
  resp <- null_model_process(incidmatrix,resanalysis,weighted)
  return(list("mres"=store_model_results(resp,nlinks,nnodes,magnitudes,weighted),"resp"=resp,
              "incidmatrix"=incidmatrix))
}

find_model_good <- function(matrix){
  if ( (all(sum(colSums(matrix)==0))) &&
       (all(sum(rowSums(matrix)==0))) )
  {
    model_full <- NULL
    found_model_full <- FALSE
  }  else {
    model_full <- matrix
    found_model_full <- TRUE
  } 
  calc_values=list("model_full"=model_full,"found_model_full"=found_model_full)
  return(calc_values)
}

create_nested_model <- function(na,nb,nlinks,weighted=FALSE,empmatrix=""){
  if (na < nb ){
    nr <- na
    nc <- nb
  } else {
    nr <- nb
    nc <- na
  }
  nl <- nlinks
  nmatrix <-  matrix(rep(0,nc*nr), nrow = nr, ncol = nc)
  nmatrix[1,] <- rep(1,nc)
  rnl <- nlinks - nc
  nmatrix[,1] <- rep(1,nr)
  rnl <- rnl - (nr-1)
  #Without loss of generality we assume that the interaction matrix has a landscape disposition
  height <- round(sqrt(2*rnl*nc/nr))
  for (i in (2:nr)){
    if (rnl==0)
      ntoks = 0
    else{
    ntoks <- max(0,min((height-i+1),rnl))
    rnl <- max(0,rnl-ntoks)
    if ((rnl>0) && (ntoks==0)){
      ntoks <- rnl
      rnl=0
      }
    }
    if (ntoks>0)
      for (j in 1:ntoks)
        nmatrix[i,1+j] <- 1
  }
  nmatrix<-nmatrix[rev(order(rowSums(nmatrix))),]
  if(na<nb)
    nmatrix <- t(nmatrix)
  if (!weighted)
    return(nmatrix)
  # If we want a weighted version, the original weights are distributed
  # according to the degrees.
  mweights <- result_analysis$matrix
  mweights <- rev(sort(mweights[mweights>0]))
  marginal_rows <- rowSums(nmatrix)
  marginal_cols <- colSums(nmatrix)
  importancematrix <- nmatrix
  for (i in 1:nb)
    for(j in 1:na)
      importancematrix[i,j] <- (marginal_rows[i]*marginal_cols[j])*nmatrix[i,j]
  for (w in (mweights)){
    mpos <- which(importancematrix==max(importancematrix),arr.ind = TRUE)
    coln <- mpos[1,2]
    rown <- mpos[1,1]
    nmatrix[rown, coln] <- w
    importancematrix[rown, coln] <- 0
  }
  return(nmatrix)
}

create_models_list <- function(lnames,init_value)
{
  mylist <- replicate(length(lnames),init_value)
  names(mylist) <- lnames
  return(mylist)
}

store_spect_values <- function(val,typem){
  return(data.frame("val"=val,"ind"=seq(1:length(val)),"type" = typem))
}

options( warn = -1 )
datadir <- "data/"

lnetw <- list()

filenames <- Sys.glob(paste0(datadir,"*SD*.csv"))
lnetw <- gsub(datadir,"",filenames)

#lnetw <- c("RA_HP_035.csv")
lnetw <- c("M_SD_008.csv")
seed <- 1000
num_experiments <- 100
ppi <- 300
plottofile <- TRUE
plotzigs <- FALSE  # Plotting ziggurats of al models is rather slow. So when TRUE magnitudes are
                   # not saved. Run the script with a big number of experiments (~1000) to compute
                   # magnitudes and plotzigs FALSE. Run it again selecting just the networks you need
                   # to plot and a small number of experiments (~10)
odir <- "plots/"
rdir <- "results/"
dirnulls <- "nullmatrix/"
NetworkMagsFile <- "NetworkMagnitudes.csv"
dir.create(odir, showWarnings = FALSE)
dir.create(rdir, showWarnings = FALSE)
dir.create(dirnulls, showWarnings = FALSE)
colormodels <- c("NETWORK"="blue","RND"="coral1", "WRND"="coral3", "DU_05"="khaki2",
                 "VAZ"="palegreen3","SHUFFLE"="bisque3","SWAP"="steelblue4","MGEN"="lightblue",
                 "BVAZ"="palegreen","BSHUFFLE"="bisque2",
                 "SYTR"="pink2","WSYTR"="darkorchid1")

fract_dummy <- c(5)
nmagnitudes <- list("spect_rad","adj_energy","lpl_energy")#,"algebraic_connectivity")
for (netw in (lnetw)){
  dfnested <- data.frame("network"=c(),"model"=c(),"NODF"=c(),"wine"=c(), "links"=c(),"totalweight"=c())
  result_analysis <- analyze_network(directory = datadir, netw, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  weighted_network <- sum(result_analysis$matrix > 1)>0
  network_nested_values <- nested(result_analysis$matrix, c("NODF","wine"))
  dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="NETWORK",
                                        "NODF"=network_nested_values["NODF"],
                                        "wine"=network_nested_values["wine"],
                                        "links"=result_analysis$links,"totalweight"=sum(result_analysis$matrix)))
  nodes_a <- result_analysis$num_guild_a
  nodes_b <- result_analysis$num_guild_b
  mtx <- sq_adjacency(result_analysis$matrix, nodes_a, nodes_b)
  adj_sq_matrix <- mtx[[1]]  # binarized matrix
  adj_sq_weighted_matrix <- mtx[[2]] # original matrix, both are equal if the network is binary
  num_links <- result_analysis$links
  num_nodes <- nodes_a+nodes_b
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
  
  #Nested model
  nstmodel <- create_nested_model(nodes_a,nodes_b,num_links)
  nstmtx <- sq_adjacency(nstmodel, nodes_a, nodes_b)
  nstadj_sq_matrix <- nstmtx[[1]]  # binarized matrix
  nstadj_spect <- eigen(nstadj_sq_matrix)
  nstspect_rad <- nstadj_spect$values[1]
  print(sprintf("Nested model spectral radius %.2f",nstspect_rad))

  nstadj_energy <- AdjEnergy(nstadj_spect$values)
  nstlapl_matrix <- 0-nstadj_sq_matrix
  nstsumweights <- rowSums(nstadj_sq_matrix)
  for (i in 1:nrow(nstlapl_matrix))
    nstlapl_matrix[i,i] <- nstsumweights[i]
  nstlpl_spect = eigen(nstlapl_matrix)
  nstlpl_energy <- LaplEnergy(nstlpl_spect$values,
                                              sum(nstmodel>0),num_nodes)
  nnst <- nested(as.matrix(nstmodel), "ALL")
  dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="NESTED",
                                        "NODF"=nnst["NODF"],
                                        "wine"=nnst["wine"],
                                        "links"= sum(nstmodel>0),"totalweight"= sum(nstmodel)))
  
  
  # Weighted magnitudes, only for weighted networks
  if (weighted_network){
    adj_weighted_spect <- eigen(adj_sq_weighted_matrix)
    spect_rad_weighted <- adj_weighted_spect$values[1]
    adj_weighted_energy <- AdjEnergy(adj_weighted_spect$values)
    lapl_weighted_matrix <- 0-adj_sq_weighted_matrix
    sumweights <- rowSums(adj_sq_weighted_matrix)
    for (i in 1:nrow(lapl_weighted_matrix))
      lapl_weighted_matrix[i,i] <- sumweights[i]
    lpl_weighted_spect = eigen(lapl_weighted_matrix)
    lpl_weighted_energy <- LaplweightedEnergy(lpl_weighted_spect$values,num_links,num_nodes)
    print(sprintf("Sum weighted Laplacian spectrum %.2f",sum(lpl_weighted_spect$values)))
    weightednstmodel <- create_nested_model(nodes_a,nodes_b,num_links,weighted = TRUE,empmatrix = result_analysis$matrix)
    weightednstmtx <- sq_adjacency(weightednstmodel, nodes_a, nodes_b)
    weightednstadj_sq_matrix <- weightednstmtx[[2]] 
    weightednstadj_spect <- eigen(weightednstadj_sq_matrix)
    weightednstspect_rad <- weightednstadj_spect$values[1]
    print(sprintf("Weighted nested model spectral radius %.2f",weightednstspect_rad))
    
    weightednstadj_energy <- AdjEnergy(weightednstadj_spect$values)
    weightednstlapl_matrix <- 0-weightednstadj_sq_matrix
    nstsumweights <- rowSums(weightednstadj_sq_matrix)
    for (i in 1:nrow(weightednstlapl_matrix))
      weightednstlapl_matrix[i,i] <- nstsumweights[i]
    weightednstlpl_spect = eigen(weightednstlapl_matrix)
    weightednstlpl_energy <- LaplweightedEnergy(weightednstlpl_spect$values,
                                                      sum(weightednstmodel>0),num_nodes)
    
    weightednnst <- nested(as.matrix(weightednstmodel), "ALL")
    dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="WNESTED",
                                          "NODF"=weightednnst["NODF"],
                                          "wine"=weightednnst["wine"],
                                          "links"= sum(weightednstmodel>0),"totalweight"= sum(weightednstmodel)))
    
  } else {
    lpl_weighted_energy = INCOMPLETE_MODEL_VALUE
  }
  
  if (weighted_network)
    networkspect <- list("spect_rad"=adj_spect$values[1],
                         "adj_energy"=adj_energy,
                       "lpl_energy"=lpl_energy,
                       "algebraic_connectivity"=lpl_spect$values[length(lpl_spect$values)-1],
                       "spect_rad_weighted"=adj_weighted_spect$values[1],
                       "adj_weighted_energy"=adj_weighted_energy,
                       "lpl_weighted_energy"=lpl_weighted_energy)
  else
    networkspect <- list("spect_rad"=adj_spect$values[1],"adj_energy"=adj_energy,
                         "lpl_energy"=lpl_energy,
                         "algebraic_connectivity"=lpl_spect$values[length(lpl_spect$values)-1],
                         "spect_rad_weighted"=INCOMPLETE_MODEL_VALUE,
                         "adj_weighted_energy"=INCOMPLETE_MODEL_VALUE,
                         "lpl_weighted_energy"=INCOMPLETE_MODEL_VALUE)
  
  mnames <- c("RND","MGEN","SHUFFLE","VAZ","SYTR")
  if (weighted_network)
    mnames <- c(mnames,c("SWAP","WSYTR","WRND","BVAZ","BSHUFFLE"))
  nmodels <- length(mnames)
  model_full <- create_models_list(mnames,NULL)
  found_model_full <- create_models_list(mnames,FALSE)
  nname <- gsub(".csv","",netw)
  print(nname)
  nullsinfo <- data.frame("spect_rad"=replicate(num_experiments,INCOMPLETE_MODEL_VALUE))
  nullsinfo$adj_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$algebraic_connectivity <- INCOMPLETE_MODEL_VALUE
  nullsinfo$spect_rad_weighted <- INCOMPLETE_MODEL_VALUE
  nullsinfo$adj_weighted_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_weighted_energy <- INCOMPLETE_MODEL_VALUE
  modresults <- lapply(1:nmodels, function(x) nullsinfo)
  names(modresults) <- mnames
  specresults <- lapply(1:nmodels, function(x) c())
  names(specresults) <- mnames
  incidmatrix <- specresults
  datamod <- stack(networkspect)
  datamod$MODEL <- "NETWORK"
  datamod <- rbind(datamod,data.frame("values"=nstspect_rad,"ind"="spect_rad","MODEL"="NESTED" ))
  datamod <- rbind(datamod,data.frame("values"=nstlpl_energy,"ind"="lpl_energy","MODEL"="NESTED" ))
  datamod <- rbind(datamod,data.frame("values"=nstadj_energy,"ind"="adj_energy","MODEL"="NESTED" ))
  
  if (weighted_network){
    datamod <- rbind(datamod,data.frame("values"=weightednstspect_rad,"ind"="spect_rad_weighted","MODEL"="WNESTED" ))
    datamod <- rbind(datamod,data.frame("values"=weightednstlpl_energy,"ind"="lpl_weighted_energy","MODEL"="WNESTED" ))
    datamod <- rbind(datamod,data.frame("values"=weightednstadj_energy,"ind"="adj_weighted_energy","MODEL"="WNESTED" ))
  }
  nodes_guild_a <- seq(1,result_analysis$num_guild_a)  # guild a top rows/left columns
  nodes_guild_b <- seq(1,result_analysis$num_guild_b)  # guild b bottom rows/right columns
  
  for (k in 1:num_experiments) {
    if (!k%%100)
      print(paste("experiment",k))
    for (tmodel in mnames)
    {
      p <- gen_null_model(tmodel,result_analysis,nodes_guild_a,nodes_guild_b,
                          num_links,num_nodes,nmagnitudes,weighted_network)
      modresults[[tmodel]][k,]<- p$mres
      specresults[[tmodel]] <- p$resp
      incidmatrix[[tmodel]] <- p$incidmatrix
      dfnested <- rbind(dfnested,data.frame("network"=netw,"model"=tmodel,"NODF"=p$resp$NODF,"wine"=p$resp$wine,
                                            "links"= sum(p$incidmatrix>0),"totalweight"= sum(p$incidmatrix)))
                                            
      
      if(!found_model_full[[tmodel]]){
        fmgr <- find_model_good(p$incidmatrix)
        if (fmgr$found_model_full){
          model_full[[tmodel]] <- fmgr$model_full
          found_model_full[[tmodel]] <- TRUE
        }
      }
    }
  }

  
  # Save model results
  for (m in names(modresults)){
    sm <- stack(modresults[[m]])
    sm$MODEL <- m
    datamod <- rbind(datamod,sm)
  }
  
  lheader <- length(networkspect)
  numvals <- (nrow(datamod)-lheader)/num_experiments
 
  ntry <- c()
  for (i in 1:lheader)
    ntry <- c(ntry,1:num_experiments)
  dtry <- replicate(lheader,0)
  for (i in 1:length(names(modresults)))
    dtry <- c(dtry,ntry)
  #datamod$TRY <- dtry
  if (!plotzigs){  
    write.csv(datamod,paste0(rdir,"MODS_",netw),row.names = FALSE)  # Model results
    write.csv(dfnested,paste0(rdir,"NESTED_",netw),row.names = FALSE)  # Model results
  }
  
  testvalues = data.frame("MODEL"=c(),"magnitude"=c(), "networkvalue"=c(),"modelmean"=c(),
                          "quantile"=c(),"distance"=c(),"reldist"=c())
  
  testvalues = data.frame("MODEL"=c(),"magnitude"=c(), "networkvalue"=c(),"modelmean"=c(),
                          "quantile"=c(),"distance"=c(),"reldist"=c())
  
  tmagnitudes <- list("spect_rad","adj_energy","lpl_energy","spect_rad_weighted","adj_weighted_energy","lpl_weighted_energy")
  for (modeln in mnames)
    for (magnitude in tmagnitudes)
    {
    datos <- datamod[(datamod$MODEL==modeln) & (datamod$ind==magnitude),]$values
    mvalue <- mean(datos)
    datored <- datamod[(datamod$MODEL=="NETWORK") & (datamod$ind==magnitude),]$values
    # w <- wilcox.test(x = datos, mu = datored, 
    #                  alternative = "two.sided") 
    percentile <- ecdf(datos)
    testvalues <- rbind(testvalues,data.frame("MODEL"=modeln,
                                              "magnitude"=magnitude, "networkvalue"=datored,
                                              "modelmean"=mvalue,
                                              "quantile"=100*percentile(datored),
                                              "distance"=datored-mvalue,"reldist"=(datored-mvalue)/datored))
    }
  if (!plotzigs)
    write.csv(testvalues,paste0(rdir,"TESTVALUES_",netw),row.names = FALSE)
  # Write Network magnitudes to global results file
  network_data <- data.frame("Network"=nname,"Weighted"=weighted_network,"NodesA"=nodes_a,"NodesB"=nodes_b,"Links"=num_links,
                 "Weighted"=weighted_network,"Model"="NETWORK")
  network_values <- cbind(network_data,as.data.frame(networkspect))
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
  if (plottofile)
  {
    mspec <- mnames
    df_adj_model <- lapply(1:length(mspec), function(x) c())
    names(df_adj_model) <- mspec
    df_lpl_model <- df_adj_model
    df_lpl_weighted_model <- df_adj_model
    
    df_adj_model[["NETWORK"]] <- store_spect_values(adj_spect$values,"NETWORK")
    df_lpl_model[["NETWORK"]] <- store_spect_values(lpl_spect$values,"NETWORK")
    if (weighted_network)
      df_lpl_weighted_model[["NETWORK"]] <- store_spect_values(lpl_weighted_spect$values,"NETWORK")
    for (m in mspec[2:length(mspec)])
    {
      df_adj_model[[m]] <- store_spect_values(specresults[[m]]$null_spect$values,m)
      df_lpl_model[[m]] <- store_spect_values(specresults[[m]]$lpl_spect_null$values,m)
      df_lpl_weighted_model[[m]] <- store_spect_values(specresults[[m]]$lpl_weighted_spect_null$values,m)
    }
    if (plotzigs)
      ziggurat_graph(datadir,netw,plotsdir=odir,print_to_file = plottofile,show_title = FALSE,weighted_links = "log10")
    pl_adj_spectrum <- plot_spectr(do.call(rbind,df_adj_model),title="Adjacency spectrum")
    pl_lpl_spectrum <- plot_spectr(do.call(rbind,df_lpl_model),title="Laplacian spectrum")
    pl_lpl_weighted_spectrum <- plot_spectr(do.call(rbind,df_lpl_weighted_model),title="Weighted Laplacian spectrum")
    
    for (tmodel in mnames)
      if (found_model_full[[tmodel]])
        if (plotzigs)
          plot_null_model(model_full[[tmodel]],tmodel,netw,dirnulls,odir,plottofile)
    if (plotzigs){
      plot_null_model(nstmodel,"NESTED",netw,dirnulls,odir,plottofile)
      if (weighted_network)
        plot_null_model(weightednstmodel,"WNESTED",netw,dirnulls,odir,plottofile)
    }
    save_nested_model(nstmodel,netw,dirnulls,"NESTED")
    if (weighted_network)
      save_nested_model(weightednstmodel,netw,dirnulls,"WNESTED")
    modresultsweighted <- modresults
    for (mname in c("RND","BSHUFFLE","BVAZ","SYTR"))
      modresultsweighted <- remove_model_data(modresultsweighted,mname)
    
    if (weighted_network)
      for (mname in c("SHUFFLE","VAZ","SWAP"))
        modresults <- remove_model_data(modresults,mname)
    else
      for (mname in c("BSHUFFLE","BVAZ"))
        modresults <- remove_model_data(modresults,mname)
    for (mname in c("WRND","WSYTR","MGEN"))
      modresults <- remove_model_data(modresults,mname)

    if (weighted_network){
      pd_weighted_adj_energy <- plot_distributions(modresultsweighted,"adj_weighted_energy","Weighted Energy",adj_weighted_energy,netw,fract_dummy)#,nestedv=weightednstadj_energy)
      pd_weighted_lpl_energy <- plot_distributions(modresultsweighted,"lpl_weighted_energy","Weighted Laplacian energy",lpl_weighted_energy,netw,fract_dummy)#,nestedv=weightednstlpl_energy)
      pd_weighted_spect_rad <- plot_distributions(modresultsweighted,"spect_rad_weighted","Weighted spectral radius",spect_rad_weighted,netw,fract_dummy,nestedv=weightednstspect_rad)
    }
    pd_spect_rad <- plot_distributions(modresults,"spect_rad","Spectral radius",adj_spect$values[1],netw,fract_dummy,nestedv=nstspect_rad)
    pd_lpl_energy <- plot_distributions(modresults,"lpl_energy","Laplacian energy",lpl_energy,netw,fract_dummy) #,nestedv=nstlpl_energy)
    pd_energy <- plot_distributions(modresults,"adj_energy","Energy",adj_energy,netw,fract_dummy)#,nestedv=nstadj_energy)
    if (weighted_network){
    #wsup <- (pl_adj_spectrum | pl_lpl_spectrum | pl_lpl_weighted_spectrum)
    wmed <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage) 
    winf <- (pd_weighted_adj_energy$pimage | pd_weighted_spect_rad$pimage | pd_weighted_lpl_energy$pimage )
    wtot <- (wmed/ winf) + plot_layout(heights = c(0.5,0.5))
             plot_annotation(title = nname,
                      theme = theme(plot.title = element_text(size = 16,hjust=0.5)))
      plwidth <- 13
      plheight <- 10
    } else {
      wsup <- (pl_adj_spectrum | pl_lpl_spectrum)
      winf <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage)
      #wtot <- (wsup / winf) + plot_layout(heights = c(0.5,0.5))
      wtot <- winf
      plot_annotation(title = nname,
                      theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
      plwidth <- 13
      plheight <- 5
    }
    sweight <- "WEIGHTED"
    if (!weighted_network)
      sweight <- "BINARY"
    nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_",sweight))
    png(paste0(nfile,".png"),width=plwidth*ppi,height=plheight*ppi,res=ppi)
    print(wtot)
    dev.off()
    
    
    
    if (weighted_network){
      pd_weighted_adj_energy <- plot_distributions(modresultsweighted,"adj_weighted_energy","Weighted Energy",adj_weighted_energy,netw,fract_dummy,nestedv=weightednstadj_energy)
      pd_weighted_lpl_energy <- plot_distributions(modresultsweighted,"lpl_weighted_energy","Weighted Laplacian energy",lpl_weighted_energy,netw,fract_dummy,nestedv=weightednstlpl_energy)
      pd_weighted_spect_rad <- plot_distributions(modresultsweighted,"spect_rad_weighted","Weighted spectral radius",spect_rad_weighted,netw,fract_dummy,nestedv=weightednstspect_rad)
    }
    pd_spect_rad <- plot_distributions(modresults,"spect_rad","Spectral radius",adj_spect$values[1],netw,fract_dummy,nestedv=nstspect_rad)
    pd_energy <- plot_distributions(modresults,"adj_energy","Energy",adj_energy,netw,fract_dummy,nestedv=nstadj_energy)
    pd_lpl_energy <- plot_distributions(modresults,"lpl_energy","Laplacian energy",lpl_energy,netw,fract_dummy,nestedv=nstlpl_energy)
    if (weighted_network){
      wsup <- (pl_adj_spectrum | pl_lpl_spectrum | pl_lpl_weighted_spectrum)
      wmed <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage) 
      winf <- (pd_weighted_adj_energy$pimage | pd_weighted_spect_rad$pimage | pd_weighted_lpl_energy$pimage )
      wtot <- (wsup / wmed/ winf) + plot_layout(heights = c(0.33,0.33,0.33))
      plot_annotation(title = nname,
                      theme = theme(plot.title = element_text(size = 16,hjust=0.5)))
      plwidth <- 13
      plheight <- 13
    } else {
      wsup <- (pl_adj_spectrum | pl_lpl_spectrum)
      winf <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage)
      wtot <- (wsup / winf) + plot_layout(heights = c(0.5,0.5))
      plot_annotation(title = nname,
                      theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
      plwidth <- 13
      plheight <- 10
    }
    sweight <- "WEIGHTED"
    if (!weighted_network)
      sweight <- "BINARY"
    nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_NESTEDMAGS_",sweight))
    png(paste0(nfile,".png"),width=plwidth*ppi,height=plheight*ppi,res=ppi)
    print(wtot)
    dev.off()
  }
}

