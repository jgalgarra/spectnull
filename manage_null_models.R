create_nullsinfo <- function(num_experiments){
  nullsinfo <- data.frame("spect_rad"=replicate(num_experiments,INCOMPLETE_MODEL_VALUE))
  nullsinfo$adj_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_spect_rad <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$algebraic_connectivity <- INCOMPLETE_MODEL_VALUE
  nullsinfo$spect_rad_weighted <- INCOMPLETE_MODEL_VALUE
  nullsinfo$adj_weighted_energy <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_spect_rad_weighted <- INCOMPLETE_MODEL_VALUE
  nullsinfo$lpl_weighted_energy <- INCOMPLETE_MODEL_VALUE
  return(nullsinfo)
}

create_datamodels <- function(weighted_network,networkspect,pnm,pnms,pnw){
  datamod <- stack(networkspect)
  datamod$MODEL <- "NETWORK"
  datamod <- rbind(datamod,data.frame("values"=pnm$nstspect_rad,"ind"="spect_rad","MODEL"="HNESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnm$nstlplspect_rad,"ind"="lpl_spect_rad","MODEL"="HNESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnm$nstlpl_energy,"ind"="lpl_energy","MODEL"="HNESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnm$nstadj_energy,"ind"="adj_energy","MODEL"="HNESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnm$nstalgebraic_connectivity,"ind"="algebraic_connectivity","MODEL"="HNESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnms$nstspect_rad,"ind"="spect_rad","MODEL"="NESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnms$nstlplspect_rad,"ind"="lpl_spect_rad","MODEL"="NESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnms$nstlpl_energy,"ind"="lpl_energy","MODEL"="NESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnms$nstadj_energy,"ind"="adj_energy","MODEL"="NESTED" ))
  datamod <- rbind(datamod,data.frame("values"=pnms$nstalgebraic_connectivity,"ind"="algebraic_connectivity","MODEL"="NESTED" ))
  
  
  if (weighted_network){
    datamod <- rbind(datamod,data.frame("values"=pnw$weightednstspect_rad,"ind"="spect_rad_weighted","MODEL"="WNESTED" ))
    datamod <- rbind(datamod,data.frame("values"=pnw$weightednstlplspect_rad,"ind"="lpl_spect_rad_weighted","MODEL"="WNESTED" ))
    datamod <- rbind(datamod,data.frame("values"=pnw$weightednstlpl_energy,"ind"="lpl_weighted_energy","MODEL"="WNESTED" ))
    datamod <- rbind(datamod,data.frame("values"=pnw$weightednstadj_energy,"ind"="adj_weighted_energy","MODEL"="WNESTED" ))
  }
return(datamod)
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

distributeweight <- function(nullmatrix,empmatrix){
  ordw <- order(nullmatrix[nullmatrix>0])
  weightsemp <- rev(sort(empmatrix[empmatrix>0]))
  rmatrix <- nullmatrix
  # empweight <- sum(empmatrix)
  # nullweight <- sum(nullmatrix)
  # for (w in weightsemp){
  #   pmax <- which(nullmatrix==max(nullmatrix))[1]
  #   #rmatrix[pmax] <- w*nullmatrix[pmax]
  #   rmatrix[pmax] <- w*empweight/nullweight
  #   
  #   nullmatrix[pmax] <- 0
  # }
  rmatrix <- rmatrix * sum(empmatrix) / sum(rmatrix)
  return(rmatrix)
}

distributeweight1 <- function(nullmatrix,empmatrix){
  meanemp <- mean(empmatrix[empmatrix>0])
  wmatrix <- unname(empmatrix)
  totweight <- sum(empmatrix)
  msynthmatrix <- nullmatrix
  msynthmatrix[msynthmatrix>0] <- 1
  restweight <- totweight - sum(msynthmatrix)
  totwnull <- sum(nullmatrix)
  for (i in 1:nrow(msynthmatrix))
    for (j in 1:ncol(msynthmatrix)){
      if (msynthmatrix[i,j]==1)
        msynthmatrix[i,j] <- 1+sum(rbinom(1,restweight,nullmatrix[i,j]/totwnull))
    }
  return(msynthmatrix)
}

# Weighted Random Model following Garlaschelli's idea for unipartite networks

WRG_build_crazy <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links)
{
  indexesmatrix <- seq(1,nrow(adjmatrix)*ncol(adjmatrix))
  weightmatrix <- sum(adjmatrix)
  tokens <- sample(indexesmatrix,weightmatrix,replace=TRUE)
  nullmatrix <- adjmatrix*0
  for (i in tokens)
    nullmatrix[i] <- nullmatrix[i]+1
  return(nullmatrix)
}


WRG_build_noise <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links)
{
  wmatrix <- sum(adjmatrix)
  rmatrix <- adjmatrix
  rmatrix <- rmatrix*0
  rows <- nrow(adjmatrix)
  cols <- ncol(adjmatrix)
  tokens <- sample(seq(1:(rows*cols)),wmatrix,replace=TRUE)
  for (i in tokens)
    #rmatrix[1+(i%/%rows)][1+(i%/%cols)]<- rmatrix[1+(i%/%rows)][1+(i%/%cols)] + 1
    rmatrix[i] <- rmatrix[i]+1
  return(rmatrix)
}

WRG_build_conserva <- function(adjmatrix,nodes_guild_a,nodes_guild_b,num_links)
{
  binmatrix <- erdosrmodel(adjmatrix,nodes_guild_a,nodes_guild_b,num_links)
  links <- which(binmatrix>0)
  remainingweight <- sum(adjmatrix)-length(links)
  tokens <- sample(links,remainingweight,replace=TRUE)
  for (j in tokens)
    binmatrix[j] <- binmatrix[j]+1
  return(binmatrix)
}

WRGGAR_build <- function(empmatrix,nodes_guild_a,nodes_guild_b,num_links)
{
  W <- sum(empmatrix)
  na <- length(nodes_guild_a)
  nb <- length(nodes_guild_b)
  N <- na*nb
  p <-  2*W/(N+2*W)
  #p <- N/W
  mmatrix <- unname(empmatrix * 0)
  ntrials <- 1.4*W/p
  rtokens <- rbinom(ntrials,1,p)
  sweights <- parsebinoutcome(rtokens)
  sweights <- sweights[1:(nrow(mmatrix)*ncol(mmatrix))]
  mmatrix <- matrix(sweights,nrow=nrow(mmatrix))
  return(mmatrix)
}


dummynullmodel_build <- function(nodes_guild_a,nodes_guild_b,num_links,re,swap_perc)
{
  mmatrix <- re$matrix
  swap_links <- max(1,round(0.01 * swap_perc * num_links))
  ones <- which(result_analysis$matrix > 0)
  zeroes <- which(result_analysis$matrix == 0)
  j = 1
  found <- FALSE
  while ((!found) & (j<10)){
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

null_model_process <- function(incidmatrix,re,weighted_network,num_links,num_nodes)
{
  mtx <- sq_adjacency(incidmatrix,re$num_guild_a,re$num_guild_b)
  model_matrix <- mtx[[1]]  # binarized matrix
  model_weighted_matrix <- mtx[[2]] # original matrix
  nullNMP_spect <- eigen(model_matrix,only.values = TRUE)
  lapl_nullNMP <- 0-model_matrix
  degreesnullNMP <- rowSums(model_matrix)
  for (i in 1:nrow(lapl_nullNMP))
    lapl_nullNMP[i,i] <- degreesnullNMP[i]
  lpl_spect_nullNMP = eigen(lapl_nullNMP,only.values = TRUE)
  if (weighted_network){
    nullNMP_weighted_spect <- eigen(model_weighted_matrix,only.values = TRUE)
    lapl_weighted_nullNMP <- 0-model_weighted_matrix
    nullweightedNMP <- rowSums(model_weighted_matrix)
    for (i in 1:nrow(lapl_weighted_nullNMP))
      lapl_weighted_nullNMP[i,i] <- nullweightedNMP[i]
    lpl_weighted_spect_nullNMP = eigen(lapl_weighted_nullNMP,only.values = TRUE)
  } else {
    lpl_weighted_spect_nullNMP = lpl_spect_nullNMP
  }
  spect_rad_nulls_NMP <- nullNMP_spect$values[1]
  adj_energy_nulls_NMP <- AdjEnergy(nullNMP_spect$values)
  lpl_spect_rad_nulls_NMP <- lpl_spect_nullNMP$values[1]
  lpl_energy_nulls_NMP <- LaplEnergy(lpl_spect_nullNMP$values,num_links,num_nodes)
  algebraic_connectivity_nulls_NMP <- lpl_spect_nullNMP$values[length(lpl_spect_nullNMP$values)-1]
  if (weighted_network){
    spect_rad_weighted_nulls_NMP <- nullNMP_weighted_spect$values[1]
    lpl_spect_rad_weighted_nulls_NMP <- lpl_spect_nullNMP$values[1]
    lpl_weighted_energy_nulls_NMP <- LaplweightedEnergy(lpl_spect_nullNMP$values,num_links,num_nodes)
    adj_weighted_energy_nulls_NMP <- AdjweightedEnergy(nullNMP_weighted_spect$values)
  } 
  
  nested_values <- nested(as.matrix(incidmatrix), c("NODF","wine","binmatnest"))
  if (!weighted_network)
    calc_values <- list("model_matrix" = model_matrix,
                        "null_spect" = nullNMP_spect,
                        "lpl_spect_nulls" = lpl_spect_nullNMP,
                        "spect_rad_nulls" = spect_rad_nulls_NMP,
                        "adj_energy_nulls" = adj_energy_nulls_NMP,
                        "lpl_spect_rad_nulls" = lpl_spect_rad_nulls_NMP,
                        "lpl_energy_nulls" = lpl_energy_nulls_NMP,
                        "algebraic_connectivity_nulls" = algebraic_connectivity_nulls_NMP,
                        "NODF"=nested_values["NODF"],"wine"=nested_values["wine"],
                        "binmatnest.temperature"=nested_values["binmatnest.temperature"])
  else
    calc_values <- list("model_matrix" = model_matrix,
                        "null_spect" = nullNMP_spect,
                        "lpl_spect_nulls" = lpl_spect_nullNMP,
                        "spect_rad_nulls" = spect_rad_nulls_NMP,
                        "lpl_spect_rad_nulls" = lpl_spect_rad_nulls_NMP,
                        "adj_energy_nulls" = adj_energy_nulls_NMP,
                        "lpl_energy_nulls" = lpl_energy_nulls_NMP,
                        "algebraic_connectivity_nulls" = algebraic_connectivity_nulls_NMP,
                        "null_weighted_spect" = nullNMP_weighted_spect,
                        "lpl_weighted_spect_nulls" = lpl_weighted_spect_nullNMP,
                        "spect_rad_weighted_nulls" = spect_rad_weighted_nulls_NMP,
                        "adj_weighted_energy_nulls" = adj_weighted_energy_nulls_NMP,
                        "lpl_spect_rad_weighted_nulls" = lpl_spect_rad_weighted_nulls_NMP,
                        "lpl_weighted_energy_nulls" = lpl_weighted_energy_nulls_NMP,
                        "NODF"=nested_values["NODF"],"wine"=nested_values["wine"],
                        "binmatnest.temperature"=nested_values["binmatnest.temperature"]
                        )
  return(calc_values)
}

binarize_matrix <- function(m){
  m[m>1]=1
  return(m)
}

sq_adjacency <- function(rmatrix, num_guild_a, num_guild_b)#, normalize=TRUE)
{
  zmatrix <- matrix(0, ncol = num_guild_a+num_guild_b, nrow = num_guild_a+num_guild_b)
  zmatrix[1:num_guild_b,(num_guild_b+1):(num_guild_b+num_guild_a)]<-rmatrix
  
  zmatrix[(num_guild_b+1):(num_guild_a+num_guild_b),1:num_guild_b]<-t(rmatrix)
  colnames(zmatrix)<-NULL
  rownames(zmatrix)<-NULL
  binarized_zmatrix <- zmatrix
  if (max(zmatrix)>1)
    binarized_zmatrix[zmatrix>1]=1
  return(list(binarized_zmatrix, zmatrix))  # They are equal for binary networks
}

create_laplacian_matrix <- function(adj_sq_matrix){
  lapl_matrix <- 0-adj_sq_matrix
  weights <- rowSums(adj_sq_matrix)
  for (i in 1:nrow(lapl_matrix))
    lapl_matrix[i,i] <- weights[i]
  return(lapl_matrix)
}


store_model_results <- function(dp,num_links,num_nodes,magnitudes,weighted){
  mr <- replicate(length(magnitudes),c())
  names(mr) <- magnitudes
  mr$spect_rad <- dp$null_spect$values[1]
  mr$adj_energy <- AdjEnergy(dp$null_spect$values)
  mr$lpl_spect_rad <- dp$lpl_spect_nulls$values[1]
  mr$lpl_energy <- LaplEnergy(dp$lpl_spect_nulls$values,num_links,num_nodes)
  mr$algebraic_connectivity <- dp$lpl_spect_nulls$values[length(dp$lpl_spect_nulls$values)-1]
  if (weighted){
    mr$spect_rad_weighted <- dp$null_weighted_spect$values[1]
    mr$adj_weighted_energy <- AdjweightedEnergy(dp$null_weighted_spect$values)
    mr$lpl_spect_rad_weighted <- dp$lpl_weighted_spect_nulls$values[1]
    mr$lpl_weighted_energy <- LaplweightedEnergy(dp$lpl_weighted_spect_nulls$values,num_links,num_nodes)
  }
  else{
    mr$spect_rad_weighted <- INCOMPLETE_MODEL_VALUE
    mr$adj_weighted_energy <- INCOMPLETE_MODEL_VALUE
    mr$lpl_spect_rad_weighted <- INCOMPLETE_MODEL_VALUE
    mr$lpl_weighted_energy <- INCOMPLETE_MODEL_VALUE
  }
  return(as.data.frame(mr))
}

gen_null_model <- function(typemodel,datamod,resanalysis,trfmatrix,nodesga,nodesgb,nlinks,nnodes,magnitudes,weighted)
{
  if (typemodel == "WRND"){
    incidmatrix <- WRG_build_conserva(trfmatrix,nodesga,nodesgb,nlinks)
  }
  if (typemodel == "WRNDGAR"){
    incidmatrix <- WRGGAR_build(trfmatrix,nodesga,nodesgb,nlinks)
  }
  if (typemodel == "RND"){
    incidmatrix <- erdosrmodel(binarize_matrix(trfmatrix),nodesga,nodesgb,nlinks)
  }
  if (typemodel == "VAZ")
    incidmatrix <- bipartite::nullmodel(trfmatrix,N=1,method="vaz")[[1]]
  if (typemodel == "BVAZ")
    incidmatrix <- bipartite::nullmodel(binarize_matrix(trfmatrix),N=1,method="vaz")[[1]]
  if (typemodel == "SWAP")
    incidmatrix <- bipartite::swap.web(trfmatrix,N=1,c.crit=1000)[[1]]
  if (typemodel == "MGEN")
    incidmatrix <- bipartite::nullmodel(trfmatrix,N=1,method="mgen")[[1]]
  if (typemodel == "PATEFIELD")
    incidmatrix <- bipartite::nullmodel(trfmatrix,N=1,method="r2dtable")[[1]]
  if (typemodel == "SHUFFLE")
    incidmatrix <- bipartite::nullmodel(trfmatrix,N=1,method="shuffle.web")[[1]]
  if (typemodel == "BSHUFFLE")
    incidmatrix <- bipartite::nullmodel(binarize_matrix(trfmatrix),N=1,method="shuffle.web")[[1]]
  if(typemodel == "SYTR"){
    repeat{
      incidmatrix <- SynthTradeNull(resanalysis)$matrix_synth
      if(sum(incidmatrix>0)==resanalysis$links)
        break
    }
  }
  if(typemodel == "WSYTR"){
    repeat{
      incidmatrix <- SynthTradeNull(resanalysis)$matrix_synth
      if(resanalysis$links <= 10)  
        break
    }
    incidmatrix <- distributeweight(incidmatrix,trfmatrix) 
    
  }
  
  if (grepl('^DU_', typemodel)){
    fdummy <- as.integer(gsub("DU_","",typemodel))
    incidmatrix <- dummynullmodel_build(nodesga,nodesgb,nlinks,resanalysis,fdummy)
  }
  resp <- null_model_process(incidmatrix,resanalysis,weighted,nlinks,nnodes)
  return(list("mres"=store_model_results(resp,nlinks,nnodes,magnitudes,weighted),"resp"=resp,
              "incidmatrix"=incidmatrix))
}

find_model_fully_connected <- function(matrix){
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



basetriangle <- function(area){
  return((-1 + sqrt(1+8*area))/2)  # Solution of b^2 + b - 2area = 0
}
area_triangle <- function(base){
  return(base*(base+1)/2)
}

create_perfect_nested_model <- function(na,nb,nlinks){
    if (nlinks > na*nb)
      stop("FATAL ERROR. Number of links greater than product of number of species")
    if (na > nb ){
      nr <- nb
      nc <- na
    } else {
      nr <- na
      nc <- nb
    }
    conn <- nlinks/(nr*nc)     # Connectance
    x <- seq(1,nc,by=1)
    k <- (1/conn)*((1-conn)/conn)
    pnr <- rev(nr*(1-(1-((x-1)/(nc-1))^(k*conn))))
    pnr <- pnr+1
    pnr[1] <- pnr[1] - 1
    pnr[2:(length(pnr)-1)]<- round(pnr[2:(length(pnr)-1)]*(nlinks-nr-1)/sum(pnr[2:(length(pnr)-1)]))
    pnr[pnr>nr] <- nr
    prevdifflinks <- 0
    filllinks <- sum(pnr)
    difflinks <- filllinks-nlinks
    if (conn<1)  {
      n = 0
      while ((n<10) && (prevdifflinks!=difflinks)){
        #print(paste("n",n,"nlinks",nlinks,"difflinks",difflinks))
        if (sum(pnr)!=nlinks){
          pnr[2:(length(pnr))] <- round(pnr[2:(length(pnr))]*(nlinks-nr)/sum(pnr[2:(length(pnr))]))

          pnr[pnr>nr] <- nr
          pnr[pnr==0]=1
        }

        n<-n+1
        prevdifflinks <- difflinks
        difflinks <- sum(pnr)-nlinks
      }

    }

    mincid <- matrix(rep(0,nc*nr),nrow=nr,ncol=nc)
    pnr <- rev(sort(pnr))
    for (col in 1:ncol(mincid)){
      mincid[1:pnr[col],col] <- 1
    }
    filllinks <- sum(mincid)
    difflinks <- filllinks-nlinks
    if (difflinks>0){
      removelinks <- difflinks
      remcol <- 2
      while (removelinks > 0){
        remrow <- max(which(mincid[,remcol]==1))
        if (remrow==1){
          remcol <- 2 
          remrow <- max(which(mincid[,remcol]==1))
        }
        mincid[remrow,remcol] <- 0
        removelinks <- removelinks - 1
        remcol <- (remcol + 1)
        if (remcol > ncol(mincid))
          rmcol <- 2
      }
    } else if (difflinks<0){
      addlinks <- -difflinks
      addrow <- 2
      while (addlinks > 0){
        addcol <- max(which(mincid[addrow,]==1))+1
        if (addcol>ncol(mincid)){
          addrow <- addrow + 1
        }
        else {
          mincid[addrow,addcol] <- 1
          addlinks <- addlinks - 1
          addrow <- (addrow + 1)
          if (addrow > nrow(mincid))
            addrow <- 2
        }
      }
    }
    if (sum(abs(rev(sort(colSums(mincid)))-colSums(mincid)))>0)
      mincid<-mincid[,rev(order(colSums(mincid)))]
    if (sum(abs(rev(sort(rowSums(mincid)))-rowSums(mincid)))>0)
      mincid<-mincid[rev(order(rowSums(mincid))),]
    if (ncol(mincid)>nrow(mincid))
      mincid <- t(mincid)
    return(mincid)
  }
  
  
create_bin_nested_model <- function(na,nb,nlinks){
  if (nlinks > na*nb)
    stop("FATAL ERROR. Number of links greater than product of number of species")
  if (na > nb ){
    nr <- nb
    nc <- na
  } else {
    nr <- na
    nc <- nb
  }
  rnl <- nlinks
  if (nc*nr==nlinks){  # Full connected network
    nmatrix <-  matrix(rep(1,nc*nr), nrow = nr, ncol = nc)
    return(nmatrix)
  }
  sqdim <- min(nr,nc)
  nmatrix <-  matrix(rep(0,nc*nr), nrow = nr, ncol = nc)
  nmatrix[1,]<-1
  nmatrix[,1]<-1
  rnl <- rnl - sum(nmatrix)
  minleft <- max(nr,nc)-sum(nmatrix[1,]) # Minimum number of tokens to fill row1
  baset <- floor(basetriangle(rnl-minleft))
  baset <- min(nr-2,baset)
  tokenst <- area_triangle(baset)

  for (i in 1:baset){
    for (j in 1:(baset-i+1))
      if (rnl>0){
        nmatrix[1+i,1+j] <- 1
        rnl <- rnl - 1
      }
  }
  if (minleft>0)
    rowfill <- 1
  else
    rowfill <- 2
  colfill <- which(nmatrix[rowfill,]==0)[1]
  while(rnl>0){
    nmatrix[rowfill,colfill] <- 1
    rnl <- rnl -1
    if (rnl==0)
      break
    # Fill the symmetrical if it's not the last row
    if (colfill <= nr-(rowfill-1)){
      if (nmatrix[colfill,rowfill]==0){
        nmatrix[colfill,rowfill]<-1
        rnl <- rnl - 1
      }
    }
    if (rowfill==nr){
      rowfill <- 2
      while (is.na(which(nmatrix[rowfill,]==0)[1])){
        rowfill <- rowfill + 1
      }
      colfill <- which(nmatrix[rowfill,]==0)[1]
    }
    else if (colfill < nc-(rowfill-1))
      colfill <- colfill + 1
    else{
      rowfill <- rowfill + 1
      colfill <- which(nmatrix[rowfill,]==0)[1]
    }
  }
  if (ncol(nmatrix)>nrow(nmatrix))
    nmatrix <- t(nmatrix)
  return(nmatrix)
}

create_hypernested_model <- function(na,nb,links){
  if (na > nb ){
    nr <- nb
    nc <- na
  } else {
    nr <- na
    nc <- nb
  }
  if (nc*nr==links){  # Full connected network
    mhyp <-  matrix(rep(1,nc*nr), nrow = nr, ncol = nc)
    return(mhyp)
  }
  mhyp <- matrix(rep(0,(nr*nc)),nrow=nr,ncol=nc)
  mhyp[1,]<-1
  mhyp[,1]<-1
  munos <- rep(1,links-sum(mhyp))
  munos<-matrix(c(munos,rep(0,(nr-1)*(nc-1)-sum(munos))),nrow=nr-1,ncol=nc-1,byrow = T)
  mhyp[2:nr,2:nc]<-munos
  if (ncol(mhyp)>nrow(mhyp))
    mhyp=t(mhyp)
  return(mhyp)
}

create_weighted_nested_model <- function(na,nb,nlinks,trfmatrix){
  #nmatrix <- create_bin_nested_model(na,nb,nlinks)
  nmatrix <- create_perfect_nested_model(na,nb,nlinks)
  
  empmatrix <- trfmatrix
  nweight <- sum(empmatrix)
  marginal_rows <- rowSums(nmatrix)/sum(nmatrix)
  marginal_cols <- colSums(nmatrix)/sum(nmatrix)
  wmatrix <- nmatrix
  for (k in 1:nrow(wmatrix))
    wmatrix[k,] <- marginal_cols*marginal_rows[k]
  wmatrix[nmatrix==0]=0
  wmatrix <-wmatrix/sum(wmatrix)
  ordw <- order(wmatrix[wmatrix>0])
  weightsemp <- rev(sort(empmatrix[empmatrix>0]))
  nestedmatrix <- wmatrix
  wweight <- sum(wmatrix)
  for (w in weightsemp){
    pmax <- which(wmatrix==max(wmatrix))[1]
    nestedmatrix[pmax] <- w
    wmatrix[pmax] <- 0
  }
  nestedmatrix <- nestedmatrix * sum(empmatrix) / sum(nestedmatrix)
  return(nestedmatrix)
}

create_weighted_nested_model_binom <- function(na,nb,nlinks){
  #nmatrix <- create_bin_nested_model(na,nb,nlinks)
  nmatrix <- create_perfect_nested_model(na,nb,nlinks)
  
  empmatrix <- result_analysis$matrix
  nweight <- sum(empmatrix)
  marginal_rows <- rowSums(nmatrix)/sum(nmatrix)
  marginal_cols <- colSums(nmatrix)/sum(nmatrix)
  wmatrix <- nmatrix
  wmatrix[nmatrix==0]=0
  wmatrix <-wmatrix/sum(wmatrix)
  for (k in 1:nrow(wmatrix))
    wmatrix[k,] <- marginal_cols*marginal_rows[k]
  #wmatrix[nmatrix==0]=0
  for (i in 1:nrow(wmatrix))
    for (j in 1:ncol(wmatrix))
      wmatrix[i,j] <- rbinom(1,nweight,wmatrix[i,j])
  return(wmatrix)
}

create_normalized_laplacian <- function(empmatrix){
  wrows <- colSums(empmatrix)
  wcols <- rowSums(empmatrix)
  wmatrix <- empmatrix
  for (i in 1:nrow(wmatrix))
    for (j in 1:ncol(wmatrix))
      if (empmatrix[i,j]>0)
        wmatrix[i,j]=-1/sqrt(wrows[i]*wcols[j])
  for (i in 1:nrow(wmatrix))
    if (wrows[i]>0)
      wmatrix[i,i]=1
  return(wmatrix)
}

create_models_list <- function(lnames,init_value)
{
  mylist <- replicate(length(lnames),init_value)
  names(mylist) <- lnames
  return(mylist)
}

remove_model_data <- function(m,mod){
  if (sum(names(m)==mod)==0)
    return(m)
  else
    return(m[-which(names(m)==mod)])
}


process_nested_model_bin <- function(nodes_a,nodes_b,num_links,hypernested=TRUE){
  num_nodes <- nodes_a + nodes_b
  if (hypernested){
    nstmodel <- create_hypernested_model(nodes_a,nodes_b,num_links)
    pref <- "Hyper"
  } else {
    #nstmodel <- create_bin_nested_model(nodes_a,nodes_b,num_links)
    nstmodel <- create_perfect_nested_model(nodes_a,nodes_b,num_links)
    pref <- "Soft"
  }
  if (nodes_a<nodes_b)
    nstadj_sq_matrix <- sq_adjacency(nstmodel, nodes_a, nodes_b)[[1]]  # binarized matrix
  else
    nstadj_sq_matrix <- sq_adjacency(nstmodel, nodes_b, nodes_a)[[1]]  # binarized matrix

  nstadj_spect <- eigen(nstadj_sq_matrix,only.values = TRUE)
  nstspect_rad <- nstadj_spect$values[1]
  secondratio <-  -log10(nstadj_spect$values[2]/nstspect_rad)
  print(sprintf("%s Nested model spectral radius %.2f Second ratio %.2f",pref,nstspect_rad,secondratio))
  nstadj_energy <- AdjEnergy(nstadj_spect$values)
  nstlapl_matrix <- 0-nstadj_sq_matrix
  nstsumweights <- rowSums(nstadj_sq_matrix)
  for (i in 1:nrow(nstlapl_matrix))
    nstlapl_matrix[i,i] <- nstsumweights[i]
  nstlpl_spect = eigen(nstlapl_matrix,only.values = TRUE)
  nstlplspect_rad <- nstlpl_spect$values[1]
  nstalgebraic_connectivity <- sort(nstlpl_spect$values)[2]
  nstlpl_energy <- LaplEnergy(nstlpl_spect$values,sum(nstmodel>0),num_nodes)
  nnst <- nested(as.matrix(nstmodel), c("NODF","wine","binmatnest"))
  calc_values <- list("nstmodel"=nstmodel,"nnst"=nnst,
                      "nstalgebraic_connectivity"=nstalgebraic_connectivity,"nstadj_spect" = nstadj_spect, 
                      "nstspect_rad" = nstspect_rad,"nstadj_energy" = nstadj_energy, 
                      "nstlpl_spect"=nstlpl_spect,"nstlplspect_rad"=nstlplspect_rad, 
                      "nstlpl_energy" = nstlpl_energy)
  return(calc_values)
}

process_nested_model_weighted <- function(nodes_a,nodes_b,num_links,trfmatrix){
  num_nodes <- nodes_a + nodes_b
  weightednstmodel <- create_weighted_nested_model(nodes_a,nodes_b,num_links,trfmatrix)
  if (nodes_a<nodes_b)
    weightednstmstx <- sq_adjacency(weightednstmodel, nodes_a, nodes_b)
  else
    weightednstmstx <- sq_adjacency(weightednstmodel, nodes_b, nodes_a)
  
  weightednstadj_sq_matrix <- weightednstmstx[[2]] 
  weightednstadj_spect <- eigen(weightednstadj_sq_matrix,only.values = TRUE)
  weightednstspect_rad <- weightednstadj_spect$values[1]
  print(sprintf("Weighted nested model spectral radius %.2f",weightednstspect_rad))
  weightednstadj_energy <- AdjweightedEnergy(weightednstadj_spect$values)
  weightednstlapl_matrix <- 0-weightednstadj_sq_matrix
  nstsumweights <- rowSums(weightednstadj_sq_matrix)
  for (i in 1:nrow(weightednstlapl_matrix))
    weightednstlapl_matrix[i,i] <- nstsumweights[i]
  weightednstlpl_spect <- eigen(weightednstlapl_matrix)
  weightednstlplspect_rad <- weightednstlpl_spect$values[1]
  weightednstlpl_energy <- LaplweightedEnergy(weightednstlpl_spect$values,
                                              sum(weightednstmodel>0),num_nodes)
  weightednnst <- nested(as.matrix(weightednstmodel), c("NODF","wine","binmatnest"))
  calc_values <- list("weightednstmodel"=weightednstmodel,"weightednnst"=weightednnst,
                      "weightednstadj_spect" = weightednstadj_spect, "weightednstspect_rad" = weightednstspect_rad,
                      "weightednstadj_energy" = weightednstadj_energy, "weightednstlpl_spect"=weightednstlpl_spect, 
                      "weightednstlplspect_rad"=weightednstlplspect_rad, "weightednstlpl_energy" = weightednstlpl_energy)
  return(calc_values)
}


# Stack null model results: value, magnitude, model
stack_models_results <- function(nnm)
{
  datamod <- nnm$datamod
  for (m in names(nnm$modresults)){
    sm <- stack(nnm$modresults[[m]])
    sm$MODEL <- m
    datamod <- rbind(datamod,sm)
    if (!nnm$weighted_network)
      datamod <- datamod[!grepl("weighted",datamod$ind),]
  }
  return(datamod)
}

save_null_model <- function(netname,dirn,incidmatrix,nullname,adj_spect_val,lpl_spect_val,
                            wadj_spect_val=NULL,wlpl_spect_val=NULL){
  filenull <- paste0(gsub(".csv","",netname),"_NULL_",nullname,".csv")
  write.csv(incidmatrix,paste0(dirnulls,filenull))
  filenullspect <- paste0(gsub(".csv","",netname),"_SPECTRA_",nullname,".csv")
  spectra <- data.frame("adjacency"=adj_spect_val$values,"laplacian"=lpl_spect_val$values)
  spectra$adjacency_weighted <- -1
  spectra$laplacian_weighted <- -1
  if (!is.null(wadj_spect_val))
    spectra$adjacency_weighted <- wadj_spect_val$values
  if (!is.null(wlpl_spect_val))
    spectra$laplacian_weighted <- wlpl_spect_val$values
  spectra[abs(spectra)<0.000000001]<-0
  write.csv(spectra,paste0(dirnulls,filenullspect),row.names=FALSE)
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
                         "lpl_spect_rad"=bm$lpl_spect$values[1],
                         "lpl_energy"=bm$lpl_energy,
                         "algebraic_connectivity"=bm$lpl_spect$values[length(bm$lpl_spect$values)-1],
                         "spect_rad_weighted"=wm$adj_weighted_spect$values[1],
                         "lpl_spect_rad_weighted"=wm$lpl_weighted_spect$values[1],
                         "adj_weighted_energy"=wm$adj_weighted_energy,
                         "lpl_weighted_energy"=wm$lpl_weighted_energy)
  else
    networkspect <- list("spect_rad"=bm$adj_spect$values[1],
                         "adj_energy"=bm$adj_energy,
                         "lpl_spect_rad"=bm$lpl_spect$values[1],
                         "lpl_energy"=bm$lpl_energy,
                         "algebraic_connectivity"=bm$lpl_spect$values[length(bm$lpl_spect$values)-1],
                         "spect_rad_weighted"=-1,
                         "adj_weighted_energy"=-1,
                         "lpl_spect_rad_weighted"=-1,
                         "lpl_weighted_energy"=-1)
  return(networkspect)
}

# Updates Networ magnitudes file
store_network_magnitudes_global_file <- function(netw,nnm){
  # Write Network magnitudes to global results file
  network_data <- data.frame("Network"=gsub(".csv","",netw),"Weighted"=nnm$weighted_network,
                             "NodesA"=nnm$nodes_a,"NodesB"=nnm$nodes_b,"Links"=nnm$num_links,
                             "Weight"=nnm$network_total_weight,"Model"="NETWORK")
  network_values <- cbind(network_data,as.data.frame(nnm$networkspect))
  
  network_values <- rbind(network_values,network_values)
  network_values[2,]$Model = "HNESTED"
  network_values[2,which(names(network_values)=="spect_rad"):length(names(network_values))] = INCOMPLETE_MODEL_VALUE
  network_values <- rbind(network_values,network_values[2,])
  nested_df <- nnm$datamod[nnm$datamod$MODEL=="HNESTED",]
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
  network_values$Connectance <- 0
  for (k in 1:nrow(network_values))
    network_values$Connectance[k] = network_values$Links[k]/(network_values$NodesA[k]*network_values$NodesB[k])
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
                                        "binmatnest.temperature"=network_nested_values["binmatnest.temperature"],
                                        "links"=result_analysis$links,"totalweight"=sum(matrix)))
  
  mtx <- sq_adjacency(matrix, nodes_a, nodes_b)
  adj_sq_matrix <- mtx[[1]]  # binarized matrix
  adj_sq_weighted_matrix <- mtx[[2]] # original matrix, both are equal if the network is binary
  adj_spect <- eigen(adj_sq_matrix)
  adj_energy <- AdjEnergy(adj_spect$values)
  print(sprintf("Sum adjacency spectrum^2 %.2f Energy %.2f",sum((adj_spect$values)^2),adj_energy))
  
  lapl_matrix <- create_laplacian_matrix(adj_sq_matrix)
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
  adj_weighted_spect <- eigen(adj_sq_weighted_matrix,only.values = TRUE)
  spect_rad_weighted <- adj_weighted_spect$values[1]
  adj_weighted_energy <- AdjweightedEnergy(adj_weighted_spect$values)
  lapl_weighted_matrix <- create_laplacian_matrix(adj_sq_weighted_matrix)
  lpl_weighted_spect = eigen(lapl_weighted_matrix,only.values = TRUE)
  lpl_spect_rad_weighted <- lpl_weighted_spect$values[1]
  lpl_weighted_energy <- LaplweightedEnergy(lpl_weighted_spect$values,num_links,num_nodes)
  print(sprintf("Sum weighted Laplacian spectrum %.2f",sum(lpl_weighted_spect$values)))
  calc_values <- list("adj_weighted_spect"=adj_weighted_spect,
                      "adj_weighted_energy"=adj_weighted_energy,
                      "spect_rad_weighted" = spect_rad_weighted,
                      "lpl_weighted_spect"=lpl_weighted_spect,
                      "lpl_spect_rad_weighted" = lpl_spect_rad_weighted,
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
    tmagnitudes <- list("spect_rad","adj_energy","lpl_spect_rad","lpl_energy",
                        "spect_rad_weighted", "adj_weighted_energy","lpl_weighted_energy")
  else
    tmagnitudes <- list("spect_rad","adj_energy","lpl_spect_rad","lpl_energy")
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
  network_nested_values <- nested(trfmatrix, c("NODF","wine","binmatnest"))
  nodes_a <- result_analysis$num_guild_a
  nodes_b <- result_analysis$num_guild_b
  num_links <- result_analysis$links
  num_nodes <- nodes_a+nodes_b
  bmag <- compute_bin_magnitudes(result_analysis,trfmatrix,network_nested_values,nodes_a,nodes_b,num_links,num_nodes,dfnested)
  #Hyper nested model
  pnm <- process_nested_model_bin(nodes_a, nodes_b, num_links, hyper = TRUE)
  dfnested <- rbind(dfnested,bmag$dfnested,
                    data.frame("network"=netw,"model"="HNESTED",
                                                      "NODF"=pnm$nnst["NODF"],
                                                      "wine"=pnm$nnst["wine"],
                                                      "binmatnest.temperature"=pnm$nnst["binmatnest.temperature"],
                                                      "links"= sum(pnm$nstmodel>0),
                                                      "totalweight"= sum(pnm$nstmodel)))
  nmatrix <- pnm$nstmodel
  if (nrow(orimatrix)!=nrow(pnm$nstmodel))
    nmatrix <- t(nmatrix)
  save_null_model(netw,dirnulls,nmatrix,"HNESTED",pnm$nstadj_spect,pnm$nstlpl_spect)
  #Soft nested model
  pnms <- process_nested_model_bin(nodes_a, nodes_b, num_links, hyper = FALSE)
  dfnested <- rbind(dfnested,#bmag$dfnested,
                    data.frame("network"=netw,"model"="NESTED",
                                                      "NODF"=pnms$nnst["NODF"],
                                                      "wine"=pnms$nnst["wine"],
                                                      "binmatnest.temperature"=pnms$nnst["binmatnest.temperature"],
                                                      "links"= sum(pnm$nstmodel>0),
                                                      "totalweight"= sum(pnms$nstmodel)))
  nmatrix <- pnms$nstmodel
  if (nrow(orimatrix)!=nrow(pnms$nstmodel))
    nmatrix <- t(nmatrix)
  save_null_model(netw,dirnulls,nmatrix,"NESTED",pnms$nstadj_spect,pnms$nstlpl_spect)
  
  # Weighted magnitudes, only for weighted networks
  if (weighted_network){
    wmag <- compute_weighted_magnitudes(bmag$adj_sq_weighted_matrix)
    pnw <- process_nested_model_weighted(nodes_a, nodes_b, num_links, trfmatrix)
    dfnested <- rbind(dfnested,data.frame("network"=netw,"model"="WNESTED",
                                          "NODF"=pnw$weightednnst["NODF"],
                                          "wine"=pnw$weightednnst["wine"],
                                          "binmatnest.temperature"=pnm$nnst["binmatnest.temperature"],
                                          "links"= sum(pnw$weightednstmodel>0),"totalweight"= sum(pnw$weightednstmodel)))
    
  } else {
    lpl_weighted_energy = INCOMPLETE_MODEL_VALUE
    wmag <- NULL
    pnw <- NULL
  }
  
  nmatrix <- pnm$nstmodel
  if (nrow(orimatrix)!=nrow(pnm$nstmodel))
    nmatrix <- t(nmatrix)
  
  save_null_model(netw,dirnulls,nmatrix,"WNESTED",pnm$nstadj_spect,pnm$nstlpl_spect,
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
  datamod <- create_datamodels(weighted_network,networkspect,pnm,pnms,pnw)
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
            return(NA)
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
      dfnested <- rbind(dfnested,data.frame("network"=netw,"model"=tmodel,"NODF"=p$resp$NODF,"wine"=p$resp$wine,"binmatnest.temperature"=p$resp$binmatnest.temperature,
                                            "links"= sum(p$incidmatrix>0),"totalweight"= sum(p$incidmatrix)))
    }
  }
  calc_values <- list("matrix"=trfmatrix,"modresults"=modresults,"specresults"=specresults,"datamod"=datamod,
                      "dfnested" = dfnested,"networkspect"=networkspect,
                      "mnames"=mnames,"weighted_network"=weighted_network,
                      "nodes_a"=nodes_a,"nodes_b"=nodes_b,"num_links"=num_links,
                      "network_total_weight"=network_total_weight,
                      "bmag"=bmag,"wmag"=wmag,"pnm"=pnm,"pnms"=pnms,"pnw"=pnw)
  return(calc_values)
}

create_empty_dfnested <- function()
{
  p <- data.frame("network"=c(),"model"=c(),"NODF"=c(),"wine"=c(), "binmatnest.temperature"=c(),
                  "links"=c(),"totalweight"=c())
  return(p)
}

store_nested_values <- function(dfnested){
  dfnested <- rbind(dfnested,data.frame("network"=netw,"model"=tmodel,"NODF"=p$resp$NODF,"wine"=p$resp$wine,
                                        "binmatnest.temperature"=p$resp$binmatnest.temperature,"links"= sum(p$incidmatrix>0),"totalweight"= sum(p$incidmatrix)))
  return(dfnested)
}
debugp<-FALSE
if(debugp){
  na=76
  np=31
  ini<-455
  p <- create_perfect_nested_model(na,np,ini)
  fin <- ini+3
  for (i in seq(ini+1,fin)){
    q<-p
    print(i)
    p <- create_perfect_nested_model(na,np,i)
    #print(p-q)
  }
}
