rm(list = ls())
library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)
library(ggrepel)
library(ggbeeswarm)
library(ggdist)

source("config_values.R")

library("data.table")
source("colormodels.R")
NESTplot <- function(ndata,yvar,ylabel="",nestmeasure="NODF"){
  if (nestmeasure=="NODF")
    myndata <- data.frame("network"=ndata$network,"NEST"=ndata$NODF,"y"=yvar)
  if (nestmeasure=="wine")
    myndata <- data.frame("network"=ndata$network,"NEST"=ndata$wine,"y"=yvar)
  myndata$interaction <- "SD"
  if (sum(grepl("HP",myndata$network))>0)
    myndata[grepl("HP",myndata$network),]$interaction <- "HP"
  if (sum(grepl("PL",myndata$network))>0)
    myndata[grepl("PL",myndata$network),]$interaction <- "PL"
  nodfp <- ggplot(data=myndata,aes(x=NEST,y=yvar,fill=as.factor(interaction)))+
    geom_point(shape=21,color="transparent")+
    geom_text_repel( aes(x=NEST, y=yvar,color=as.factor(interaction),
                         label=sprintf("%s%d",interaction,as.numeric(gsub('\\D',"",network)))),
                     size=2.5)+ ylab(ylabel)+xlab(nestmeasure)+
    ggtitle(sprintf("Correlation %.2f",cor(myndata$NEST,yvar)))+theme_bw()+
    theme(legend.position = "none")
  return(nodfp)
}

normdistplot <- function(plotdata,disttype,rvalue){
  plotdatanetwork <- plotdata[plotdata$MODEL %in% c("NETWORK","NESTED"),]
  plotdataNULLS <- plotdata[!(plotdata$MODEL %in% c("NETWORK","NESTED")),]
  clist <- get_colors(c("NETWORK","NESTED",unique(plotdataNULLS$MODEL)),colormodels) 
  q <- ggplot(data=plotdataNULLS) +
    geom_jitter(aes( network ,normdist,fill=MODEL),position=position_jitter(0.1),
                color="transparent",shape=21, size=1,  alpha = 0.3)+
    geom_point(data=plotdatanetwork,aes( network ,normdist,fill=MODEL),
               color="transparent",shape=21, size=2,  alpha = 0.8)+
    geom_hline(data = data.frame("val"=rvalue), aes(yintercept = val), 
               color = "black", size=0.75,alpha=0.5, linetype="dotted")+
    geom_hline(data = data.frame("val"=1-rvalue), aes(yintercept = val), 
               color = "red", size=0.75,alpha=0.5)+
    xlab("Network")+ylab("Normalized distance")+
    ggtitle(disttype)+ guides(fill = guide_legend(override.aes = list(size=8,alpha=0.8)))+
    scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
    theme_bw()+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1),
                     panel.grid = element_line(size=0.1, linetype="dotted", colour="darkgrey"))
    return(q)
}

reldistplot <- function(avg_data,disttype){
  clist <- get_colors(unique(avg_data$MODEL),colormodels) 
  s <- ggplot(data=avg_data) +
    geom_point(aes( network ,avg,fill=MODEL),
               color="transparent",shape=21, size=2.5,  alpha = 0.9)+
    scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
    geom_hline(data = data.frame("val"=0), aes(yintercept = val), 
               color = "blue", size=0.5,alpha=0.4)+
    xlab("Network")+ylab("Normalized relative distance null model-network")+
    ggtitle(disttype)+ guides(fill = guide_legend(override.aes = list(size=8,alpha=0.8)))+
    theme_bw()+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1),
    )
  return(s)
}

normandreldistplot <- function(alldistances,ldisttype,lm,lnetworktype,binmags=TRUE){
  for (disttype in ldisttype){
    plotdata <- alldistances[!is.na(alldistances$network),]
    plotdata <- plotdata[plotdata$ind %in% disttype & plotdata$MODEL %in% lm,]
    plotdata$ntype <- "BINARY"
    for (i in 1:nrow(networktypes))
      plotdata[plotdata$network == networktypes$network[i],]$ntype <- networktypes$type[i]
    plotdata <- plotdata[plotdata$ntype %in% lnetworktype,]
    
    if (binmags)
      plotdata <- plotdata[!(plotdata$ntype=="WEIGHTED" & (plotdata$MODEL %in% c("VAZ","SHUFFLE"))),]
    plotdatanor <- plotdata
    if (binmags){
      plotdatanor[plotdatanor$MODEL %in% c("SHUFFLE","BSHUFFLE"),]$MODEL="B/SHUFFLE"
      plotdatanor[plotdatanor$MODEL %in% c("VAZ","BVAZ"),]$MODEL="B/VAZ"
    }
    
    datarel <- plotdatanor[plotdatanor$ind==disttype,]
    datarel <- datarel[datarel$MODEL!="NETWORK",]
    
    p <- datarel %>%
      group_by(MODEL)%>%
      summarise(medn = mean(reldist))
    datarel$MODEL <- factor(datarel$MODEL,levels=p[order(p$medn),]$MODEL)

    avgrels <- datarel %>%
      group_by(network,MODEL)%>%
      summarise(medn = mean(reldist))
    
    
    clist <- get_colors(unique(avgrels$MODEL),colormodels)
    distrpl <- ggplot(data=datarel[datarel$MODEL!="NETWORK"],aes(y=reldist,x=MODEL,fill=MODEL))+
      ggdist::stat_halfeye( adjust = 1, width = 0.4, .width = 0, alpha=0.5,
                            justification = -.5,point_colour = NA ) + 
      geom_boxplot(width = 0.1,alpha = 0,outlier.shape = NA) +
    geom_point(data=avgrels,aes(y=medn,x=MODEL,color=MODEL),
      shape = 95,size = 10, alpha = .2)+
     geom_hline(data = data.frame("val"=0), aes(yintercept = val), 
                 color = "blue", size=0.5,alpha=0.4,linetype="dotted")+
     scale_fill_manual(name='Null model', breaks=clist$lnames, values=clist$lcolors)+
    scale_color_manual(name='Null model', breaks=clist$lnames, values=clist$lcolors)+
      guides(color=guide_legend("MODEL"), fill = "none")+
      ylab(paste("Relative distance",disttype))+theme_bw()+xlab("")+
      theme(legend.position="bottom")#+coord_flip()
    if (binmags){
      w = 9
      distrpl <- distrpl + ggtitle("All networks")
    }
    else{
      w =12
      distrpl <- distrpl + ggtitle("Weighted networks")
    }
    
    
    png(paste0(pldir,"RELDIST_MODELS_",paste(disttype,collapse="-"),".png"),width=w*ppi,height=8*ppi,res=ppi)
    print(distrpl)
    dev.off()
    
    mdata <- plotdata[plotdata$MODEL==lm[1],]
    p <- mdata %>%
      group_by(network)%>%
      summarise(medn = mean(normdist))
    plotdatanor$network <- factor(plotdatanor$network,levels=p[order(p$medn),]$network)
    
    rvalue=0
    if ((disttype=="adj_energy")||(disttype=="adj_weighted_energy"))
      rvalue=1
    ndistall <- normdistplot(plotdatanor,disttype,rvalue)
    png(paste0(pldir,"NORMDIST_DOTS_",paste(disttype,collapse="-"),".png"),width=20*ppi,height=10*ppi,res=ppi)
    print(ndistall)
    dev.off()
    
    if (binmags){
      plotdata[plotdata$MODEL %in% c("SHUFFLE","BSHUFFLE"),]$MODEL="B/SHUFFLE"
      plotdata[plotdata$MODEL %in% c("VAZ","BVAZ"),]$MODEL="B/VAZ"
    }
    else
      plotdata <- plotdata[!(plotdata$MODEL %in% c("BVAZ","BSHUFFLE")),]
    
    avg_data <- plotdata %>%
      group_by(network,MODEL)%>%
      summarise(avg = mean(reldist))
    if (binmags)
      orderdata <- avg_data[avg_data$MODEL=="B/SHUFFLE",]
    else
      orderdata <- avg_data[avg_data$MODEL=="VAZ",]
    avg_data$network <- factor(avg_data$network,orderdata[order(orderdata$avg),]$network)
    s <- reldistplot(avg_data,disttype)
    png(paste0(pldir,"RELDIST_DOT_",paste(disttype,collapse="-"),".png"),width=20*ppi,height=10*ppi,res=ppi)
    print(s)
    dev.off()
  }
}

plot_corr_nor_magnitude <- function(normdists,ndist,ymag,yscale="",cutoff_links=0){
  dp <- normdists[normdists$disttype==ndist,]
  if (cutoff_links != 0)
    dp <- dp[dp$Links<cutoff_links,]
  p <- ggplot(data=dp,aes(x=meannormdist,y=dp[[ymag]]))+geom_point(size=3)+
    xlab("Normalized distance")+ylab(ymag)+xlim(c(min(dp$meannormdist),max(dp$meannormdist)))+
    #scale_x_sqrt()+
    ggtitle(sprintf("Normalized %s Correlation %.2f",ndist,cor(dp[[ymag]],dp$meannormdist,method="spearman")))+theme_bw()
  if(yscale=="log")
    p <- p + scale_y_log10()
  return(p)
}

plot_corr_magnitude <- function(nmags,xmag,ymag,xscale="",yscale="",cutoff_links=0){
  if (cutoff_links != 0)
    nmags <- nmags[nmags$Links<cutoff_links,]
  p <- ggplot(data=nmags,aes(x=nmags[[xmag]],y=nmags[[ymag]]))+
    geom_point(size=3)+
    xlab("Magnitude")+ylab(ymag)+
    ggtitle(sprintf("%s Correlation %.2f",ndist,cor(nmags[[xmag]],nmags[[ymag]],method="spearman")))+theme_bw()
  if(xscale=="log")
    p <- p + scale_x_log10()
  if(yscale=="log")
    p <- p + scale_y_log10()
  return(p)
}

save_corr_plots <- function(ckp,cnp,filetext="",w=18,h=12)
{
  pd <- ( (ckp[[1]] | ckp[[2]] | ckp[[3]]) / (cnp[[1]] | cnp[[2]] | cnp[[3]]) )
  png(paste0(pldir,filetext,".png"),width=w*ppi,height=h*ppi,res=ppi)
  print(pd)
  dev.off()
}

set.seed(122)
for (weightrf in lweightrf)
{
  print(paste("Transformation",weightrf))
  rdir <- paste0(debugpref,dbaseb,weightrf,"/",rdirb)
  resultsdir <- rdir
  normresults <- paste0(rdir,"normalized/")
  odir <- paste0(debugpref,dbaseb,weightrf,"/",odirb)
  networkmags <- fread(paste0(resultsdir,"networkmagnitudes.csv"))
  normdistavgs <- fread(paste0(normresults,"meannormdistances.csv"))
  alldistances <- fread(paste0(normresults,"ALLNORMALIZED.csv"))
  networktypes <- fread(paste0(normresults,"networktypes.csv"))
  nestedmeasures <- fread(paste0(normresults,"ALLNESTEDMEASURES.csv"))
  alldistances$network <- gsub("MINMAX_MODS_","",alldistances$network)
  lbinmags <- c("spect_rad","adj_energy","lpl_energy")
  lweightmags <- c("spect_rad_weighted","adj_weighted_energy","lpl_weighted_energy")
  # Compute relative distance to netork
  alldistances$reldist <- -5
  for (myn in unique(alldistances$network))
    alldistances[alldistances$network==myn,]$reldist <- alldistances[alldistances$network==myn,]$normdist - alldistances[alldistances$network==myn & alldistances$MODEL=="NETWORK",]$normdist
  
  
  # Global boxplots
  
  z <- vector(mode='list', length=6)
  j = 1
  
  # Binarized magnitudes
  for (nt in c("WEIGHTED"))#c("BINARY","WEIGHTED"))
  {
  
    if (nt=="BINARY"){
      #mags <- c("adj_energy","lpl_energy","spect_rad","adj_weighted_energy","lpl_weighted_energy","spect_rad_weighted")
      mags <- lbinmags
      vmodels <- c("NETWORK","HNESTED","NESTED","RND", "SHUFFLE","VAZ")
    }
    else {
      mags <- c(lbinmags,lweightmags)
      vmodels <- c("NETWORK","WNESTED", "VAZ","SHUFFLE","WRND","WSYTR","SWAP","PATEFIELD","MGEN")
      binmags <- lbinmags
      binmodels <- c("NETWORK","HNESTED","NESTED","RND", "BSHUFFLE","VAZ")
    }
    for (ms in mags){
      print(paste(ms,j))
      lmodels = vmodels
      if ((nt=="WEIGHTED") & (ms %in% binmags))
        lmodels = binmodels
      clist <- get_colors(lmodels,colormodels)    
      ndata <- normdistavgs[(normdistavgs$networktype==nt) & (normdistavgs$MODEL %in% lmodels) & (normdistavgs$disttype==ms),]
      p <- ndata %>%
        group_by(MODEL)%>%
        summarise(medn = median(meannormdist))
      ndata$MODEL <- factor(ndata$MODEL,levels=p[order(p$medn),]$MODEL)
      z[[j]] <- ggplot(data=ndata,aes(MODEL,meannormdist,col=MODEL)) +
        #geom_boxplot(fill="blue", alpha = 0.5)+
        geom_beeswarm(alpha = 0.5,pch=18,size=2,corral="omit")+
        geom_hline(yintercept=0, color = "pink",linewidth=1, alpha=0.3)+
        ylim(c(-0.1,1.2))+
        geom_hline(yintercept=1,  color = "pink",linewidth=1, alpha=0.3)+
        scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
        ggtitle(paste(nt,ms))+theme_bw()+ theme(
          legend.position =  "none",
          panel.grid.minor.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major.y = element_line(color="ivory3",linewidth=0.1),
          panel.grid.major.x = element_line(color="ivory3",linewidth=0.1))
      j <- j + 1
    } 
  }
  #allz <- (z[[1]] | z[[2]] | z[[3]]) / (z[[4]] | z[[5]] | z[[6]])
  
  allbinz <- (z[[1]] | z[[2]] | z[[3]])
  
  pldir <- paste0(odir,"/analysis/")
  if (!dir.exists(pldir))
    dir.create(pldir)
  ppi=300
  png(paste0(pldir,"ALLBINZ.png"),width=14*ppi,height=6*ppi,res=ppi)
  print(allbinz)
  dev.off()
  
  allweightz <- (z[[4]] | z[[5]] | z[[6]])
  png(paste0(pldir,"ALLWEIGHTZ.png"),width=14*ppi,height=6*ppi,res=ppi)
  print(allweightz)
  dev.off()
  
  # Nestedness plots
  
  networksnesteddata <- nestedmeasures[nestedmeasures$model=="NETWORK",]
  networksnesteddata$binary_spect_rad <- 0
  networksnesteddata$adj_energy <- 0
  networksnesteddata$lpl_energy <- 0
  #networksnesteddata$algebraic_connectivity <- 0
  for (i in 1:nrow(networksnesteddata)){
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="spect_rad",]
    networksnesteddata$binary_spect_rad[i] <- recordnet$meannormdist
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="adj_energy",]
    networksnesteddata$adj_energy[i] <- recordnet$meannormdist
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="lpl_energy",]
    networksnesteddata$lpl_energy[i] <- recordnet$meannormdist
  }
  binarynets <- networktypes[networktypes$type=="BINARY",]$network
  
  nodfplt_s <- NESTplot(networksnesteddata,networksnesteddata$binary_spect_rad,ylabel="Normalized spectral radius distance")
  nodfplt_adj <- NESTplot(networksnesteddata,networksnesteddata$adj_energy,ylabel="Normalized adjacency energy distance")
  nodfplt_lpl <- NESTplot(networksnesteddata,networksnesteddata$lpl_energy,ylabel="Normalized laplacian energy distance")
  
  nodfplots <- (nodfplt_s | nodfplt_adj | nodfplt_lpl)
  png(paste0(pldir,"NODFPLOTS.png"),width=21*ppi,height=7*ppi,res=ppi)
  
  print(nodfplots)
  dev.off()
  
  weightednets <- networktypes[networktypes$type=="WEIGHTED",]$network
  networksnestedweighted <- networksnesteddata[networksnesteddata$network %in% weightednets]
  
  wineplt_s <- NESTplot(networksnestedweighted,networksnestedweighted$binary_spect_rad,ylabel="Normalized spectral radius distance",nestmeasure = "wine")
  wineplt_adj <-  NESTplot(networksnestedweighted,networksnestedweighted$adj_energy,ylabel="Normalized adjacency energy distance",nestmeasure = "wine")
  wineplt_lpl <- NESTplot(networksnestedweighted,networksnestedweighted$lpl_energy,ylabel="Normalized laplacian energy distance",nestmeasure = "wine")
  
  wineplots <- (wineplt_s | wineplt_adj | wineplt_lpl)
  png(paste0(pldir,"WINEPLOTS.png"),width=21*ppi,height=7*ppi,res=ppi)
  
  print(wineplots)
  dev.off()
  
  # Correlation with size plots
  normdists <- normdistavgs[normdistavgs$MODEL=="NETWORK"]
  nmags <- networkmags[networkmags$Model=="NETWORK",]
  normdists$Links <- 0
  normdists$Nodes <- 0
  normdists$Weight <- 0
  for (i in 1:nrow(normdists)){
    dnet <- nmags[nmags$Network==normdists$network[i],][1]
    normdists$Links[i] <- dnet$Links
    normdists$Nodes[i] <- dnet$NodesA+dnet$NodesB
    normdists$Weight[i] <- dnet$Weight
  }
  nmags$Nodes <- nmags$NodesA + nmags$NodesB
  cnl <- vector(mode='list', length=length(lbinmags))
  cml <- vector(mode='list', length=length(lbinmags))
  cnn <- vector(mode='list', length=length(lbinmags))
  cmn <- vector(mode='list', length=length(lbinmags))
  
  i <- 1
  for (ndist in lbinmags){
    cnl[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Links")
    cml[[i]] <- plot_corr_magnitude(nmags,ndist,"Links")
    cnn[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Nodes")
    cmn[[i]] <- plot_corr_magnitude(nmags,ndist,"Nodes")
    i <- i+1
  }
  save_corr_plots(cml,cnl,filetext="LINKSDISTBIN_LINKS")
  save_corr_plots(cmn,cnn,filetext="LINKSDISTBIN_NODES")
  
  i <- 1
  minlinks <- 300
  for (ndist in lbinmags){
    cnl[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Links",cutoff_links=minlinks)
    cml[[i]] <- plot_corr_magnitude(nmags,ndist,"Links",cutoff_links=minlinks)
    i <- i+1
  }
  save_corr_plots(cml,cnl,filetext=paste0("LINKSDISTBIN_LINKS_CUTOFF_",minlinks))

  i <- 1
  for (ndist in lweightmags){
    cnn[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Weight",yscale="log")
    cmn[[i]] <- plot_corr_magnitude(nmags,ndist,"Weight",xscale="log",yscale="log")
    cnl[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Links",yscale="log")
    cml[[i]] <- plot_corr_magnitude(nmags,ndist,"Links",xscale="log",yscale="log")
    i <- i+1
  }
  save_corr_plots(cml,cnl,filetext="LINKSDISTWEIGHT_LINKS")
  save_corr_plots(cmn,cnn,filetext="LINKSDISTWEIGHT_WEIGHT")
  
  
  # Model comparison plots. All plots
  lm = c("RND","HNESTED")  # The first model will set the order
  lm = c("NETWORK","NESTED","SHUFFLE","VAZ","BVAZ","BSHUFFLE",lm)
  lnetworktype = c("BINARY","WEIGHTED")
  ldisttype <- lbinmags
  normandreldistplot(alldistances,ldisttype,lm,lnetworktype)
  
  # Model comparison, weighted networks
  
  lm = c("WNESTED","WRND")  # The first model will set the order
  lm = c("NETWORK","MGEN","VAZ","PATEFIELD","SWAP","SHUFFLE",lm)
  lnetworktype = c("WEIGHTED")
  ldisttype <- c("spect_rad_weighted","lpl_weighted_energy","adj_weighted_energy")
  normandreldistplot(alldistances,ldisttype,lm,lnetworktype,binmags = FALSE)
  
  # Linear plot
  disttype <- c("lpl_energy","adj_energy")
  plotdata <- alldistances[alldistances$ind %in% disttype,]
  plotdata$ntype <- "BINARY"
  for (i in 1:nrow(networktypes)){
    plotdata[plotdata$network == networktypes$network[i],]$ntype <- networktypes$type[i]
  }
  
  plotdata %>%
         group_by(network,ind,MODEL)%>%
         summarise(avg = mean(normdist))
  disttype <- c("lpl_energy","adj_energy")
  plotdata <- alldistances[alldistances$ind %in% disttype,]
  plotdata$ntype <- ""
  for (i in 1:nrow(networktypes)){
    plotdata[plotdata$network == networktypes$network[i],]$ntype <- networktypes$type[i]
  }
  plotdata <- plotdata[plotdata$MODEL == "NETWORK",]
  plotdata <- plotdata %>%
    group_by(network,ind,MODEL)%>%
    summarise(avg = mean(normdist))
  dx <- plotdata[plotdata$ind == "adj_energy",]
  dy <- plotdata[plotdata$ind == "lpl_energy",]
  
  pld <- data.frame("network"=dx$network,"model"=dx$MODEL,"adj_energy"=dx$avg,
                   "lpl_energy"=dy$avg)
  
  scatterpl <- ggplot(data=pld[pld$model %in% lm,]) +
    geom_point(aes( adj_energy ,lpl_energy,fill=as.factor(model)),
                 color="transparent",shape=21, size=2,  alpha = 0.3)+
    geom_text_repel( aes(x=adj_energy, y=lpl_energy,
                   label=gsub("SD0","SD",gsub("M0","M",gsub("HP0","HP",gsub("PL0","PL",gsub("M","",gsub("RA","",gsub("00","0",gsub("_","",network)))))))),
                   color=as.factor(model)),max.overlaps=6,size=4)+
    geom_hline(yintercept=1,  
               color = "lightblue",linewidth=1, linetype="dotted", alpha=0.8)+
    geom_vline(xintercept=0,  
               color = "lightblue",linewidth=1, linetype="dotted", alpha=0.8)+
    geom_hline(yintercept=0,  
               color = "orange",linewidth=1, linetype="dotted", alpha=0.8)+
    geom_vline(xintercept=1,  
               color = "orange",linewidth=1, linetype="dotted", alpha=0.8)+
    ggtitle(disttype)+ coord_fixed()+
    theme_bw()#+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1))
  
  png(paste0(pldir,"SCATTER_",paste(lm,collapse="-"),"_",paste(disttype,collapse="-"),".png"),width=12*ppi,height=12*ppi,res=ppi)
  print(scatterpl)
  dev.off()
}