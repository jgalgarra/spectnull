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

shortenlabel <- function(network){
 return(gsub("SD0","SD",gsub("M0","M",gsub("HP0","HP",gsub("PL0","PL",gsub("M","",gsub("RA","",gsub("00","0",gsub("_","",network)))))))))
}
NESTplot <- function(ndata,yvar,ylabel="",nestmeasure="NODF"){
  if (nestmeasure=="NODF")
    myndata <- data.frame("network"=ndata$network,"NEST"=ndata$NODF,"y"=yvar)
  if (nestmeasure=="binmatnest.temperature")
    myndata <- data.frame("network"=ndata$network,"NEST"=ndata$binmatnest.temperature,"y"=yvar)
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
               color="transparent",shape=21, size=2.5,  alpha = 0.8)+
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
  avg_data <- avg_data[!(avg_data$MODEL %in% c("VAZ","B/VAZ","SHUFFLE","B/SHUFFLE","SWAP","MGEN","PATEFIELD")),]
  clist <- get_colors(unique(avg_data$MODEL),colormodels) 
  data_network <- avg_data[avg_data$MODEL=="NETWORK",]
  data_network$Connectance <- 0
  for (i in 1:nrow(data_network)){
    data_network$Connectance[i]<-nmags[nmags$Network==data_network$network[i],][1]$Connectance
  }
  data_network$Nodes <- 0
  for (i in 1:nrow(data_network)){
    data_network$Nodes[i]<-nmags[nmags$Network==data_network$network[i],][1]$Nodes
  }
  data_models <- avg_data[avg_data$MODEL!="NETWORK",]
  s <- ggplot(data=data_models) +
    geom_point(aes( network ,avg,fill=MODEL),
               color="transparent",shape=22, size=2.5,  alpha = 0.7)+
    # geom_point(data=data_network,aes( x=network ,y=avg,size=Connectance),
    #            fill= unname(clist$lcolors["NETWORK"]),color="transparent",shape=21, alpha = 0.9)+
    scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
    xlab("Network")+ylab("Normalized relative distance null model-network")+
    ggtitle(disttype)+ guides(fill = guide_legend(override.aes = list(size=8,alpha=0.8)))+
    theme_bw()+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1),
    )
  sc <- s + geom_point(data=data_network,aes( x=network ,y=avg,size=Connectance),
                       fill= unname(clist$lcolors["NETWORK"]),color="transparent",shape=21, alpha = 0.9)
  sn <- s + geom_point(data=data_network,aes( x=network ,y=avg,size=Nodes),
                       fill= unname(clist$lcolors["NETWORK"]),color="transparent",shape=21, alpha = 0.9)
  return(list(sc,sn))
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
      summarise(medn = mean(reldistnetwork))
    datarel$MODEL <- factor(datarel$MODEL,levels=p[order(p$medn),]$MODEL)

    avgrels <- datarel %>%
      group_by(network,MODEL)%>%
      summarise(medn = mean(reldistnetwork))
    
    
    clist <- get_colors(unique(avgrels$MODEL),colormodels)
    distrpl <- ggplot(data=datarel[datarel$MODEL!="NETWORK"],aes(y=reldistnetwork,x=MODEL,fill=MODEL))+
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
    
    png(paste0(pldir,"RELDIST_NETWORK_MODELS_",paste(disttype,collapse="-"),".png"),width=w*ppi,height=8*ppi,res=ppi)
    print(distrpl)
    dev.off()
    
    mdata <- plotdata[plotdata$MODEL==lm[1],]
    p <- mdata %>%
      group_by(network)%>%
      summarise(medn = mean(normdist))
    plotdatanor$network <- factor(plotdatanor$network,levels=p[order(p$medn),]$network)
    plotdatanor$Connectance <- 0
    for (i in 1:nrow(plotdatanor)){
      plotdatanor$Connectance[i]<-nmags[nmags$Network==plotdatanor$network[i],][1]$Connectance
    }
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
      summarise(avg = mean(reldistnetwork))
    if (binmags)
      orderdata <- avg_data[avg_data$MODEL=="B/SHUFFLE",]
    else
      orderdata <- avg_data[avg_data$MODEL=="VAZ",]
    avg_data$network <- factor(avg_data$network,orderdata[order(orderdata$avg),]$network)
    s <- reldistplot(avg_data,disttype)
    png(paste0(pldir,"RELDIST_NETWORK_DOT_",paste(disttype,collapse="-"),"_CONNECTANCE.png"),width=20*ppi,height=10*ppi,res=ppi)
    print(s[[1]])
    dev.off()
    png(paste0(pldir,"RELDIST_NETWORK_DOT_",paste(disttype,collapse="-"),"_NODES.png"),width=20*ppi,height=10*ppi,res=ppi)
    print(s[[2]])
    dev.off()
    
    
    avg_data <- plotdata %>%
      group_by(network,MODEL)%>%
      summarise(avg = mean(reldistnested))
    orderdata <- avg_data[avg_data$MODEL=="NETWORK",]
    avg_data$network <- factor(avg_data$network,orderdata[order(orderdata$avg),]$network)
    t <- reldistplot(avg_data,disttype)
    png(paste0(pldir,"RELDIST_NESTED_DOT_",paste(disttype,collapse="-"),"_CONNECTANCE.png"),width=20*ppi,height=10*ppi,res=ppi)
    print(t[[1]])
    dev.off()
    png(paste0(pldir,"RELDIST_NESTED_DOT_",paste(disttype,collapse="-"),"_NODES.png"),width=20*ppi,height=10*ppi,res=ppi)
    print(t[[2]])
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
    ggtitle(sprintf("Normalized %s Correlation Spearman: %.2f",ndist,
                    cor(dp[[ymag]],dp$meannormdist,method="spearman")))+theme_bw()
  if(yscale=="log")
    p <- p + scale_y_log10()
  
  p <- p+coord_flip()
  return(p)
}

plot_corr_magnitude <- function(nmags,xmag,ymag,xscale="",yscale="",cutoff_links=0){
  if (cutoff_links != 0)
    nmags <- nmags[nmags$Links<cutoff_links,]
  p <- ggplot(data=nmags,aes(x=nmags[[xmag]],y=nmags[[ymag]]))+
    geom_point(size=3)+
    xlab("Magnitude")+ylab(ymag)+
    ggtitle(sprintf("%s Correlation Spearman: %.2f",xmag,
                    cor(nmags[[xmag]],nmags[[ymag]],method="spearman")))+theme_bw()
  if(xscale=="log")
    p <- p + scale_x_log10()
  if(yscale=="log")
    p <- p + scale_y_log10()
  
  p <- p+coord_flip()
  return(p)
}

save_corr_plots <- function(ckp,cnp,filetext="",w=24,h=12,ignore_second_row=FALSE)
{
  if (ignore_second_row)
    pd <- (ckp[[1]] | ckp[[2]] | ckp[[3]] )
  else
    pd <- ( (ckp[[1]] | ckp[[2]] | ckp[[3]] | ckp[[4]])/ (cnp[[1]] | cnp[[2]] | cnp[[3]] | cnp[[4]]) )
  png(paste0(pldir,filetext,".png"),width=w*ppi,height=h*ppi,res=ppi)
  print(pd)
  dev.off()
}

plot_magnitudes_relationship <- function(dfall){
  dfrads <- data.frame("spect_rad"=dfall[dfall$disttype=="spect_rad",]$meannormdist,
                       "lpl_spect_rad"=dfall[dfall$disttype=="lpl_spect_rad",]$meannormdist,
                       "adj_energy"=dfall[dfall$disttype=="adj_energy",]$meannormdist,
                       "lpl_energy"=dfall[dfall$disttype=="lpl_energy",]$meannormdist,
                       "algebraic_connectivity"=dfall[dfall$disttype=="algebraic_connectivity",]$meannormdist,
                       "network"=dfall[dfall$disttype=="spect_rad",]$network)
  plr <- ggplot(data=dfrads,aes(x=spect_rad,y=lpl_spect_rad))+geom_point()+
    geom_text_repel(data=dfrads, aes(x=spect_rad, y=lpl_spect_rad,
                         label=shortenlabel(network)),
                         max.overlaps=6,size=3)+
        ylab("Norm. Lapl. spectral radius")+xlab("Norm. spectral radius")+
             ggtitle(sprintf("Corr Spearman: %.2f ",
                    cor(dfrads$spect_rad,dfrads$lpl_spect_rad,method="spearman")))+theme_bw()
  ple <- ggplot(data=dfrads,aes(x=adj_energy,y=lpl_energy))+geom_point()+
    geom_text_repel(data=dfrads, aes(x=adj_energy, y=lpl_energy,
                                     label=shortenlabel(network)),
                    max.overlaps=6,size=3)+
    ylab("Laplacian Energy")+xlab("Adjacency Energy")+
    ggtitle(sprintf("Corr Spearman: %.2f ",
                    cor(dfrads$adj_energy,dfrads$lpl_energy,method="spearman")))+theme_bw()
  dfalnonnull <- dfrads[dfrads$algebraic_connectivity>0.000001,]
  plac <- ggplot(data=dfalnonnull,aes(x=lpl_spect_rad,y=algebraic_connectivity))+geom_point()+
    geom_text_repel(data=dfalnonnull, aes(x=lpl_spect_rad, y=algebraic_connectivity,
                                     label=shortenlabel(network)),
                    max.overlaps=6,size=3)+
    xlab("Algebraic connectivity")+ylab("Norm. Lpl. spectral radius")+
    ggtitle(sprintf("Corr Spearman: %.2f",
                  cor(dfrads$lpl_spect_rad,dfrads$algebraic_connectivity,method="spearman")))+theme_bw()
  
  
  return(list(plr,plac,ple))
}

create_references_df <-function(distnet,lmagnitudes,refmodel){
  vrefs <- data.frame("dummy"=1)
  for (mag in lmagnitudes){
   vrefs <-cbind(vrefs,(distnet[(distnet$MODEL==refmodel) & (distnet$ind==mag)]$normdist))
    # vrefs <-data.frame("spect_rad"=distnet[(distnet$MODEL=="NETWORK") & (distnet$ind=="spect_rad")]$normdist,
    #        "lpl_spect_rad"=distnet[(distnet$MODEL=="NETWORK") & (distnet$ind=="lpl_spect_rad")]$normdist,
    #        "adj_energy"=distnet[(distnet$MODEL=="NETWORK") & (distnet$ind=="adj_energy")]$normdist,
    #        "lpl_energy"=distnet[(distnet$MODEL=="NETWORK") & (distnet$ind=="lpl_energy")]$normdist)
  }
  vrefs <- vrefs[,2:ncol(vrefs)]
  names(vrefs)<-lmagnitudes
  return(vrefs)
}

plot_hist_spect<-function(dfspectrad,networktype=c("PL","SD","HP"),bcolor="orange")
{  
  dfspectrad$spect_rad <- dfspectrad$meannormdist
  dfspectrad$type <-"NONE"
  for (t in networktype){
    dfspectrad[grepl(paste0("_",t,"_"),dfspectrad$network),]$type <- t
  }
  dfspectrad <- dfspectrad[dfspectrad$type!="NONE",]
  plspectrad <- ggplot(data=dfspectrad,aes(x=spect_rad))+geom_histogram(fill=bcolor,alpha=0.4)+stat_bin(fill="orange",alpha=0.4)+
    geom_vline(xintercept=mean(dfspectrad$spect_rad),  
               color = "black",linewidth=0.75, linetype="dashed", alpha=0.8)+
    theme_bw()+theme(legend.position="none")
  
  return(list(plspectrad,dfspectrad))
} 


set.seed(122)
ppi=300
for (weightrf in lweightrf)
{
  rdir <- paste0(debugpref,dbaseb,weightrf,"/",rdirb)
  resultsdir <- rdir
  normresults <- paste0(rdir,"normalized/")
  odir <- paste0(debugpref,dbaseb,weightrf,"/",odirb)
  pldir <- paste0(odir,"/analysis/")
  if (!dir.exists(pldir))
    dir.create(pldir)
  print(paste("Transformation",weightrf))

  networkmags <- fread(paste0(resultsdir,"networkmagnitudes.csv"))
  normdistavgs <- fread(paste0(normresults,"meannormdistances.csv"))
  alldistances <- fread(paste0(normresults,"ALLNORMALIZED.csv"))
  networktypes <- fread(paste0(normresults,"networktypes.csv"))
  nestedmeasures <- fread(paste0(normresults,"ALLNESTEDMEASURES.csv"))
  dataconnection <- fread(paste0(dataprocessed,"dataconnection.csv"))
  alldistances$network <- gsub("MINMAX_MODS_","",alldistances$network)

  if (ignore_GC_results)    # Use raw networks, even if they are not fully connected
  {
    alldistances <- alldistances[!grepl("_GC",alldistances$network),]
    normdistavgs <- normdistavgs[!grepl("_GC",normdistavgs$network),]
    nestedmeasures <- nestedmeasures[!grepl("_GC",nestedmeasures$network),]
  }
  else                      # Use processed data files with removed disconnected components 
  {
    discnets <- dataconnection[dataconnection$FullyConnected==0,]$Network
    alldistances <- alldistances[!(alldistances$network %in% discnets),]
    normdistavgs <- normdistavgs[!(normdistavgs$network %in% discnets),]
    nestedmeasures <- nestedmeasures[!(nestedmeasures$network %in% discnets),]
  }
  
  lbinmags <- c("spect_rad","adj_energy","lpl_spect_rad","lpl_energy")
  lweightmags <- c("spect_rad_weighted","adj_weighted_energy","lpl_spect_rad_weighted","lpl_weighted_energy")
  
  # Compute relative distance to netork
  alldistances$reldistnetwork <- -5
  alldistances$reldistnested <- -5
  for (myn in unique(alldistances$network)){
    print(myn)
    distnet <- alldistances[alldistances$network==myn,]
    valnetwork <- create_references_df(distnet,lbinmags,"NETWORK")
    valnested <- create_references_df(distnet,lbinmags,"NESTED")
    for (m in names(valnetwork)){
         distnet[distnet$ind==m,]$reldistnetwork <- distnet[distnet$ind==m,]$normdist - valnetwork[[m]]
         distnet[distnet$ind==m,]$reldistnested <- distnet[distnet$ind==m,]$normdist - valnested[[m]]
         
    }
    if (networkmags[networkmags$Network==myn,]$Weighted[1]){
      valnetwork <- create_references_df(distnet,lweightmags,"NETWORK")
      valnested <- create_references_df(distnet,lweightmags,"WNESTED")
      for (m in names(valnetwork)){
        distnet[distnet$ind==m,]$reldistnetwork <- distnet[distnet$ind==m,]$normdist - valnetwork[[m]]
        distnet[distnet$ind==m,]$reldistnested <- distnet[distnet$ind==m,]$normdist - valnested[[m]]
        
      }
    }
    alldistances[alldistances$network==myn,] <- distnet
  }
  
  # Algebraic connectivity histogram
  
  dfalg <- networkmags[networkmags$Model=="NETWORK" & networkmags$algebraic_connectivity>0.000001,]
  plalgconn <- ggplot(data=dfalg,aes(x=algebraic_connectivity,fill = cut(algebraic_connectivity, 100)))+geom_histogram()+stat_bin()+
    geom_vline(xintercept=median(dfalg$algebraic_connectivity),  
               color = "black",linewidth=0.75, linetype="dashed", alpha=0.8)+
    theme_bw()+theme(legend.position="none")
  
  png(paste0(pldir,"ALGCONN_HISTOGRAM.png"),width=5*ppi,height=5*ppi,res=ppi)
  print(plalgconn)
  dev.off()

    
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
  
  allbinz <- (z[[1]] | z[[2]] | z[[3]])
  
  png(paste0(pldir,"ALLBINZ.png"),width=14*ppi,height=6*ppi,res=ppi)
  print(allbinz)
  dev.off()
  
  allweightz <- (z[[4]] | z[[5]] | z[[6]])
  png(paste0(pldir,"ALLWEIGHTZ.png"),width=14*ppi,height=6*ppi,res=ppi)
  print(allweightz)
  dev.off()
  
  # Normalized spectral radius histogram
  dfspectrad <- normdistavgs[normdistavgs$disttype=="spect_rad" & normdistavgs$MODEL=="NETWORK",]
  pllist<- plot_hist_spect(dfspectrad)
  plspectrad <- pllist[[1]]
  png(paste0(pldir,"NORM_SPECTRAD_HISTOGRAM_ALL.png"),width=5*ppi,height=5*ppi,res=ppi)
  print(plspectrad)
  dev.off()
  print("Normalized spectral radius")
  print(ks.test(dfspectrad$meannormdist, "pnorm", mean = median(dfspectrad$meannormdist), sd = sd(dfspectrad$meannormdist)))
  print(shapiro.test(dfspectrad$meannormdist))
  
  pllist<- plot_hist_spect(dfspectrad,networktype = c("PL","SD"),bcolor="green")
  plspectrad <- pllist[[1]]
  dfspectrad <- pllist[[2]]
  png(paste0(pldir,"NORM_SPECTRAD_HISTOGRAM_MUTUALISTC.png"),width=5*ppi,height=5*ppi,res=ppi)
  print(plspectrad)
  dev.off()
  print("Normalized spectral radius mutualistic networks")
  print(ks.test(dfspectrad$meannormdist, "pnorm", mean = median(dfspectrad$meannormdist), sd = sd(dfspectrad$meannormdist)))
  print(shapiro.test(dfspectrad$meannormdist))
  
  
  # Normalized adjacency energy histogram
  dfnetw <- normdistavgs[normdistavgs$disttype=="adj_energy" & normdistavgs$MODEL=="NETWORK",]
  dfnetw$adj_energy <- dfnetw$meannormdist
  pls <- ggplot(data=dfnetw,aes(x=adj_energy))+geom_histogram(fill="pink",alpha=0.4)+stat_bin(fill="pink",alpha=0.4)+
    geom_vline(xintercept=mean(dfnetw$adj_energy),  
               color = "magenta",linewidth=0.75, linetype="dashed", alpha=0.8)+
    theme_bw()+theme(legend.position="none")
  
  png(paste0(pldir,"NORM_ADJ_ENERGY_HISTOGRAM.png"),width=5*ppi,height=5*ppi,res=ppi)
  print(pls)
  dev.off()
  print("Normalized adjacency energy")
  print(ks.test(dfnetw$adj_energy , "pnorm", mean = mean(dfnetw$adj_energy ), sd = sd(dfnetw$adj_energy )))
  print(shapiro.test(dfnetw$adj_energy ))
  
  # Nestedness plots
  
  networksnesteddata <- nestedmeasures[nestedmeasures$model=="NETWORK",]
  networksnesteddata$binary_spect_rad <- 0
  networksnesteddata$lpl_spect_rad <- 0
  networksnesteddata$adj_energy <- 0
  networksnesteddata$lpl_energy <- 0
  #networksnesteddata$algebraic_connectivity <- 0
  for (i in 1:nrow(networksnesteddata)){
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="spect_rad",]
    networksnesteddata$binary_spect_rad[i] <- recordnet$meannormdist
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="lpl_spect_rad",]
    networksnesteddata$lpl_spect_rad[i] <- recordnet$meannormdist
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="adj_energy",]
    networksnesteddata$adj_energy[i] <- recordnet$meannormdist
    recordnet <- normdistavgs[normdistavgs$network==networksnesteddata$network[i] &normdistavgs$MODEL=="NETWORK" &normdistavgs$disttype=="lpl_energy",]
    networksnesteddata$lpl_energy[i] <- recordnet$meannormdist
  }
  binarynets <- networktypes$network
  
  nodfplt_s <- NESTplot(networksnesteddata,networksnesteddata$binary_spect_rad,ylabel="Normalized spectral radius distance")
  nodfplt_lpl_s <- NESTplot(networksnesteddata,networksnesteddata$lpl_spect_rad,ylabel="Normalized Laplacian spectral radius distance")
  nodfplt_adj <- NESTplot(networksnesteddata,networksnesteddata$adj_energy,ylabel="Normalized adjacency energy distance")
  nodfplt_lpl <- NESTplot(networksnesteddata,networksnesteddata$lpl_energy,ylabel="Normalized laplacian energy distance")
  
  nodfplots <- (nodfplt_s | nodfplt_adj) / (nodfplt_lpl_s | nodfplt_lpl)
  png(paste0(pldir,"NODFPLOTS.png"),width=14*ppi,height=14*ppi,res=ppi)
  
  print(nodfplots)
  dev.off()
  
  binmatnestplt_s <- NESTplot(networksnesteddata,networksnesteddata$binary_spect_rad,ylabel="Normalized spectral radius distance",nestmeasure = "binmatnest.temperature")
  binmatnestplt_lpl_s <- NESTplot(networksnesteddata,networksnesteddata$lpl_spect_rad,ylabel="Normalized Laplacian spectral radius distance",nestmeasure = "binmatnest.temperature")
  binmatnestplt_adj <- NESTplot(networksnesteddata,networksnesteddata$adj_energy,ylabel="Normalized adjacency energy distance",nestmeasure = "binmatnest.temperature")
  binmatnestplt_lpl <- NESTplot(networksnesteddata,networksnesteddata$lpl_energy,ylabel="Normalized laplacian energy distance",nestmeasure = "binmatnest.temperature")
  
  binmatnestplots <- (binmatnestplt_s | binmatnestplt_adj) / (binmatnestplt_lpl_s | binmatnestplt_lpl)
  png(paste0(pldir,"BINMMATNESTPLOTS.png"),width=14*ppi,height=14*ppi,res=ppi)
  
  print(binmatnestplots)
  dev.off()
  
  weightednets <- networktypes[networktypes$type=="WEIGHTED",]$network
  networksnestedweighted <- networksnesteddata[networksnesteddata$network %in% weightednets]
  wineplt_s <- NESTplot(networksnestedweighted,networksnestedweighted$binary_spect_rad,ylabel="Normalized spectral radius distance",nestmeasure = "wine")
  wineplt_lpl_s <- NESTplot(networksnestedweighted,networksnestedweighted$lpl_spect_rad,ylabel="Normalized Laplacian spectral radius distance",nestmeasure = "wine")
  wineplt_adj <-  NESTplot(networksnestedweighted,networksnestedweighted$adj_energy,ylabel="Normalized adjacency energy distance",nestmeasure = "wine")
  wineplt_lpl <- NESTplot(networksnestedweighted,networksnestedweighted$lpl_energy,ylabel="Normalized laplacian energy distance",nestmeasure = "wine")
  
  wineplots <- ((wineplt_s | wineplt_adj)/( wineplt_lpl_s| wineplt_lpl))
  png(paste0(pldir,"WINEPLOTS.png"),width=14*ppi,height=14*ppi,res=ppi)
  
  print(wineplots)
  dev.off()
  
  normdists <- normdistavgs[normdistavgs$MODEL=="NETWORK"]
  # Magnitudes plot
  mpl <- plot_magnitudes_relationship(normdists)
  plmag <- (mpl[[1]] | mpl[[2]] | mpl[[3]])
  png(paste0(pldir,"MAGNITUDES_RELATIONSHIP.png"),width=24*ppi,height=8*ppi,res=ppi)
  print(plmag)
  dev.off()
  
  # Correlation with size plots
  normdists <- normdistavgs[normdistavgs$MODEL=="NETWORK"]
  nmags <- networkmags[networkmags$Model=="NETWORK",]
  normdists$Links <- 0
  normdists$Nodes <- 0
  normdists$Weight <- 0
  normdists$Connectance <- 0
  for (i in 1:nrow(normdists)){
    dnet <- nmags[nmags$Network==normdists$network[i],][1]
    normdists$Links[i] <- dnet$Links
    normdists$Nodes[i] <- dnet$NodesA+dnet$NodesB
    normdists$Weight[i] <- dnet$Weight
    normdists$Connectance[i] <- dnet$Connectance
  }
  nmags$Nodes <- nmags$NodesA + nmags$NodesB
  cnl <- vector(mode='list', length=length(lbinmags))
  cml <- vector(mode='list', length=length(lbinmags))
  cnn <- vector(mode='list', length=length(lbinmags))
  cmn <- vector(mode='list', length=length(lbinmags))
  cnc <- vector(mode='list', length=length(lbinmags))
  cmc <- vector(mode='list', length=length(lbinmags))
  
  i <- 1
  for (ndist in lbinmags){
    cnl[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Links")
    cml[[i]] <- plot_corr_magnitude(nmags,ndist,"Links")
    cnn[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Nodes")
    cmn[[i]] <- plot_corr_magnitude(nmags,ndist,"Nodes")
    cnc[[i]] <- plot_corr_nor_magnitude(normdists,ndist,"Connectance")
    cmc[[i]] <- plot_corr_magnitude(nmags,ndist,"Connectance")   
    i <- i+1
  }
  save_corr_plots(cml,cnl,filetext="LINKSDISTBIN_LINKS")
  save_corr_plots(cmn,cnn,filetext="LINKSDISTBIN_NODES")
  save_corr_plots(cmc,cnc,filetext="LINKSDISTBIN_CONNECTANCE")
  
  calgc <- vector(mode='list', length=3)
  algmags <- nmags[nmags$algebraic_connectivity>0.00001,]
  calgc[[1]] <- plot_corr_magnitude(algmags,"algebraic_connectivity","Links")
  calgc[[2]] <- plot_corr_magnitude(algmags,"algebraic_connectivity","Nodes")
  calgc[[3]] <- plot_corr_magnitude(algmags,"algebraic_connectivity","Connectance")
  
  save_corr_plots(calgc,NULL,filetext=paste0("CORR_ALGCONN"),w=24,h=8,ignore_second_row=TRUE)
  
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
                   label=shortenlabel(network),
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