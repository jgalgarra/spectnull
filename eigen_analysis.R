library(ggplot2)
library(kcorebip)
library("zoo")
source("manage_null_models.R")
source("config_values.R")
library(patchwork)


create_rnd_matrix <- function(na,np,nlinks){
  rndlinks <- matrix(rep(0,na*np),nrow=na,ncol=np)
  rndlinks[sample(seq(1,na*np),nlinks)] <- 1
  randomm <- matrix(rndlinks,nrow=na,ncol=np)
  if (ncol(randomm)>nrow(randomm))
    randomm <- t(randomm)
  return(randomm)
}

prettyfy <- function(pl){
  return(pl+theme_bw()+
        theme(legend.position = "bottom",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 13),
        axis.title = element_text(size = 13)))
}

filenames <- Sys.glob(paste0(datadir,"*_HP*.csv"))
lnetw <- gsub(".csv","",gsub(datadir,"",filenames))
weightrf <- "none"
liminfplconn <- 0.1

create_dirs(weightrf)
                                          
NFile <- paste0(rdir,NetworkMagsFile)
NMags <- read.csv(NFile)
for (netw in lnetw){
  print(netw)
    ndata <- NMags[NMags$Network==netw & NMags$Model=="NETWORK",][1,]
    if (is.na(ndata$Links))
      next    
    dfnest <- data.frame("species"=c(),"links"=c(),"Lplspectrad"=c())
    dfhyper <- dfnest
    dfrandom <- dfnest
    dfnestbin <- dfnest
    na=ndata$NodesA # poner na mayor que np
    np=ndata$NodesB
    networklinks=ndata$Links
    
    upperlinks <- na*np    
    
    steplinks <- (max(1,1+(upperlinks %/% 200)))
    lowerlinks <- (na+np-1)
    
    samplelinks <- seq(lowerlinks,upperlinks,by=steplinks)
    if (!ndata$Links %in% samplelinks){                    # Include the empirical network number of links
      samplelinks <- c(samplelinks,ndata$Links)
      samplelinks <- sort(samplelinks)
    }

    
    for (nlinks in samplelinks){  
      print(paste(nlinks,"out of",upperlinks,"connectance",nlinks/(na*np)))
      
      mincidhyper <- create_hypernested_model(na,np,nlinks)
      if (na>np)
        mhyper <- sq_adjacency(mincidhyper, np, na)[[1]]
      
      else
        mhyper <- sq_adjacency(mincidhyper, na, np)[[1]]
      mlaplhyper <- create_laplacian_matrix(mhyper)
      
      nhyper <- nested(mincidhyper,c("binmatnest","NODF"))
      dfhyper <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(mhyper)$values[1],
                                           "binmatnest"=nhyper[[1]],"NODF"=nhyper[[2]],
                                           "AdjEnergy"=AdjEnergy(eigen(mhyper)$values),
                                           "Lplspectrad"=eigen(mlaplhyper)$values[1],
                                           "AlgConnectivity"=rev(eigen(mlaplhyper)$values)[2],
                                           "LplEnergy"=LaplEnergy(eigen(mlaplhyper)$values,nlinks,na+np)))

      mincid <- create_perfect_nested_model(na,np,nlinks)
      if (na>np)
        m <- sq_adjacency(mincid,np,na)[[1]]
      else
        m <- sq_adjacency(mincid,na,np)[[1]]
      mlapl <- create_laplacian_matrix(m)
      nm <- nested(mincid,c("binmatnest","NODF"))
      dfnest <- rbind(dfnest, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(m)$values[1],
                                         "binmatnest"=nm[[1]],"NODF"=nm[[2]],
                                         "AdjEnergy"=AdjEnergy(eigen(m)$values),
                                        "Lplspectrad"=eigen(mlapl)$values[1],
                                        "AlgConnectivity"=rev(eigen(mlapl)$values)[2],
                                        "LplEnergy"=LaplEnergy(eigen(mlapl)$values,nlinks,na+np)))
    
      mincidrandom <- create_rnd_matrix(na,np,nlinks)
      if (na>np)
        mrandom <- sq_adjacency(mincidrandom,np,na)[[1]]
      else
        mrandom <- sq_adjacency(mincidrandom,na,np)[[1]]
      mlaplrandom <- create_laplacian_matrix(mrandom)
      nmrandom <- nested(mincidrandom,c("binmatnest","NODF"))
      dfrandom <- rbind(dfrandom, data.frame("species"=na+np,"links"=nlinks,
                                             "binmatnest"=nmrandom[[1]],"NODF"=nmrandom[[2]],
                                             "spectrad"=eigen(mrandom)$values[1],
                                             "AdjEnergy"=AdjEnergy(eigen(mrandom)$values),
                                             "Lplspectrad"=eigen(mlaplrandom)$values[[1]],
                                             "AlgConnectivity"=rev(eigen(mlaplrandom)$values)[2],
                                             "LplEnergy"=LaplEnergy(eigen(mlaplrandom)$values,nlinks,na+np)))
                                             
    }
    dfnest$Model="PERFNESTED"
    dfhyper$Model="HYPERNESTED"
    dfrandom$Model="RND"
    # dfblocked$Model="BLOCKED"
    dfall <- rbind(dfhyper,dfnest,dfrandom)
    
    dfrand <- dfall[dfall$Model=="RND",]
    dfhyper <- dfall[dfall$Model=="HYPERNESTED",]
    
    lrsr <- lm(dfrand$spectrad ~ dfrand$links)
    lrlplsr <- lm(dfrand$Lplspectrad ~ dfrand$links)
    lrlplsrhyp <- lm(dfhyper$Lplspectrad ~ dfhyper$links)
    lr_lplsr_sr <- lm(dfrand$Lplspectrad ~ dfrand$spectrad)
    predsrrnd <- unname(predict(lrsr, data.frame(dfrand$links)))
    predlplsrrnd <- unname(predict(lrlplsr, data.frame(dfrand$links)))
    predlplsrhyper <- unname(predict(lrlplsrhyp, data.frame(dfhyper$links)))
    
    dfall$predsrrnd <- 0
    dfall[dfall$Model=="HYPERNESTED",]$predsrrnd <- predsrrnd 
    dfall[dfall$Model=="PERFNESTED",]$predsrrnd <- predsrrnd 
    dfall[dfall$Model=="RND",]$predsrrnd <- predsrrnd 
    dfall$predlplsrrnd <- 0
    dfall[dfall$Model=="HYPERNESTED",]$predlplsrrnd <- predlplsrrnd 
    dfall[dfall$Model=="PERFNESTED",]$predlplsrrnd <- predlplsrrnd 
    dfall[dfall$Model=="RND",]$predlplsrrnd <- predlplsrrnd 
    dfall$predlplsrhyper <- 0
    dfall[dfall$Model=="HYPERNESTED",]$predlplsrhyper <- predlplsrhyper 
    dfall[dfall$Model=="PERFNESTED",]$predlplsrhyper <- predlplsrhyper 
    dfall[dfall$Model=="RND",]$predlplsrhyper <- predlplsrhyper 
    
    
    print("Lpl RDN")
    print(summary(lrlplsr))
    print("Lpl hyper")
    print(summary(lrlplsrhyp))
    
    dfnormalized <- dfall
    dfnormalized$normspecrad <- 0
    for (i in 1:nrow(dfnormalized)){
      minspec <- dfnormalized[(dfnormalized$Model=="RND")&(dfnormalized$links==dfnormalized$links[i]),]$spectrad
      maxspec <- dfnormalized[(dfnormalized$Model=="HYPERNESTED")&(dfnormalized$links==dfnormalized$links[i]),]$spectrad
      if (minspec != maxspec)
        dfnormalized$normspecrad[i]<-minmaxnorm(dfnormalized$spectrad[i],minspec,maxspec)
    }
    sr <- ggplot(data=dfall,aes(x=links/(na*np),y=spectrad,color=Model))+geom_point(alpha=0.5)+
          geom_point(data=dfall,aes(x=links/(na*np),y=sqrt(links)),color="blue",size=0.5)+
      geom_point(data=dfall,aes(x=links/(na*np),y=predsrrnd),color="magenta",size=0.5)+
          geom_point(data=ndata,aes(x=ndata$Connectance,y=ndata$spect_rad),color="black",size=2.5)+
          xlab("Connectance")+ggtitle(paste("NA:",na,"NB:",np))+theme_bw()

    adje <- ggplot(data=dfall,aes(x=links/(na*np),y=AdjEnergy,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=ndata,aes(x=ndata$Connectance,y=ndata$adj_energy),color="black",size=2.5)+
      xlab("Connectance")+ggtitle(paste("NA:",na,"NB:",np))+theme_bw()
    
    lplsr <- ggplot(data=dfall,aes(x=links/(na*np),y=Lplspectrad,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=ndata,aes(x=ndata$Connectance,y=ndata$lpl_spect_rad),color="black",size=2.5)+
      geom_point(data=dfall,aes(x=links/(na*np),y=predlplsrrnd),color="magenta",size=0.5)+
      geom_point(data=dfall,aes(x=links/(na*np),y=predlplsrhyper),color="blue",size=0.5)+
            xlab("Connectance")+ggtitle(paste("NA:",na,"NB:",np))+theme_bw()
    lple <- ggplot(data=dfall,aes(x=links/(na*np),y=LplEnergy,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=ndata,aes(x=ndata$Connectance,y=ndata$lpl_energy),color="black",size=2.5)+
            xlab("Connectance")+ggtitle(paste("NA:",na,"NB:",np))+theme_bw()
    dirprop <- "PropertiesPlots/"
    dir.create(dirprop, showWarnings = FALSE)
    ppi <- 100
    plsize <- 5
    psup <- (prettyfy(sr) | prettyfy(adje))
    psdown <- (prettyfy(lplsr)| prettyfy(lple))
    pslinks <- (psup / psdown) + plot_layout(heights = c(0.5,0.5))
    png(paste0(eiplotsdir,netw,"_","LINKPLOTS_",na,"_",np,".png"),width=3*plsize*ppi,height=3*plsize*ppi,res=ppi)
    print(pslinks)
    dev.off()
    EnEn <- ggplot(data=dfall,aes(x=AdjEnergy,y=LplEnergy,color=Model))+geom_point(alpha=0.5)+
            geom_point(data=ndata,aes(x=ndata$adj_energy,y=lpl_energy),color="black",size=2.5)+
            ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))
    SrSr <- ggplot(data=dfall,aes(x=spectrad,y=Lplspectrad,color=Model))+geom_point(alpha=0.5)+
            geom_point(data=ndata,aes(x=spect_rad,y=lpl_spect_rad),color="black",size=2.5)+
            ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))
    EnEn <- prettyfy(EnEn)
    SrSr <- prettyfy(SrSr)
    psAdjLpl <- (SrSr | EnEn)
    png(paste0(eiplotsdir,netw,"_","ADJLPLPLOTS_",na,"_",np,".png"),width=3*plsize*ppi,height=1.5*plsize*ppi,res=ppi)
    print(psAdjLpl)
    dev.off()
    datanestnetwork <- read.csv(paste0(rdir,"NESTED_",netw,".csv"))
    datanetwork <- read.csv(paste0(rdir,"MODS_",netw,".csv"))
    datanestnetwork <- datanestnetwork[datanestnetwork$model=="NETWORK",]
    datanestnetwork$connectance <- ndata$Connectance
    datanestnetwork$spectrad <- ndata$spect_rad
    datanestnetwork$algconn <-datanetwork[datanetwork$MODEL=="NETWORK" & datanetwork$ind=="algebraic_connectivity",]$values
    dfall$normindex <- 0
    for (k in 1:nrow(dfall)){
      vmax <- dfall[dfall$Model =="HYPERNESTED" & dfall$links == dfall$links[k],]$spectrad
      vmin <- dfall$predsrrnd[k]
      dfall$normindex[k] = (dfall$spectrad[k]-vmin)/(vmax - vmin)
    }
    dfall$lplnormindex <- 0
    for (k in 1:nrow(dfall)){
      vmax <- dfall[dfall$Model =="HYPERNESTED" & dfall$links == dfall$links[k],]$Lplspectrad
      vmin <- dfall$predlplsrrnd[k]
      dfall$lplnormindex[k] = (dfall$Lplspectrad[k]-vmin)/(vmax - vmin)
    }
    
    nmin <- dfall[dfall$Model =="RND" & dfall$links == ndata$Links,]$spectrad
    nmax <- dfall[dfall$Model =="HYPERNESTED" & dfall$links == ndata$Links,]$spectrad
    datanestnetwork$normindex <- (ndata$spect_rad - nmin) /
                                 (nmax- nmin ) 
    
    nmin <- dfall[dfall$Model =="RND" & dfall$links == ndata$Links,]$Lplspectrad
    nmax <- dfall[dfall$Model =="HYPERNESTED" & dfall$links == ndata$Links,]$Lplspectrad
    datanestnetwork$lplnormindex <- (ndata$lpl_spect_rad - nmin) / (nmax- nmin ) 
    
    bnm <- ggplot(data=dfall,aes(x=links/(na*np),y=binmatnest,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=dfall,aes(x=ndata$Connectance,y=datanestnetwork$binmatnest.temperature),color="black",size=2.5)+
      xlim(c(min(liminfplconn,1-as.numeric(ndata$Connectance<liminfplconn)),1))+xlab("Connectance")+ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))+theme_bw()
    pnodf <- ggplot(data=dfall,aes(x=links/(na*np),y=NODF,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=datanestnetwork,aes(x=connectance,y=NODF),color="black",size=2.5)+
      xlim(c(min(liminfplconn,1-as.numeric(ndata$Connectance<liminfplconn)),1))+xlab("Connectance")+ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))+theme_bw()
    pnormindex <- ggplot(data=dfall[dfall$links<0.95*na*np,],aes(x=links/(na*np),y=normindex,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=datanestnetwork,aes(x=connectance,y=normindex),color="black",size=2.5)+
      xlim(c(min(liminfplconn,1-as.numeric(ndata$Connectance<liminfplconn)),1))+xlab("Connectance")+ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))+theme_bw()
    pnormlpl <- ggplot(data=dfall[dfall$links<0.95*na*np,],aes(x=links/(na*np),y=lplnormindex,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=datanestnetwork,aes(x=connectance,y=lplnormindex),color="black",size=2.5)+
      xlim(c(min(liminfplconn,1-as.numeric(ndata$Connectance<liminfplconn)),1))+xlab("Connectance")+ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))+theme_bw()
    
    # psr <- ggplot(data=dfall,aes(x=links/(na*np),y=spectrad,color=Model))+geom_point(alpha=0.5)+
    #   geom_point(data=datanestnetwork,aes(x=connectance,y=spectrad),color="black",size=2.5)+
    #   xlim(c(min(liminfplconn,1-as.numeric(ndata$Connectance<liminfplconn)),1))+xlab("Connectance")+ggtitle(sprintf("NA: %d NB: %d conn: %.2f",na,np,ndata$Connectance))+theme_bw()
    
    bnm <- prettyfy(bnm)
    pnodf <- prettyfy(pnodf)
    pnormindex <- prettyfy(pnormindex)
    pnormlpl <- prettyfy(pnormlpl)
    # psr <- prettyfy(psr)
    pup <- (bnm | pnodf)
    pdown <- (pnormindex | pnormlpl)
    pnested <- (pup / pdown)
    png(paste0(eiplotsdir,netw,"_","NESTEDNESS_",na,"_",np,".png"),width=4*plsize*ppi,height=3*plsize*ppi,res=ppi)
    print(pnested)
    dev.off()
    
    algconn <- ggplot(data=dfall[dfall$links<=max(ceiling(datanestnetwork$connectance*max(dfall$links)),
                                                  max(dfall$links)/2),],
                                                  aes(x=links/(na*np),y=AlgConnectivity,color=Model))+geom_point(alpha=0.5)+
      geom_point(data=datanestnetwork,aes(x=connectance,y=algconn),color="black",size=2.5)+
      xlab("Connectance")+ggtitle(paste("NA:",na,"NB:",np))+theme_bw()
    
    plsize <- 5
    psup <- (prettyfy(algconn))
    png(paste0(eiplotsdir,netw,"_","ALGCONNECTIVITY_",na,"_",np,".png"),width=2*plsize*ppi,height=2*plsize*ppi,res=ppi)
    print(psup)
    dev.off()
}