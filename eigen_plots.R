library(ggplot2)
library(kcorebip)
source("manage_null_models.R")
library(patchwork)


mtopower <- function(m,pow){
  mres <- m
  if (pow>1){
    for (i in 2:pow){
      mres <- mres %*% t(m)
    }
  }
  return(mres)
}

create_rnd_matrix <- function(na,np,nlinks){
  
  rndlinks <- matrix(rep(0,na*np),nrow=na,ncol=np)
  rndlinks[sample(seq(1,na*np),nlinks)] <- 1
  randomm <- matrix(rndlinks,nrow=na,ncol=np)
  if (ncol(randomm)>nrow(randomm))
    randomm <- t(randomm)
  return(randomm)
}

create_blocked_model <- function(na,nb,links){
    nr <- nb
    nc <- na

  if (nc*nr==links){  # Full connected network
    mhyp <-  matrix(rep(1,nc*nr), nrow = nr, ncol = nc)
    return(mhyp)
  }
  mhyp <- matrix(rep(0,(nr*nc)),nrow=nr,ncol=nc)
  mhyp[1,1:min(nc,links)]<-1
  if ((links-sum(mhyp))>0){
    munos <- rep(1,links-sum(mhyp))
    munos<-matrix(c(munos,rep(0,(nr-1)*(nc)-sum(munos))),nrow=nr-1,ncol=nc,byrow = T)
    mhyp[2:nr,1:nc]<-munos
  }
  if (ncol(mhyp)<nrow(mhyp))
    mhyp <- t(mhyp)
  return(mhyp)
}


dfnest <- data.frame("species"=c(),"links"=c(),"Lplspectrad"=c())#,"binmatnest"=c())
dfhyper <- dfnest
dfrandom <- dfnest
dfnestbin <- dfnest
na=35 # poner na mayor que np
np=27 
softnest=TRUE
upperlinks = round(na*np)
for (nlinks in seq((na+np-1),upperlinks)){  
  if (nlinks %% 25 == 0)
    print(paste(nlinks,"out of",upperlinks))
  mincidblocked <- create_blocked_model(na,np,nlinks)
  mblocked <- sq_adjacency(mincidblocked, na, np)[[1]]
  mlaplblocked <- create_laplacian_matrix(mblocked)
  
  mincidhyper <- create_hypernested_model(na,np,nlinks)
  if (na>np)
    mhyper <- sq_adjacency(mincidhyper, np, na)[[1]]
  
  else
    mhyper <- sq_adjacency(mincidhyper, na, np)[[1]]
  mlaplhyper <- create_laplacian_matrix(mhyper)
  
  # dfblocked <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(mblocked)$values[1],
  #                                        "binmatnest"=nested(mblocked)[[1]],
  #                                        "AdjEnergy"=AdjEnergy(eigen(mblocked)$values),
  #                                        "Lplspectrad"=eigen(mlaplblocked)$values[1],
  #                                        "LplEnergy"=LaplEnergy(eigen(mlaplblocked)$values,nlinks,na+np)))
  nhyper <- nested(mincidhyper,c("binmatnest","NODF"))
  dfhyper <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(mhyper)$values[1],
                                       "binmatnest"=nhyper[[1]],"NODF"=nhyper[[2]],
                                       "AdjEnergy"=AdjEnergy(eigen(mhyper)$values),
                                       "Lplspectrad"=eigen(mlaplhyper)$values[1],
                                       "LplEnergy"=LaplEnergy(eigen(mlaplhyper)$values,nlinks,na+np)))
                                       #"binmatnest"=nested(mhyper,c("binmatnest"))[[1]]))
                                
  #print(paste("links",nlinks,"Laplacian radius",eigen(mlaplhyper)$values[[1]],"sqrt links",sqrt(nlinks)))
  mincid <- create_perfect_nested_model(na,np,nlinks)
  if (na>np)
    m <- sq_adjacency(mincid,np,na)[[1]]
  else
    m <- sq_adjacency(mincid,na,np)[[1]]
 
  nm <- nested(mincid,c("binmatnest","NODF"))
  dfnest <- rbind(dfnest, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(m)$values[1],
                                     "binmatnest"=nm[[1]],"NODF"=nm[[2]],
                                     "AdjEnergy"=AdjEnergy(eigen(m)$values),
                                    "Lplspectrad"=eigen(mlapl)$values[1],
                                    "LplEnergy"=LaplEnergy(eigen(mlapl)$values,nlinks,na+np)))
                                    #"binmatnest"=nested(m,c("binmatnest"))[[1]]))
  
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
                                         "LplEnergy"=LaplEnergy(eigen(mlaplrandom)$values,nlinks,na+np)))
                                         #"binmatnest"=nested(mrandom,c("binmatnest"))[[1]]))
  
}
dfnest$label="PERFNESTED"
dfhyper$label="HYPERNESTED"
dfrandom$label="RND"
# dfblocked$label="BLOCKED"
dfall <- rbind(dfrandom,dfnest,dfnestbin,dfhyper)

dfnormalized <- dfall
dfnormalized$normspecrad <- 0
for (i in 1:nrow(dfnormalized)){
  minspec <- dfnormalized[(dfnormalized$label=="RND")&(dfnormalized$links==dfnormalized$links[i]),]$spectrad
  maxspec <- dfnormalized[(dfnormalized$label=="HYPERNESTED")&(dfnormalized$links==dfnormalized$links[i]),]$spectrad
  if (minspec != maxspec)
    dfnormalized$normspecrad[i]<-minmaxnorm(dfnormalized$spectrad[i],minspec,maxspec)
}

sr <- ggplot(data=dfall,aes(x=links/(na*np),y=spectrad,color=label))+geom_point(alpha=0.5)+
      geom_point(data=dfall,aes(x=links/(na*np),y=sqrt(links)),color="black",size=0.5)+
      xlab("Connectance")+ggtitle(paste("na:",na,"np:",np))+theme_bw()
linksmed <- 1.0*na*np/2
adje <- ggplot(data=dfall,aes(x=links/(na*np),y=AdjEnergy,color=label))+geom_point(alpha=0.5)+
  #geom_point(data=dfall,aes(x=links,y=linksmed-((links-linksmed)^2/linksmed)),color="black")+
  xlab("Connectance")+ggtitle(paste("na:",na,"np:",np))+theme_bw()

lplsr <- ggplot(data=dfall,aes(x=links/(na*np),y=Lplspectrad,color=label))+geom_point(alpha=0.5)+
        xlab("Connectance")+ggtitle(paste("na:",na,"np:",np))+theme_bw()
lple <- ggplot(data=dfall,aes(x=links/(na*np),y=LplEnergy,color=label))+geom_point(alpha=0.5)+
        xlab("Connectance")+ggtitle(paste("na:",na,"np:",np))+theme_bw()
dirprop <- "PropertiesPlots/"
dir.create(dirprop, showWarnings = FALSE)
ppi <- 100
plsize <- 5
psup <- (sr | adje)
psdown <- (lplsr| lple)
pslinks <- (psup / psdown) + plot_layout(heights = c(0.5,0.5))
png(paste0(dirprop,"LINKPLOTS_",na,"_",np,".png"),width=3*plsize*ppi,height=2*plsize*ppi,res=ppi)
print(pslinks)
dev.off()

EnEn <- ggplot(data=dfall,aes(x=AdjEnergy,y=LplEnergy,color=label))+geom_point()+
        ggtitle(paste("na:",na,"np:",np))+theme_bw()
SrSr <- ggplot(data=dfall,aes(x=spectrad,y=Lplspectrad,color=label))+geom_point()+
        ggtitle(paste("na:",na,"np:",np))+theme_bw()

psAdjLpl <- (SrSr | EnEn)
png(paste0(dirprop,"ADJLPLPLOTS_",na,"_",np,".png"),width=3*plsize*ppi,height=plsize*ppi,res=ppi)
print(psAdjLpl)
dev.off()

bnm <- ggplot(data=dfall,aes(x=links/(na*np),y=binmatnest,color=label))+geom_point(alpha=0.5)+
 xlab("Connectance")+ggtitle(paste("na:",na,"np:",np))+theme_bw()
pnodf <- ggplot(data=dfall,aes(x=links/(na*np),y=NODF,color=label))+geom_point(alpha=0.5)+
 xlab("Connectance")+ggtitle(paste("na:",na,"np:",np))+theme_bw()

pnested <- (bnm | pnodf)

png(paste0(dirprop,"NESTEDNESS_",na,"_",np,".png"),width=4*plsize*ppi,height=1.5*plsize*ppi,res=ppi)
print(pnested)
dev.off()
