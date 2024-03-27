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
  if (na<np)
    randomm=t(randomm)
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

  return(mhyp)
}



dfnest <- data.frame("species"=c(),"links"=c(),"Lplspectrad"=c())
dfhyper <- dfnest
dfrandom <- dfnest
na=50 # poner na mayor que np
np=50                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
frows=5
softnest=TRUE
#for (nlinks in seq((na+np-1)+frows*(na-1),(na+np-1)+frows*(na-1))){#(na*np),by=1)){

#for (nlinks in seq(frows*(na),frows*(na))){#(na*np),by=1)){

for (nlinks in seq((na+np-1),na*np)){  
  mincidblocked <- create_blocked_model(na,np,nlinks)
  mblocked <- sq_adjacency(mincidblocked, na, np)[[1]]
  mlaplblocked <- create_laplacian_matrix(mblocked)
  
  mincidhyper <- create_hypernested_model(na,np,nlinks)
  mhyper <- sq_adjacency(mincidhyper, na, np)[[1]]
  mlaplhyper <- create_laplacian_matrix(mhyper)
  
  dfblocked <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(mblocked)$values[1],
                                         "AdjEnergy"=AdjEnergy(eigen(mblocked)$values),
                                         "Lplspectrad"=eigen(mlaplblocked)$values[1],
                                         "LplEnergy"=LaplEnergy(eigen(mlaplblocked)$values,nlinks,na+np)))
  dfhyper <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(mhyper)$values[1],
                                       "AdjEnergy"=AdjEnergy(eigen(mhyper)$values),
                                       "Lplspectrad"=eigen(mlaplhyper)$values[1],
                                         "LplEnergy"=LaplEnergy(eigen(mlaplhyper)$values,nlinks,na+np)))

  #print(paste("links",nlinks,"Laplacian radius",eigen(mlaplhyper)$values[[1]],"sqrt links",sqrt(nlinks)))
  mincid <- create_bin_nested_model(na,np,nlinks) 
  m <- sq_adjacency(mincid,na,np)[[1]]
  mlapl <- create_laplacian_matrix(m)
    
  dfnest <- rbind(dfnest, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(m)$values[1],
                                     "AdjEnergy"=AdjEnergy(eigen(m)$values),
                                    "Lplspectrad"=eigen(mlapl)$values[1],
                                    "LplEnergy"=LaplEnergy(eigen(mlapl)$values,nlinks,na+np)))
  mincidrandom <- create_rnd_matrix(na,np,nlinks)
  mrandom <- sq_adjacency(mincidrandom,na,np)[[1]]
  mlaplrandom <- create_laplacian_matrix(mrandom)
  dfrandom <- rbind(dfrandom, data.frame("species"=na+np,"links"=nlinks,
                                         "spectrad"=eigen(mrandom)$values[1],
                                         "AdjEnergy"=AdjEnergy(eigen(mrandom)$values),
                                         "Lplspectrad"=eigen(mlaplrandom)$values[[1]],
                                         "LplEnergy"=LaplEnergy(eigen(mlaplrandom)$values,nlinks,na+np)))
  
}
dfnest$label="NESTED"
dfhyper$label="HYPERNESTED"
dfrandom$label="RND"
dfblocked$label="BLOCKED"
dfall <- rbind(dfrandom,dfnest,dfhyper)

sr <- ggplot(data=dfall,aes(x=links,y=spectrad,color=label))+geom_point(alpha=0.5)+
  geom_point(data=dfall,aes(x=links,y=sqrt(links)),color="black",size=0.5)+
  ggtitle(paste("na:",na,"np:",np))+theme_bw()
linksmed <- 1.0*na*np/2
adje <- ggplot(data=dfall,aes(x=links,y=AdjEnergy,color=label))+geom_point(alpha=0.5)+
  #geom_point(data=dfall,aes(x=links,y=linksmed-((links-linksmed)^2/linksmed)),color="black")+
  ggtitle(paste("na:",na,"np:",np))+theme_bw()

lplsr <- ggplot(data=dfall,aes(x=links,y=Lplspectrad,color=label))+geom_point(alpha=0.5)+
     ggtitle(paste("na:",na,"np:",np))+theme_bw()
lple <- ggplot(data=dfall,aes(x=links,y=LplEnergy,color=label))+geom_point(alpha=0.5)+
  ggtitle(paste("na:",na,"np:",np))+theme_bw()
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