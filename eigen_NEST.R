library(ggplot2)
library(kcorebip)
source("manage_null_models.R")



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

create_block_model <- function(na,nb,links){
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



dfnest <- data.frame("species"=c(),"links"=c(),"spectrad"=c())
dfhyper <- dfnest
dfbad <- dfnest
na=30 # poner na mayor que np
np=30                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
frows=5


#for (nlinks in seq(frows*(na),frows*(na))){#(na*np),by=1)){

#for (nlinks in seq(frows*na,frows*na)){  

for (nlinks in seq((na+np-1),na*np)){#(na*np),by=1)){

  mincidblocks <- create_block_model(na,np,nlinks)
  
  mincidhyper <- create_hypernested_model(na,np,nlinks)
  mincidbin <- create_bin_nested_model(na,np,nlinks)
  mincidrnd <- create_rnd_matrix(na,np,nlinks)
  
  mblocks <- sq_adjacency(mincidblocks, na, np)[[1]]
  dfblocks <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,
                                       "spectrad"=eigen(mblocks)$values[1],
                                       "AdjEnergy"=AdjEnergy(eigen(mblocks)$values)))
  
  mhyper <- sq_adjacency(mincidhyper, na, np)[[1]]
  dfhyper <- rbind(dfhyper, data.frame("species"=na+np,"links"=nlinks,
                                       "spectrad"=eigen(mhyper)$values[1],
                                       "AdjEnergy"=AdjEnergy(eigen(mhyper)$values)))
  rspec <- eigen(mhyper)$values[[1]]
  print(paste("links",nlinks,"radius",rspec,"sqrt links",sqrt(nlinks)))
  
  

    mincid <- create_bin_nested_model(na,np,nlinks) 
    m <- sq_adjacency(mincid,na,np)[[1]]
    dfnest <- rbind(dfnest, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(m)$values[1],
                                       "AdjEnergy"=AdjEnergy(eigen(m)$values)))
    nfil <- na + np
    
    mincidbad <- create_rnd_matrix(na,np,nlinks)
    mbad <- sq_adjacency(mincidbad,na,np)[[1]]
    # unom <- matrix(rep(1,nfil),nrow=nfil,ncol=1)
    rspec <- eigen(mbad)$values[[1]]
    #print(paste("links",nlinks,"radius",rspec,"sqrt links",sqrt(nlinks)))
    dfbad <- rbind(dfbad, data.frame("species"=na+np,"links"=nlinks,"spectrad"=eigen(mbad)$values[1],
                                     "AdjEnergy"=AdjEnergy(eigen(mbad)$values)))
    # mcol <- m%*%unom
    #print(paste(max(mcol),min(mcol)))
    #print(norm(m, type = c("I")))
    #mpow = mtopower(m,1)
    # print(mpow)
    # print("Gershgorin disks")
    # for (k in 1:nrow(mpow)){
    #    print(sprintf("c: %d r %d",mpow[k,k],colSums(abs(mpow))[k]-mpow[k,k]))
    # }
    # print((eigen(mpow)$values))
    # print((eigen(m)$values)^0.25)
    # print("norm m4")
    # print(rowSums(m4))
    # print(max(rowSums(m4)))
    # print(max(rowSums(m4))^0.25)

   }
dfnest$label="nested"
dfhyper$label="hypernested"
dfbad$label="RND"
dfblocks$label="Block"
dfall <- rbind(dfnest,dfbad,dfblocks)
# if(sum(dfall$links==na*(1+na)/2)>0)
#   dfall[dfall$links==na*(1+na)/2,]$label="TRIANGLE"
#p <- ggplot(data=dfall,aes(x=species,y=spectrad,color=label))+geom_point()
q <- ggplot(data=dfall,aes(x=links,y=spectrad,color=label))+geom_point()+
  geom_point(data=dfall,aes(x=links,y=sqrt(links)),color="black",size=0.5)+
  theme_bw()
r <- ggplot(data=dfall,aes(x=links,y=AdjEnergy,color=label))+geom_point(alpha=0.3)+theme_bw()
  #geom_point(data=dfall,aes(x=links,y=linksmed-((links-linksmed)^2/linksmed)),color="black")