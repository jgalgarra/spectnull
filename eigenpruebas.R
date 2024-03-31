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
if (na>np)
  randomm=t(randomm)
return(randomm)
}


dfnest <- data.frame("species"=c(),links=c(),"rowfillius"=c())
dfbad <- dfnest

for (na in 8:8){
#np=3
#na=3
np=na+6
#nestedm <- create_nested_matrix(na,np,triangle=TRUE)
# m <- nestedm$adjmatrix
# nlinks <- nestedm$nlinks

#for (links in seq(na+np-1,floor(na*(np+1)/2))) {# seq(na+np-1,na*np)){

    
  #nlinks <- na*(1+np)/2



result_analysis <- analyze_network(directory = "data/", "M_PL_004.csv", only_NODF = TRUE)
unmatrix <- unname(result_analysis$matrix)
rs <- rowSums(unmatrix)
positions <- c()
for (i in 1:length(rs)){
  pm <- which(rs==max(rs))[1]
  positions <- c(positions,pm)
  rs[pm] <- 0
}
unmatrix <- unmatrix[positions,]
rs <- colSums(unmatrix)
positions <- c()
for (i in 1:length(rs)){
  pm <- which(rs==max(rs))[1]
  positions <- c(positions,pm)
  rs[pm] <- 0
}
unmatrix <- unmatrix[,positions]
nlinks <- result_analysis$links
na <- result_analysis$num_guild_a
np <- result_analysis$num_guild_b
mor <- sq_adjacency(unmatrix, na, np)[[1]]

  mincidhyper <- create_hypernested_model(na,np,nlinks)
  
  mhyper <- sq_adjacency(mincidhyper, na, np)[[1]]
  #nlinks <- links
  print(paste("links",nlinks))
  mincid <- create_perfect_nested_model(na,np,nlinks) 
  m <- sq_adjacency(mincid,na,np)[[1]]
  # mincid[5,5] <- 1
  # mincid[4,7] <- 0

  # nstmtx <- sq_adjacency(mincid, na, np)
  # m <- nstmtx[[1]]  
  

  dfnest <- rbind(dfnest, data.frame("species"=na+np,"links"=nlinks,"rowfillius"=eigen(m)$values[1]))
  nfil <- na + np
  
  mincidbad <- create_rnd_matrix(na,np,nlinks)
  mbad <- sq_adjacency(mincidbad,na,np)[[1]]
  unom <- matrix(rep(1,nfil),nrow=nfil,ncol=1)
  dfbad <- rbind(dfbad, data.frame("species"=na+np,"links"=nlinks,"rowfillius"=eigen(mbad)$values[1]))
  mcol <- m%*%unom
  #print(paste(max(mcol),min(mcol)))
  #print(norm(m, type = c("I")))
  mpow = mtopower(m,1)
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
#}
dfnest$label="nested"
dfbad$label="rnd"
dfall <- rbind(dfnest,dfbad)
if(sum(dfall$links==na*(1+na)/2)>0)
  dfall[dfall$links==na*(1+na)/2,]$label="TRIANGLE"
p <- ggplot(data=dfall,aes(x=species,y=rowfillius,color=label))+geom_point()
q <- ggplot(data=dfall,aes(x=links,y=rowfillius,color=label))+geom_point()
