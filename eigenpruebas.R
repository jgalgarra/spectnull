library(ggplot2)

create_nested_matrix <- function(na,np,triangle=FALSE){
nfil = na+np
ncol = nfil
nunos = 0
m<-matrix(rep(0,nfil*ncol), nrow = nfil, ncol = ncol)
if(triangle){
for (i in 1:na)
  for(j in 1:(np-i+1)){
    if (j>0)
      m[i,na+j]=1
  }
}
else{
  m[1,(np+1):(np+na)] <- 1
  m[2:np,np+1]<-1
}
m[(na+1):(np+na),1:na] <-  t(m[1:na,(na+1):(np+na)])
numberones=sum(m)/2
unom <- matrix(rep(1,nfil),nrow=nfil,ncol=1)
return(list("adjmatrix"=m,"numberones"=numberones))
}

mtopower <- function(m,pow){
  mres <- m
  if (pow>1){
    for (i in 2:pow){
      mres <- mres %*% t(m)
    }
  }
  return(mres)
}
create_rnd_matrix <- function(na,np,numberones){
mbad <- matrix(rep(0,(na+np)^2),nrow=na+np,ncol=na+np)
rndlinks <- matrix(rep(0,na*np),nrow=na,ncol=np)
rndlinks[sample(seq(1,na*np),numberones)] <- 1
randomm <- matrix(rndlinks,nrow=na,ncol=np)
mbad[1:na,(na+1):(na+np)]<-randomm
mbad[(na+1):(na+np),1:na]<-t(randomm)
return(mbad)
}

dfnest <- data.frame("species"=c(),links=c(),"radius"=c())
dfbad <- dfnest

for (na in 5:200){
#np=3
#na=3
np=na
nestedm <- create_nested_matrix(na,np,triangle=TRUE)
m <- nestedm$adjmatrix
numberones <- nestedm$numberones
  
dfnest <- rbind(dfnest, data.frame("species"=na+np,"links"=numberones,"radius"=eigen(m)$values[1]))

nfil <- na + np
mbad <- create_rnd_matrix(na,np,numberones)
unom <- matrix(rep(1,nfil),nrow=nfil,ncol=1)

dfbad <- rbind(dfbad, data.frame("species"=na+np,"links"=numberones,"radius"=eigen(mbad)$values[1]))



# for (i in (1:nfil))
#   for (j in (1:ncol))
#     print(paste(i,j))
#     if(i<j)
#       m[i,j]=m[j,i]
#print(m)
#print(eigen(m)$values)
mcol <- m%*%unom
print(paste(max(mcol),min(mcol)))
print(norm(m, type = c("I")))
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
dfnest$label="nested"
dfbad$label="rnd"
dfall <- rbind(dfnest,dfbad)
p <- ggplot(data=dfall,aes(x=species,y=radius,color=label))+geom_point()
