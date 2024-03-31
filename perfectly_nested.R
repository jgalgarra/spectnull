#source("manage_null_models.R")

create_myperfect_nested_model <- function(na,nb,nlinks){
  if (nlinks > na*nb)
    stop("FATAL ERROR. Number of links greater than product of number of species")
  if (na > nb ){
    nr <- nb
    nc <- na
  } else {
    nr <- na
    nc <- nb
  }
  x <- seq(0,1,by=1/(nc))
  linkstr <- nc*(1+nr)/2
  ftr <- nlinks/linkstr
  y <- 1-(1-x^(1/ftr))^ftr
  y <- rev(y)
  pnr <- ceiling((y*nr))
  pnr[pnr==0]=1
  # n = 0
  # if (nlinks < (nc*nr/2))
  #   while (n<2){
  #     if (sum(pnr)!=nlinks){
  #       pnr[2:(length(pnr)-1)] <- min(nr,pnr[2:(length(pnr)-1)]*nlinks/sum(pnr))
  #       pnr <- round(pnr)
  #     }
  #     n<-n+1
  #   }
  mincid <- matrix(rep(0,nc*nr),nrow=nr,ncol=nc)
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
  if (ncol(mincid)>nrow(mincid))
    mincid <- t(mincid)
  return(mincid)
}



na = 50
np = 50
linkstr <- na*(1+np)/2

nlinks <- 170

nestmatrix <- create_myperfect_nested_model(na,np,0.9*na*np)
print(nested(nestmatrix,c("binmatnest","NODF","wine")))

# binnestmatrix <- create_bin_nested_model(na,np,nlinks)
# print(nested(binnestmatrix,c("binmatnest","NODF","wine")))