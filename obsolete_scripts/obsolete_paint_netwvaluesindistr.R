library(ggplot2)
library(patchwork)
library("data.table")

plot_distr_null <- function(df,nvalue,networkname="",title="",xl="",yl=""){
  plot <- ggplot(data=df)+geom_histogram(aes(x=values,fill=MODEL),bins=30,alpha=0.5)+
    scale_y_continuous(expand=c(0,0))+xlab("")+theme_bw()+ 
    theme(legend.position = "none",legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, 'cm'))+xlab(xl)+ylab(yl)
    #guides(fill=guide_legend(nrow=2, byrow=TRUE, title= "Null model"))
  plot <- plot+geom_vline(data = data.frame("val"=nvalue), aes(xintercept = val), 
                          color = "blue", size=0.5,alpha=0.5)+ggtitle(title)
  return(plot)
}


rdir <- "results/"
networkmags <- fread("results/networkmagnitudes.csv")
lnetw <- unique(networkmags$Network)

lmagnitudes <- c("adj_energy","lpl_energy","spect_rad")

binnetworks <- data.frame("network"=unique(networkmags[!networkmags$Weighted,]$Network))
nbin <- nrow(binnetworks)
binnetworks$Magnitude <- "adj_energy"
binnetworks$DU_05 <- 0
binnetworks$MGEN <- 0
binnetworks$RND <- 0
binnetworks$SHUFFLE <- 0
binnetworks$VAZ <- 0
binnetworks <- rbind(binnetworks,binnetworks,binnetworks)
for(i in ((nbin+1):(nbin*2))){
  binnetworks$Magnitude[i] <- "lpl_energy"
}
for(i in ((2*nbin+1):nrow(binnetworks))){
  binnetworks$Magnitude[i] <- "spect_rad"
}

weightednetworks <- data.frame("network"=unique(networkmags[networkmags$Weighted,]$Network))
nweight <- nrow(weightednetworks)
weightednetworks$Magnitude <- "adj_energy"
weightednetworks$SWAP <- 0
weightednetworks$SYTR <- 0
weightednetworks$RND <- 0
weightednetworks$SHUFFLE <- 0
weightednetworks$VAZ <- 0
weightednetworks <- rbind(weightednetworks,weightednetworks,weightednetworks)
for(i in ((nweight+1):(nweight*2))){
  weightednetworks$Magnitude[i] <- "lpl_energy"
}
for(i in ((2*nweight+1):nrow(weightednetworks))){
  weightednetworks$Magnitude[i] <- "spect_rad"
}

mbin <- c("SHUFFLE","VAZ","DU_05","RND","MGEN")
mweight <- c("SYTR","VAZ","SWAP","RND","SHUFFLE")
#models <- unique(dnulls$MODEL)
for (netw in lnetw){
  dnulls <-read.csv(paste0(rdir,"MODS_",netw,".csv"))
  
  if (networkmags[networkmags$Network==netw,]$Weighted){
    ntype = "WEIGHTED"
    models= mweight
  } else {
    ntype = "BINARY"
    models= mbin
  }
  z <- vector(mode='list', length=9)
  j = 1
  for (magnitude in lmagnitudes){
    for (model in models){
      magmoddata <- dnulls[dnulls$ind==magnitude & dnulls$MODEL==model,]
      nvalue <- dnulls[dnulls$ind==magnitude & dnulls$MODEL=="NETWORK",]$values
      percentile <- ecdf(magmoddata$values)
      if (ntype=="BINARY")
        binnetworks[binnetworks$network==netw & binnetworks$Magnitude==magnitude,][[model]] <- 100*percentile(nvalue)
      else
        weightednetworks[weightednetworks$network==netw & weightednetworks$Magnitude==magnitude,][[model]] <- 100*percentile(nvalue)
      if (((ntype=="BINARY") && (model %in% c("SHUFFLE","VAZ","DU_05"))) ||
          ((ntype=="WEIGHTED") && (model %in% c("SYTR","VAZ","SWAP")))){
        print(paste(netw,ntype,magnitude,model,percentile(nvalue)))
        z[[j]] <- plot_distr_null(magmoddata,nvalue,xl=model,yl=magnitude,
                                  title=sprintf("Percentile %.2f",100*percentile(nvalue)))#,shapiro.test(summary(magmoddata$values))$p.value))
        j = j+1
      }
    }
  }
  allz <- (z[[1]] | z[[2]] | z[[3]] )/ ( z[[4]] | z[[5]] | z[[6]]) /
          (z[[7]] | z[[8]] | z[[9]] )
  allz <- allz + plot_annotation(title = paste(netw,ntype))
  
  pldir <- "plots/percentiles/"
  if (!dir.exists(pldir))
    dir.create(pldir)
  ppi=300
  png(paste0(pldir,paste0("DISTS_MODEL_",netw,"_",ntype,".png")),width=20*ppi,height=12*ppi,res=ppi)
  print(allz)
  dev.off()
}
write.csv(binnetworks,paste0(rdir,"PERC_BINARY_NETWORKS.csv"),row.names=FALSE)
write.csv(weightednetworks,paste0(rdir,"PERC_WEIGHTED_NETWORKS.csv"),row.names=FALSE)