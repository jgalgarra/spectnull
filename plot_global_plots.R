library(ggplot2)
library(reshape2)
library(forcats)

plot_histogram <- function(ndata){
  ndata <- ndata[ndata$Weighted,]
  #mods <- c("NETWORK","SYTR","RND","VAZ")
  mods <- unique(ndata$Model)
  magn <- "lpl_energy"
  ndata <- ndata[ndata$Magnitude==magn & is.element(ndata$Model,mods),]
  networks <- unique(ndata$Network)
  ndata$rel <- 0
  for (i in 1:nrow(ndata)){
    valnet <- ndata[ndata$Network==ndata$Network[i] & ndata$Model=="NETWORK",]$value
    ndata$rel[i] <- (ndata$value[i]-ndata[ndata$Network==ndata$Network[i],]$value)
  }
  norder <- data.frame("Network"=NULL,"rel"=NULL)
  for (n in networks){
    norder <- rbind(norder,data.frame("Network"=n,"rel"=max(abs(ndata[ndata$Network==n,]$rel))))
  }
  #norder <- ndata[(ndata$Model=="NETWORK") & (ndata$Magnitude==magn),]
  ndata$Network<- fct_relevel(ndata$Network, norder[order(norder$rel),]$Network)
  plot <- ggplot(data=ndata[ndata$Model!="NETWORK",],aes(x=Network,y=rel,fill=Model))+
    scale_y_sqrt()+
    #geom_density(aes(x=value,fill=Model),color="transparent",alpha=0.3)+
    geom_point(size=3,shape=21,color="transparent",alpha=0.5)+xlab("")+theme_bw()+ 
    theme(legend.position = "bottom",legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, 'cm'))+coord_flip()
    guides(fill=guide_legend(nrow=2, byrow=TRUE, title= "Null model"))
  # plot <- plot+geom_vline(data = data.frame("val"=nvalue), aes(xintercept = val), 
  #                         color = "blue", size=0.5,alpha=0.5)+ggtitle(title)
  return(plot)
  
}

plot_magnitude <- function(ndata,magn="lpl_energy",weight = TRUE,scale="none")
{

  ndata <- ndata[ndata$Weighted == weight,]
  ndata <- ndata[is.element(ndata$Model,c("NETWORK","VAZ","SYTR","RND")) & ndata$Magnitude==magn,]
  if (scale=="none")
    ndata$pval <- ndata$value
  if (scale=="nodes")
    ndata$pval <- ndata$value/(ndata$NodesA+ndata$NodesB)
  if (scale=="links")
    ndata$pval <- ndata$value/ndata$Links
  #ndata <- ndata[order(ndata[ndata$Model=="NETWORK",]$pval),]
  
  ndata$Network<- fct_reorder(ndata$Network,ndata$pval)
  plot <- ggplot(data=ndata)+
    geom_point(aes(y=pval,x=Network,fill=Model),shape=21,alpha=0.5,size=3,color="transparent")+theme_bw()+    
    theme(legend.position = "bottom",legend.text = element_text(size = 8),
          legend.key.size = unit(0.3, 'cm'))+coord_flip()
  return(plot)
}

rdir <- "results/"
raw_data <- read.csv(paste0(rdir,"NetworkMagnitudes.csv"),stringsAsFactors = TRUE)
models_data <-melt(raw_data, id.var = names(raw_data)[1:which(names(raw_data)=="Model")], variable.name = 'Magnitude')
p <- plot_histogram(models_data)
q <- plot_magnitude(models_data,weight = FALSE, magn="lpl_energy", scale="nodes")