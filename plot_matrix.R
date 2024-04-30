library(reshape2)
library(ggplot2)
library(kcorebip)
library(patchwork)
source("config_values.R")
source("manage_null_models.R")
source("colormodels.R")

lmodels <- c("NETWORK","NESTED","HNESTED","RND","SYTR","SHUFFLE","VAZ","MGEN")
weightrf <- "none"
create_dirs(weightrf)
filenames <- Sys.glob(paste0(dirnulls,"*_*_NULL_RND.csv"))
lnetw <- gsub("_NULL_RND","",gsub(".csv","",gsub(dirnulls,"",filenames)))

for (netw in lnetw) 
{
  print(netw)
  mplots <- vector(mode='list', length=length(lmodels))
  i = 1
  for (model in lmodels){
    print(model)
    if (model=="NETWORK"){
      result_analysis <- analyze_network(directory = "data/", paste0(netw,".csv"), guild_a = "A", guild_b = "P", only_NODF = TRUE)
    } else
      result_analysis <- analyze_network(directory = "smodels/none/nullmatrix/", paste0(netw,"_NULL_",model,".csv"), guild_a = "A", guild_b = "P", only_NODF = TRUE)
    
    ppi<-300
    A<-unname(result_analysis$matrix)
    A[A>0]=1
    A<-A[order(rowSums(A)),rev(order(colSums(A)))]
    
    
    longData<-melt(A)
    clist <- get_colors(list(model),colormodels)
    if (model=="NETWORK"){
      ncolor <- "blue"
    } else
      ncolor = unname(clist$lcolors[model])
      mplots[[i]]<- ggplot(longData, aes(x = Var2, y = Var1,fill=as.factor(value))) + 
                    geom_tile(color = ncolor,alpha=1) +
                    scale_fill_manual(values=c("white",ncolor))+
                    ylab("")+xlab(model)+
                    theme_void()+
                    theme(legend.position="none",
                          axis.title.x=element_text(size=14,face="bold",color="grey40",hjust=0.5),
                          plot.margin = unit(c(0.1,0.1,1,0.1), "cm"))
  
    aratio <- ncol(A)/nrow(A)
    i = i + 1
  }
  if (ncol(A)<nrow(A)) {
    pd <- ( (mplots[[1]] | mplots[[2]] | mplots[[3]])/ (mplots[[4]] | mplots[[5]] | mplots[[6]])/ (mplots[[7]] | mplots[[8]] | ggplot() + theme_void())  )
  } else
    pd <- ( (mplots[[1]] | mplots[[2]] ) /( mplots[[3]]| mplots[[4]] )/( mplots[[5]] | mplots[[6]])/ (mplots[[7]] | mplots[[8]] )  )
  pd <- pd +  plot_annotation(title = netw,theme = theme(plot.title = element_text(size = 20,face="bold",color="gray20",hjust=0.5)))
  plsize=15
  nfile <- paste0(matrixplotsdir,"/",netw,"_MATRIX",".png")
  png(nfile,width=plsize*ppi,height=plsize*ppi,res=ppi)
  print(pd)
  dev.off()
}