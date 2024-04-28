library(reshape2)
library(ggplot2)
library(kcorebip)
library(patchwork)
source("manage_null_models.R")
source("colormodels.R")

netw<-"M_PL_010"
lmodels <- c("EMPIRICAL","NESTED","HNESTED","RND","SYTR","SHUFFLE","VAZ","MGEN")
mplots <- vector(mode='list', length=length(lmodels))
i = 1
for (model in lmodels){
  if (model=="EMPIRICAL"){
    result_analysis <- analyze_network(directory = "data/", paste0(netw,".csv"), guild_a = "A", guild_b = "P", only_NODF = TRUE)
  } else
    result_analysis <- analyze_network(directory = "smodels/none/nullmatrix/", paste0(netw,"_NULL_",model,".csv"), guild_a = "A", guild_b = "P", only_NODF = TRUE)
  
  ppi<-300
  A<-unname(result_analysis$matrix)
  A[A>0]=1
  A<-A[order(rowSums(A)),rev(order(colSums(A)))]
  
  
  longData<-melt(A)
  clist <- get_colors(list(model),colormodels)
  if (model=="EMPIRICAL"){
    ncolor <- "black"
  } else
    ncolor = unname(clist$lcolors[model])
    mplots[[i]]<- ggplot(longData, aes(x = Var2, y = Var1,fill=as.factor(value))) + 
                  geom_tile(color = ncolor,alpha=1) +
                  scale_fill_manual(values=c("white",ncolor))+
                  ggtitle(model)+
                  theme_void()+theme(legend.position="none",plot.title = element_text(size=20,hjust = 0.5))

  aratio <- ncol(A)/nrow(A)
  i = i + 1
}
if (ncol(A)<nrow(A)) {
  pd <- ( (mplots[[1]] | mplots[[2]] | mplots[[3]])/ (mplots[[4]] | mplots[[5]] | mplots[[6]])/ (mplots[[7]] | mplots[[8]] | ggplot() + theme_void())  )
} else
  pd <- ( (mplots[[1]] | mplots[[2]] ) /( mplots[[3]]| mplots[[4]] )/( mplots[[5]] | mplots[[6]])/ (mplots[[7]] | mplots[[8]] )  )
plsize=15
nfile <- paste0("matrix/",netw,"_MATRIX",".png")
png(nfile,width=plsize*ppi,height=plsize*ppi,res=ppi)
print(pd)
dev.off()