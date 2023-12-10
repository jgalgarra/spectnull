library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)
library(ggrepel)
set.seed(122)
normresults <- "results/normalized/"
results <- "results/"
networkmags <- read.csv(paste0(results,"networkmagnitudes.csv"))
normdistavgs <- read.csv(paste0(normresults,"meannormdistances.csv"))
alldistances <- read.csv(paste0(normresults,"ALLNORMALIZED.csv"))
networktypes <- read.csv(paste0(normresults,"networktypes.csv"))
alldistances$network <- gsub("MINMAX_MODS_","",alldistances$network)

# Global boxplots

z <- vector(mode='list', length=6)
j = 1
for (n in c("WEIGHTED","BINARY"))
{
  for (s in c("adj_norm_energy","lpl_norm_energy","spect_rad_norm")){
    z[[j]] <- ggplot(data=normdistavgs[normdistavgs$networktype==n & normdistavgs$disttype==s,],aes(MODEL,meannormdist)) +
      geom_boxplot(fill="blue", alpha = 0.5)+geom_hline(yintercept=0,  
                                                        color = "pink",linewidth=1, alpha=0.3)+
      ggtitle(paste(n,s))+theme_bw()+ theme(
        legend.key = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(color="ivory3",linewidth=0.1),
        panel.grid.major.x = element_line(color="ivory3",linewidth=0.1))
    j <- j + 1
  }
}
allz <- (z[[1]] | z[[2]] | z[[3]]) / (z[[4]] | z[[5]] | z[[6]])

pldir <- "plots/analysis/"
if (!dir.exists(pldir))
  dir.create(pldir)
ppi=300
png(paste0(pldir,"ALLZ.png"),width=12*ppi,height=6*ppi,res=ppi)
print(allz)
dev.off()

# Linear plot
disttype <- c("lpl_norm_energy")
plotdata <- alldistances[alldistances$ind %in% disttype,]
plotdata$ntype <- "BINARY"
for (i in 1:nrow(networktypes)){
  plotdata[plotdata$network == networktypes$network[i],]$ntype <- networktypes$type[i]
}

plotdata %>%
       group_by(network,ind,MODEL)%>%
       summarise(avg = mean(normdist))

lm = c("VAZ","RND")
networktype = "WEIGHTED"
disttype <- c("lpl_norm_energy","adj_norm_energy")
plotdata <- alldistances[alldistances$ind %in% disttype,]
plotdata$ntype <- ""
for (i in 1:nrow(networktypes)){
  plotdata[plotdata$network == networktypes$network[i],]$ntype <- networktypes$type[i]
}
plotdata <- plotdata[plotdata$ntype == networktype,]
plotdata <- plotdata %>%
  group_by(network,ind,MODEL)%>%
  summarise(avg = mean(normdist))
dx <- plotdata[plotdata$ind == "adj_norm_energy",]
dy <- plotdata[plotdata$ind == "lpl_norm_energy",]
pld <- data.frame("network"=dx$network,"model"=dx$MODEL,"adj_norm_energy"=dx$avg,"lpl_norm_energy"=dy$avg)

scatterpl <- ggplot(data=pld[pld$model %in% lm,]) +
  geom_point(aes( adj_norm_energy ,lpl_norm_energy,fill=as.factor(model)),
               color="transparent",shape=21, size=2,  alpha = 0.3)+
  geom_text_repel( aes(x=adj_norm_energy, y=lpl_norm_energy,
                 label=gsub("SD0","SD",gsub("M0","M",gsub("HP0","HP",gsub("PL0","PL",gsub("M","",gsub("RA","",gsub("00","0",gsub("_","",network)))))))),
                 color=as.factor(model)),max.overlaps=Inf,size=4)+
  geom_hline(yintercept=0,  
             color = "blue",linewidth=0.5, alpha=0.3)+
  ggtitle(networktype)+ coord_fixed()+
  theme_bw()#+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1))

png(paste0(pldir,"SCATTER",networktype,"_",paste(lm,collapse="-"),"_",paste(disttype,collapse="-"),".png"),width=12*ppi,height=12*ppi,res=ppi)
print(scatterpl)
dev.off()

# Model comparison plots
lm = c("SYTR","RND")  # The first model will set the order
networktype = "WEIGHTED"
disttype <- c("lpl_norm_energy")
#plotdata <- plotdata[grepl("HP",plotdata$network) ,]
plotdata <- alldistances[!is.na(alldistances$network),]
plotdata <- plotdata[plotdata$ind %in% disttype & plotdata$MODEL %in% lm,]

plotdata$ntype <- "BINARY"
for (i in 1:nrow(networktypes)){
  plotdata[plotdata$network == networktypes$network[i],]$ntype <- networktypes$type[i]
}

plotdata <- plotdata[plotdata$ntype==networktype,]

mdata <- normdistavgs[normdistavgs$MODEL == lm[1] & normdistavgs$disttype %in% disttype & normdistavgs$networktype==networktype,]
plotdata$network <- factor(plotdata$network,levels=mdata[order(mdata$meannormdist),]$network)
r <- ggplot(data=plotdata,aes(network,normdist,fill=MODEL)) +
  geom_boxplot(alpha = 0.5)+theme_bw()+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1))

q <- ggplot(data=plotdata) +
  geom_jitter(aes( network ,normdist,fill=MODEL),position=position_jitter(0.2),
              color="transparent",shape=21, size=0.5,  alpha = 0.1)+
  ggtitle(networktype)+ guides(fill = guide_legend(override.aes = list(size=8,alpha=0.5)))+
  theme_bw()+theme(axis.text.x = element_text( size=8, angle = 85, hjust = 1))

png(paste0(pldir,"NORMDIST_",networktype,"_",paste(disttype,collapse="-"),".png"),width=10*ppi,height=6*ppi,res=ppi)
print(q)
dev.off()
# 
# q <- ggplot(data=weightdata) +
#   geom_point(aes(network,meannormdist,fill=MODEL),shape=21,  alpha = 0.5)+coord_flip()

