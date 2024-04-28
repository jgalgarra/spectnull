library(ggplot2)
library(kcorebip)
source("manage_null_models.R")

netw<-"M_PL_042_GC.csv"
result_analysis <- analyze_network(directory = "data/", netw, guild_a = "A", guild_b = "P", only_NODF = TRUE)
# 
# if (result_analysis$num_guild_a>result_analysis$num_guild_b){
#   m <- sq_adjacency(result_analysis$matrix,result_analysis$num_guild_a,result_analysis$num_guild_b)[[1]]
# } else
#   m <- sq_adjacency(result_analysis$matrix,result_analysis$num_guild_b,result_analysis$num_guild_a)[[1]]

m <- sq_adjacency(result_analysis$matrix,result_analysis$num_guild_a,result_analysis$num_guild_b)[[1]]
mlapl <- create_laplacian_matrix(m)
lspect <- eigen(mlapl)

algconn<-lspect$values[length(lspect$values)-1]
fiedlervect <- lspect$vectors[,length(lspect$values)-1]
clustnet <- data.frame("nodeindex"=c(),"fiedlervect"=c())
for (i in 1:length(fiedlervect))
  clustnet <- rbind(clustnet,data.frame("nodeindex"=i,"fiedlervect"=fiedlervect[i]))
for (j in 1:result_analysis$num_guild_b)
  clustnet$nodeindex[j]<-paste0("B",clustnet$nodeindex[j])
for (k in 1:result_analysis$num_guild_a)
  clustnet$nodeindex[result_analysis$num_guild_b+k]<-paste0("A",as.numeric(clustnet$nodeindex[result_analysis$num_guild_b+k])-result_analysis$num_guild_b)
clustnet <- clustnet[order(clustnet$fiedlervect),]