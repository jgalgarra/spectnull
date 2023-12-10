library(kcorebip)
library(ggplot2)
library(patchwork)
library(statnet.common)
library(network)
library(sna)
source("../synthrade_nullmodel.R")

PegaRes <- function(nombre,datinp){
  return(data.frame("Objeto"=nombre,"Enlaces"=sum(datinp>0),
                    "Suma Filas"=paste(as.numeric(rowSums(datinp)), collapse = ","),
                    "Sep"=" ",
                    "Suma Cols"=paste(as.numeric(colSums(datinp)),  collapse = ",")))
}

inputdir <- "../data/"
outputdir <- "matrices/"
dir.create(outputdir, showWarnings = FALSE)

netw <- "modelitobin"
printtofile <- TRUE
odir <- "grafresults_null/"
result_analysis <- analyze_network(directory = inputdir, paste0(netw,".csv"), guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)

webs <- list(result_analysis$matrix)
redpesada <- any(webs[[1]]>1)
print(paste(netw,"A",result_analysis$num_guild_a," B",result_analysis$num_guild_b,
            " Links",result_analysis$links, " Pesada",redpesada))

# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'nestedness') 

# Calculate network metric links per species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species') 


Nredes = 1
comparativa = data.frame("Objeto"=c(),"Enlaces"=c(),"Suma Filas"=c(),"Sep"=c(),"Suma Cols"=c())
comparativa=rbind(comparativa,data.frame("Objeto"=netw,"Enlaces"=sum(webs[[1]]>0),
                                         "Suma Filas"=paste(as.numeric(rowSums(webs[[1]])), collapse = ","),
                                         "Sep"=" ",
                                         "Suma Cols"=paste(as.numeric(colSums(webs[[1]])),  collapse = ",")))
if (redpesada){
print("net.nulls.r2d")
# Make null models for all sites using the r2dtable null
net.nulls.r2d <- lapply(webs, nullmodel, method = "r2dtable", N = Nredes) 

# Make null models for all sites using the vaznull null
print("net.nulls.vaz")
net.nulls.vaz <- lapply(webs, nullmodel, method = "vaznull", N = Nredes) 

print("net.nulls.swap")
net.nulls.swap <- lapply(webs, nullmodel, method = "swap.web", N = Nredes) 

print("net.nulls.shuffle")
net.nulls.shuffle <- lapply(webs, nullmodel, method = "shuffle.web", N = Nredes) 

print(paste(netw,"r2d A",ncol(net.nulls.r2d[[1]][[1]])," B",nrow(net.nulls.r2d[[1]][[1]]),
            " Links",sum(net.nulls.r2d[[1]][[1]]>0),"Pesada",sum(net.nulls.r2d[[1]][[1]]>1)))

print(paste(netw,"vaz A",ncol(net.nulls.vaz[[1]][[1]])," B",nrow(net.nulls.vaz[[1]][[1]]),
            " Links",sum(net.nulls.vaz[[1]][[1]]>0),"Pesada",sum(net.nulls.vaz[[1]][[1]]>1)))

print(paste(netw,"swap A",ncol(net.nulls.swap[[1]][[1]])," B",nrow(net.nulls.swap[[1]][[1]]),
            " Links",sum(net.nulls.swap[[1]][[1]]>0),"Pesada",sum(net.nulls.swap[[1]][[1]]>1)))

comparativa = rbind(comparativa,PegaRes("r2d",net.nulls.r2d[[1]][[1]]))
comparativa = rbind(comparativa,PegaRes("vaz",net.nulls.vaz[[1]][[1]]))
comparativa = rbind(comparativa,PegaRes("swap",net.nulls.swap[[1]][[1]]))

write.csv(net.nulls.r2d[[1]][[1]],paste0(outputdir,netw,"_r2d.csv"))
#pgr <- ziggurat_graph(inputdir,paste0(netw,"_r2d.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE,weighted_links = "ln")
write.csv(net.nulls.vaz[[1]][[1]],paste0(outputdir,netw,"_vaz.csv"))
#pgr <- ziggurat_graph(inputdir,paste0(netw,"_vaz.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE,weighted_links = "ln")


incidmatrix <- SynthTradeNull(result_analysis)$matrix_synth
write.csv(incidmatrix,paste0(outputdir,netw,"_sytr.csv"))
comparativa = rbind(comparativa,PegaRes("sytr",incidmatrix))
}

net.nulls.shuffle <- lapply(webs, nullmodel, method = "shuffle.web", N = Nredes)
print(paste(netw,"shuffle A",ncol(net.nulls.shuffle[[1]][[1]])," B",nrow(net.nulls.shuffle[[1]][[1]]),
            " Links",sum(net.nulls.shuffle[[1]][[1]]>0),"Pesada",sum(net.nulls.shuffle[[1]][[1]]>1)))
write.csv(net.nulls.shuffle[[1]][[1]],paste0(outputdir,netw,"_shuffle.csv"))
#pgr <- ziggurat_graph(inputdir,paste0(netw,"_shuffle.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE)
comparativa = rbind(comparativa,PegaRes("shuffle",net.nulls.shuffle[[1]][[1]]))


if (!redpesada){
print("net.mgen")
net.nulls.mgen <- lapply(webs, nullmodel, method = "mgen", N = Nredes)
write.csv(net.nulls.mgen[[1]][[1]],paste0(outputdir,netw,"_mgen.csv"))

print(paste(netw,"shuffle A",ncol(net.nulls.shuffle[[1]][[1]])," B",nrow(net.nulls.shuffle[[1]][[1]]),
            " Links",sum(net.nulls.shuffle[[1]][[1]]>0),"Pesada",sum(net.nulls.shuffle[[1]][[1]]>1)))
#pgr <- ziggurat_graph(inputdir,paste0(netw,"_mgen.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE)
comparativa = rbind(comparativa,PegaRes("mgen",net.nulls.mgen[[1]][[1]]))
}


write.csv(comparativa,paste0(outputdir,"comp_",netw,".csv"),row.names = FALSE)
