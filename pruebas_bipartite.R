# Load required packages for bipartite
# https://fukamilab.github.io/BIO202/09-B-networks.html

# library(permute)
# library(lattice)
# library(vegan)
#library("bipartite")
library(kcorebip)
library(ggplot2)
library(patchwork)
library(statnet.common)
library(network)
library(sna)

PegaRes <- function(nombre,datinp){
  return(data.frame("Objeto"=nombre,"Enlaces"=sum(datinp>0),
                    "Suma Filas"=paste(as.numeric(rowSums(datinp)), collapse = ","),
                    "Sep"=" ",
                    "Suma Cols"=paste(as.numeric(colSums(datinp)),  collapse = ",")))
}

inputdir <- "fakedata/"

datosdebipartite <- TRUE
printzigs <- FALSE

if (datosdebipartite){
  
  webs <- list(barrett1987,Safariland,  elberling1999, 
               memmott1999, motten1982, olesen2002aigrettes)
  webs.names <- c( "barrett1987","Safariland", "elberling1999", 
                  "memmott1999", "motten1982", "olesen2002aigrettes") 
  
  webs <- list(motten1982)
  webs.names <- c("motten1982") 
  
  # Re-name the datasets according to the sites for each plant-pollinator network
  names(webs) <- webs.names
  
  netw <- webs.names[[1]]
    write.csv(as.matrix(webs[[1]]),paste0(inputdir,netw,".csv"))
} else {
  netw <- "M_PL_024"
}



printtofile <- TRUE
odir <- "grafresults_null/"
result_analysis <- analyze_network(directory = inputdir, paste0(netw,".csv"), guild_a = "Plant", guild_b = "Pollinator", plot_graphs = FALSE)

if (!datosdebipartite)
  webs <- list(result_analysis$matrix)
if (printzigs)
  pgr <- ziggurat_graph(inputdir,paste0(netw,".csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE,weighted_links = "ln")

redpesada <- sum(webs[[1]]>1)
#redpesada <- TRUE
print(paste(netw,"A",result_analysis$num_guild_a," B",result_analysis$num_guild_b,
            " Links",result_analysis$links, " Pesada",redpesada))

#plotweb(webs$Argentina, text.rot=90, col.low = "green", col.high = "blue")
# Calculate network metric nestedness for all plant-pollinator sites
net.metrics.nest <- lapply(webs, networklevel, index = 'nestedness') 

# Calculate network metric links per species for all plant-pollinator sites
net.metrics.links <- lapply(webs, networklevel, index = 'links per species') 

# Time consuming step!

# Load environment (already saved objects)
#load("data/network_analysis_example.RData")

Nredes = 1

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

print(paste(netw,"shuffle A",ncol(net.nulls.shuffle[[1]][[1]])," B",nrow(net.nulls.shuffle[[1]][[1]]),
            " Links",sum(net.nulls.shuffle[[1]][[1]]>0),"Pesada",sum(net.nulls.shuffle[[1]][[1]]>1)))

# # Save null objects
# save(net.nulls.r2d, net.nulls.vaz, net.nulls.swap, file = "data/network_analysis_example.RData")

comparativa = data.frame("Objeto"=c(),"Enlaces"=c(),"Suma Filas"=c(),"Sep"=c(),"Suma Cols"=c())

comparativa=rbind(comparativa,data.frame("Objeto"=netw,"Enlaces"=sum(webs[[1]]>0),
                                        "Suma Filas"=paste(as.numeric(webs[[1]]), collapse = ","),
                                        "Sep"=" ",
                                        "Suma Cols"=paste(as.numeric(colSums(webs[[1]])),  collapse = ",")))
comparativa = rbind(comparativa,PegaRes("r2d",net.nulls.r2d[[1]][[1]]))
comparativa = rbind(comparativa,PegaRes("vaz",net.nulls.vaz[[1]][[1]]))
comparativa = rbind(comparativa,PegaRes("swap",net.nulls.swap[[1]][[1]]))
# visweb(webs$Argentina)
# visweb(net.nulls.r2d[[1]][[1]])
# visweb(net.nulls.vaz[[1]][[1]])
# visweb(net.nulls.swap[[1]][[1]])

if (printzigs){
  write.csv(net.nulls.r2d[[1]][[1]],paste0(inputdir,netw,"_r2d.csv"))
  pgr <- ziggurat_graph(inputdir,paste0(netw,"_r2d.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE,weighted_links = "ln")
  write.csv(net.nulls.vaz[[1]][[1]],paste0(inputdir,netw,"_vaz.csv"))
  pgr <- ziggurat_graph(inputdir,paste0(netw,"_vaz.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE,weighted_links = "ln")
 }


print("net.mgen")
net.nulls.mgen <- lapply(webs, nullmodel, method = "mgen", N = Nredes)
write.csv(net.nulls.mgen[[1]][[1]],paste0(inputdir,netw,"_mgen.csv"))
print("mgen")
net.nulls.mgen <- lapply(webs, nullmodel, method = "mgen", N = Nredes)
print(paste(netw,"shuffle A",ncol(net.nulls.shuffle[[1]][[1]])," B",nrow(net.nulls.shuffle[[1]][[1]]),
            " Links",sum(net.nulls.shuffle[[1]][[1]]>0),"Pesada",sum(net.nulls.shuffle[[1]][[1]]>1)))
if (printzigs)
  pgr <- ziggurat_graph(inputdir,paste0(netw,"_mgen.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE)
comparativa = rbind(comparativa,PegaRes("mgen",net.nulls.mgen[[1]][[1]]))

print("shuffle.web")
net.nulls.shuffle <- lapply(webs, nullmodel, method = "shuffle.web", N = Nredes)
print(paste(netw,"shuffle A",ncol(net.nulls.shuffle[[1]][[1]])," B",nrow(net.nulls.shuffle[[1]][[1]]),
            " Links",sum(net.nulls.shuffle[[1]][[1]]>0),"Pesada",sum(net.nulls.shuffle[[1]][[1]]>1)))

write.csv(net.nulls.shuffle[[1]][[1]],paste0(inputdir,netw,"_shuffle.csv"))
if (printzigs)
  pgr <- ziggurat_graph(inputdir,paste0(netw,"_shuffle.csv"),plotsdir=odir,print_to_file = printtofile,show_title = TRUE)
comparativa = rbind(comparativa,PegaRes("shuffle",net.nulls.shuffle[[1]][[1]]))
write.csv(comparativa,paste0(inputdir,"comp_",netw,".csv"))