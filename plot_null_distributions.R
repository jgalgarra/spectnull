# Functions called from main script null_models_spectral_analysis


plot_distr_null <- function(df,nvalue,nlinks=0,networkname="",title="",nestedvalue="",hypernestedvalue=""){
  lmodels <- unique(df$type)
  clist <- get_colors(lmodels,colormodels)
  plot <- ggplot(data=df)+geom_histogram(aes(x=vals,fill=type),bins=30,position="identity",alpha=0.6)+
     xlab("")+ scale_y_continuous(expand=c(0,0))+
     scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
     theme_bw()+
     theme(legend.position = "bottom",legend.text = element_text(size = 8),
           panel.grid = element_blank(),
           legend.key.size = unit(0.3, 'cm'))+
     guides(fill=guide_legend(nrow=2, byrow=TRUE))
  plot <- plot+geom_vline(xintercept = nvalue,color = "blue", size=0.75,alpha=0.33)+
          ggtitle(title)
  if (nestedvalue!="")
    plot <- plot+geom_vline(xintercept = nestedvalue,color = "magenta", size=0.75,alpha=0.33)
  if (hypernestedvalue!="")
    plot <- plot+geom_vline(xintercept = hypernestedvalue,color = "red4", size=0.75,alpha=0.5,linetype="dotted")
  # if (nlinks>0)
  #   plot <- plot+geom_vline(xintercept = sqrt(nlinks),color = "grey7", size=0.75,alpha=0.5,linetype="dotted")
  
  return(plot)
}

plot_spectr <- function(df,title=""){
  df$type <- fct_relevel(df$type, "RND", "VAZ", "NETWORK")
  if (nrow(df)<20)
    psize <- 3
  else
    psize <- 1.5
  lmodels <- unique(df$type)
  clist <- get_colors(lmodels,colormodels)
  plot <- ggplot(data=df)+geom_point(aes(x=ind,y=val,color=type),size=psize,
                                     position="identity",alpha=0.5)+
    xlab("")+ylab("Eigen values")+theme_bw()+ggtitle(title)+
    scale_y_continuous(expand=c(0,0))+
    scale_fill_manual(name='Null model ', breaks=clist$lnames, values=clist$lcolors)+
    theme_bw()+ theme(legend.position = "bottom",legend.text = element_text(size = 8),
                      legend.key.size = unit(0.3, 'cm'),
                      axis.title.x = element_text(size = 1))+
    guides(color=guide_legend(title= "Null model"))
  return(plot)
}

plot_bipartite <- function(bg, aspect_ratio = 9/35, vframecolor = "grey70", vlabelcex = 4,
                           vsize = 4, vcolor = c("lightblue","pink2"), labelcolor = c("blue","red"),
                           framedisp = FALSE,  color_link = "grey50", vertical = FALSE)
{
  l <- layout.bipartite(bg)
  if (vertical)
    la <- l[, c(2,1)]
  else
    la <- l
  plot.igraph(bg, layout= la,asp=aspect_ratio,vertex.frame.color=vframecolor,
              vertex.label.cex=vlabelcex,vertex.label.color="black",
              vertex.size=vsize, edge.color= color_link, frame = framedisp,
              vertex.color=vcolor[V(bg)$type+1])
}

# Plots the probability distributions of all models for a particular magnitude and network
plot_distributions <- function(mresults,magnitude,stitle,vtitle,nname,fdummy,nlinks=0,nestedv="",hypernestedv=""){
 print(magnitude)
  mdls <- names(mresults)
  col <- which(names(mresults[[1]])==magnitude)
  df_nulls <- data.frame("vals"=mresults[[1]][,col],"type"=mdls[[1]])
  for (i in 2:length(mdls))
    df_nulls <- rbind(df_nulls,data.frame("vals"= mresults[[mdls[i]]][,col],"type"=mdls[i]))
  pimage <- plot_distr_null(df_nulls,vtitle,networkname=nname,nlinks=nlinks,
                            title=sprintf("%s %.2f",stitle,vtitle),nestedvalue = nestedv,
                            hypernestedvalue=hypernestedv)
  calc_values <- list("df_nulls"= df_nulls, "pimage" = pimage)
  return(calc_values)
}

plot_all_distr <- function(nname,plotzigs,nnm){
  mspec <- nnm$mnames
  df_adj_model <- lapply(1:length(mspec), function(x) c())
  names(df_adj_model) <- mspec
  df_lpl_model <- df_adj_model
  df_lpl_weighted_model <- df_adj_model
  df_adj_model[["NETWORK"]] <- store_spect_values(nnm$bmag$adj_spect$values,"NETWORK")
  df_lpl_model[["NETWORK"]] <- store_spect_values(nnm$bmag$lpl_spect$values,"NETWORK")
  if (nnm$weighted_network)
    df_lpl_weighted_model[["NETWORK"]] <- store_spect_values(nnm$wmag$lpl_weighted_spect$values,"NETWORK")
  for (m in mspec[2:length(mspec)])
  {
    df_adj_model[[m]] <- store_spect_values(nnm$specresults[[m]]$null_spect$values,m)
    df_lpl_model[[m]] <- store_spect_values(nnm$specresults[[m]]$lpl_spect_null$values,m)
    if (nnm$weighted_network)
      df_lpl_weighted_model[[m]] <- store_spect_values(nnm$specresults[[m]]$lpl_weighted_spect_null$values,m)
  }
  pl_adj_spectrum <- plot_spectr(do.call(rbind,df_adj_model),title="Adjacency spectrum")
  pl_lpl_spectrum <- plot_spectr(do.call(rbind,df_lpl_model),title="Laplacian spectrum")
  if (nnm$weighted_network)
    pl_lpl_weighted_spectrum <- plot_spectr(do.call(rbind,df_lpl_weighted_model),title="Weighted Laplacian spectrum")
  nnm$modresultsweighted <- nnm$modresults
  for (mname in c("RND","BSHUFFLE","BVAZ","SYTR"))
    nnm$modresultsweighted <- remove_model_data(nnm$modresultsweighted,mname)
  
  if (nnm$weighted_network)
    for (mname in c("SHUFFLE","VAZ","SWAP"))
      nnm$modresults <- remove_model_data(nnm$modresults,mname)
  else
    for (mname in c("BSHUFFLE","BVAZ"))
      nnm$modresults <- remove_model_data(nnm$modresults,mname)
  for (mname in c("WRND","WSYTR","MGEN","PATEFIELD"))
    nnm$modresults <- remove_model_data(nnm$modresults,mname)
  
  if (nnm$weighted_network){
    if (weightrf=="none")
      trsuff <- ""
    else
      trsuff <- paste0("[",weightrf,"]")
    pd_weighted_adj_energy <- plot_distributions(nnm$modresultsweighted,"adj_weighted_energy",paste("Weighted Adj. Energy",trsuff),nnm$wmag$adj_weighted_energy,netw,fract_dummy,nestedv=nnm$pnw$weightednstadj_energy)
    pd_weighted_lpl_energy <- plot_distributions(nnm$modresultsweighted,"lpl_weighted_energy",paste("Weighted Laplacian energy",trsuff),nnm$wmag$lpl_weighted_energy,netw,fract_dummy,nestedv=nnm$pnw$weightednstlpl_energy)
    pd_weighted_lpl_spect_rad <- plot_distributions(nnm$modresultsweighted,"lpl_spect_rad_weighted",paste("Lpl Weighted spectral radius",trsuff),nnm$wmag$lpl_spect_rad_weighted,netw,fract_dummy,nestedv=nnm$pnw$weightednstlplspect_rad)
    
    pd_weighted_spect_rad <- plot_distributions(nnm$modresultsweighted,"spect_rad_weighted",paste("Weighted spectral radius",trsuff),nnm$wmag$spect_rad_weighted,netw,fract_dummy,nestedv=nnm$pnw$weightednstspect_rad)
  }
  
  pd_spect_rad <- plot_distributions(nnm$modresults,"spect_rad","Spectral radius",nnm$bmag$adj_spect$values[1],netw,fract_dummy,nlinks=nnm$num_links,nestedv=nnm$pnms$nstspect_rad,hypernestedv=nnm$pnm$nstspect_rad)
  pd_lpl_spect_rad <- plot_distributions(nnm$modresults,"lpl_spect_rad","Lpl Spectral radius",nnm$bmag$lpl_spect$values[1],netw,fract_dummy,nestedv=nnm$pnms$nstlplspect_rad,hypernestedv=nnm$pnm$nstlplspect_rad)
  pd_lpl_energy <- plot_distributions(nnm$modresults,"lpl_energy","Laplacian energy",nnm$bmag$lpl_energy,netw,fract_dummy,nestedv=nnm$pnms$nstlpl_energy,hypernestedv=nnm$pnm$nstlpl_energy)
  pd_energy <- plot_distributions(nnm$modresults,"adj_energy","Adj. Energy",nnm$bmag$adj_energy,netw,fract_dummy,nestedv=nnm$pnms$nstadj_energy,hypernestedv=nnm$pnm$nstadj_energy)
  pd_algebraic_connectivity <- plot_distributions(nnm$modresults,"algebraic_connectivity","Alg. connectivity",sort(nnm$bmag$lpl_spect$values)[2],netw,fract_dummy,nestedv=nnm$pnms$nstalgebraic_connectivity,hypernestedv=nnm$pnm$nstalgebraic_connectivity)
  modresults <- nnm$modresults
  if (nnm$weighted_network){
    wmed <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_spect_rad$pimage | pd_lpl_energy$pimage | pd_algebraic_connectivity$pimage) 
    winf <- (pd_weighted_adj_energy$pimage |  pd_weighted_spect_rad$pimage | pd_weighted_lpl_spect_rad$pimage | pd_weighted_lpl_energy$pimage )
    wtot <- (wmed/ winf) + plot_layout(heights = c(0.5,0.5))
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5)))
    plwidth <- 21
    plheight <- 10
  } else {
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum)
    winf <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_spect_rad$pimage | pd_lpl_energy$pimage| pd_algebraic_connectivity$pimage)
    wtot <- winf
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
    plwidth <- 21
    plheight <- 5
  }
  sweight <- "WEIGHTED"
  if (!nnm$weighted_network)
    sweight <- "BINARY"
  nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_",sweight))
  png(paste0(nfile,".png"),width=plwidth*ppi,height=plheight*ppi,res=ppi)
  print(wtot)
  dev.off()
  if (nnm$weighted_network){
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum | pl_lpl_weighted_spectrum)
    
    wmed <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_spect_rad$pimage | pd_lpl_energy$pimage | pd_algebraic_connectivity$pimage) 
    winf <- (pd_weighted_adj_energy$pimage | pd_weighted_spect_rad$pimage | pd_weighted_lpl_spect_rad$pimage | pd_weighted_lpl_energy$pimage  )
    wtot <- (wsup / wmed/ winf) + plot_layout(heights = c(0.33,0.33,0.33))
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5)))
    plwidth <- 17
    plheight <- 13
  } else {
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum)
    winf <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_spect_rad$pimage | pd_lpl_energy$pimage)
    wtot <- (wsup / winf) + plot_layout(heights = c(0.5,0.5))
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
    plwidth <- 17
    plheight <- 10
  }
  sweight <- "WEIGHTED"
  if (!nnm$weighted_network)
    sweight <- "BINARY"
  nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_NESTEDMAGS_",sweight))
  png(paste0(nfile,".png"),width=plwidth*ppi,height=plheight*ppi,res=ppi)
  print(wtot)
  dev.off()
}

# Plots the ziggurat diagram of the model
plot_null_model_zigg <- function(modeldata,modeltype,netname,dirn,plotdir){
  print(paste("Plotting ziggurat",netname,modeltype))
  filenull <- paste0(gsub(".csv","",netname),"_NULL_",modeltype,".csv")
  result_analysis <- analyze_network(directory = dirn, filenull, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  pgr_null <- ziggurat_graph(dirn,filenull,plotsdir=plotdir,print_to_file = TRUE,show_title = FALSE,weighted_links = "log10")
}

plot_ziggurats <- function(nnm,netw,dirnulls,dirdata,odir){
  pdir <- paste0(odir,"zigs/")
  # Plot original network
  
  filenull <- paste0(gsub(".csv","",netw),".csv")
  result_analysis <- analyze_network(directory = dirdata, filenull, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  pgr_null <- ziggurat_graph(dirdata,filenull,plotsdir=pdir,print_to_file = TRUE,show_title = FALSE,weighted_links = "log10")
  
  plot_null_model_zigg(nnm$pnm$nstmodel,"HNESTED",netw,dirnulls,pdir)
  plot_null_model_zigg(nnm$pnms$nstmodel,"NESTED",netw,dirnulls,pdir)
  for (tmodel in nnm$mnames)
    plot_null_model_zigg(nnm$model_full[[tmodel]],tmodel,netw,dirnulls,pdir)
  if (nnm$weighted_network)
    plot_null_model_zigg(ceiling(weightednstmodel),"WNESTED",netw,dirnulls,pdir)
}
