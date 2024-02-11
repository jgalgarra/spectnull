

plot_distr_null <- function(df,nvalue,networkname="",title="",nestedvalue=""){
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
  plot <- plot+geom_vline(xintercept = nvalue,color = "blue", size=0.5,alpha=0.5)+
    ggtitle(title)
  
  if (nestedvalue!="")
    plot <- plot+geom_vline(xintercept = nestedvalue,color = "red", size=0.5,alpha=0.5,linetype="dotted")+
    ggtitle(title)
  
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

# Plots the ziggurat diagram of the model
plot_null_model <- function(modeldata,modeltype,netname,dirn,plotdir,ptfile){
  filenull <- paste0(gsub(".csv","",netname),"_NULL_",modeltype,".csv")
  print(paste("FILENULL",filenull))
  write.csv(modeldata,paste0(dirnulls,filenull))
  result_analysis <- analyze_network(directory = dirn, filenull, guild_a = "Plant", guild_b = "Pollinator", only_NODF = TRUE)
  pgr_null <- ziggurat_graph(dirn,filenull,plotsdir=plotdir,print_to_file = ptfile,show_title = FALSE,weighted_links = "log10")
}

# Plots the probability distributions of all models for a particular magnitude and network
plot_distributions <- function(mresults,magnitude,stitle,vtitle,nname,fdummy,nestedv=""){
  mdls <- names(mresults)
  col <- which(names(mresults[[1]])==magnitude)
  df_nulls <- data.frame("vals"=mresults[[1]][,col],"type"=mdls[[1]])
  for (i in 2:length(mdls))
    df_nulls <- rbind(df_nulls,data.frame("vals"= mresults[[mdls[i]]][,col],"type"=mdls[i]))
  
  pimage <- plot_distr_null(df_nulls,vtitle,networkname=nname,title=sprintf("%s %.2f",stitle,vtitle),nestedvalue = nestedv)
  calc_values <- list("df_nulls"= df_nulls, "pimage" = pimage)
  return(calc_values)
}

plot_all_distr <- function(mnames,bmag,wmag){
  mspec <- mnames
  df_adj_model <- lapply(1:length(mspec), function(x) c())
  names(df_adj_model) <- mspec
  df_lpl_model <- df_adj_model
  df_lpl_weighted_model <- df_adj_model
  
  df_adj_model[["NETWORK"]] <- store_spect_values(bmag$adj_spect$values,"NETWORK")
  df_lpl_model[["NETWORK"]] <- store_spect_values(bmag$lpl_spect$values,"NETWORK")
  if (weighted_network)
    df_lpl_weighted_model[["NETWORK"]] <- store_spect_values(wmag$lpl_weighted_spect$values,"NETWORK")
  for (m in mspec[2:length(mspec)])
  {
    df_adj_model[[m]] <- store_spect_values(specresults[[m]]$null_spect$values,m)
    df_lpl_model[[m]] <- store_spect_values(specresults[[m]]$lpl_spect_null$values,m)
    if (weighted_network)
      df_lpl_weighted_model[[m]] <- store_spect_values(specresults[[m]]$lpl_weighted_spect_null$values,m)
  }
  if (plotzigs)
    ziggurat_graph(datadir,netw,plotsdir=odir,print_to_file = plottofile,show_title = FALSE,weighted_links = "log10")
  pl_adj_spectrum <- plot_spectr(do.call(rbind,df_adj_model),title="Adjacency spectrum")
  pl_lpl_spectrum <- plot_spectr(do.call(rbind,df_lpl_model),title="Laplacian spectrum")
  if (weighted_network)
    pl_lpl_weighted_spectrum <- plot_spectr(do.call(rbind,df_lpl_weighted_model),title="Weighted Laplacian spectrum")
  
  for (tmodel in mnames)
    #if (found_model_full[[tmodel]])
      if (plotzigs)
        plot_null_model(model_full[[tmodel]],tmodel,netw,dirnulls,odir,plottofile)
  if (plotzigs){
    plot_null_model(nstmodel,"NESTED",netw,dirnulls,odir,plottofile)
    if (weighted_network)
      plot_null_model(ceiling(weightednstmodel),"WNESTED",netw,dirnulls,odir,plottofile)
  }
  # save_nested_model(pnm$nstmodel,netw,dirnulls,"NESTED")
  # if (weighted_network)
  #   save_nested_model(pnw$weightednstmodel,netw,dirnulls,"WNESTED")
  modresultsweighted <- modresults
  for (mname in c("RND","BSHUFFLE","BVAZ","SYTR"))
    modresultsweighted <- remove_model_data(modresultsweighted,mname)
  
  if (weighted_network)
    for (mname in c("SHUFFLE","VAZ","SWAP"))
      modresults <- remove_model_data(modresults,mname)
  else
    for (mname in c("BSHUFFLE","BVAZ"))
      modresults <- remove_model_data(modresults,mname)
  for (mname in c("WRND","WSYTR","MGEN","PATEFIELD"))
    modresults <- remove_model_data(modresults,mname)
  
  if (weighted_network){
    if (weightrf=="none")
      trsuff <- ""
    else
      trsuff <- paste0("[",weightrf,"]")
    pd_weighted_adj_energy <- plot_distributions(modresultsweighted,"adj_weighted_energy",paste("Weighted Energy",trsuff),wmag$adj_weighted_energy,netw,fract_dummy,nestedv=pnw$weightednstadj_energy)
    pd_weighted_lpl_energy <- plot_distributions(modresultsweighted,"lpl_weighted_energy",paste("Weighted Laplacian energy",trsuff),wmag$lpl_weighted_energy,netw,fract_dummy,nestedv=pnw$weightednstlpl_energy)
    pd_weighted_spect_rad <- plot_distributions(modresultsweighted,"spect_rad_weighted",paste("Weighted spectral radius",trsuff),wmag$spect_rad_weighted,netw,fract_dummy,nestedv=pnw$weightednstspect_rad)
  }
  pd_spect_rad <- plot_distributions(modresults,"spect_rad","Spectral radius",bmag$adj_spect$values[1],netw,fract_dummy,nestedv=pnm$nstspect_rad)
  pd_lpl_energy <- plot_distributions(modresults,"lpl_energy","Laplacian energy",bmag$lpl_energy,netw,fract_dummy,nestedv=pnm$nstlpl_energy)
  pd_energy <- plot_distributions(modresults,"adj_energy","Energy",bmag$adj_energy,netw,fract_dummy,nestedv=pnm$nstadj_energy)
  if (weighted_network){
    wmed <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage) 
    winf <- (pd_weighted_adj_energy$pimage | pd_weighted_spect_rad$pimage | pd_weighted_lpl_energy$pimage )
    wtot <- (wmed/ winf) + plot_layout(heights = c(0.5,0.5))
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5)))
    plwidth <- 13
    plheight <- 10
  } else {
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum)
    winf <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage)
    wtot <- winf
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
    plwidth <- 13
    plheight <- 5
  }
  sweight <- "WEIGHTED"
  if (!weighted_network)
    sweight <- "BINARY"
  nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_",sweight))
  png(paste0(nfile,".png"),width=plwidth*ppi,height=plheight*ppi,res=ppi)
  print(wtot)
  dev.off()
  # if (weighted_network){
  #   if (weightrf=="none")
  #     trsuff <- ""
  #   else
  #     trsuff <- paste0("[",weightrf,"]")
  #   pd_weighted_adj_energy <- plot_distributions(modresultsweighted,"adj_weighted_energy",paste("Weighted Energy",trsuff),wmag$adj_weighted_energy,netw,fract_dummy,nestedv=pnw$weightednstadj_energy)
  #   pd_weighted_lpl_energy <- plot_distributions(modresultsweighted,"lpl_weighted_energy",paste("Weighted Laplacian energy",trsuff),wmag$lpl_weighted_energy,netw,fract_dummy,nestedv=pnw$weightednstlpl_energy)
  #   pd_weighted_spect_rad <- plot_distributions(modresultsweighted,"spect_rad_weighted",paste("Weighted spectral radius",trsuff),wmag$spect_rad_weighted,netw,fract_dummy,nestedv=pnw$weightednstspect_rad)
  # }
  # pd_spect_rad <- plot_distributions(modresults,"spect_rad","Spectral radius",bmag$adj_spect$values[1],netw,fract_dummy,nestedv=pnm$nstspect_rad)
  # pd_energy <- plot_distributions(modresults,"adj_energy","Energy",bmag$adj_energy,netw,fract_dummy,nestedv=pnm$nstadj_energy)
  # pd_lpl_energy <- plot_distributions(modresults,"lpl_energy","Laplacian energy",bmag$lpl_energy,netw,fract_dummy,nestedv=pnm$nstlpl_energy)
  if (weighted_network){
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum | pl_lpl_weighted_spectrum)
    wmed <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage) 
    winf <- (pd_weighted_adj_energy$pimage | pd_weighted_spect_rad$pimage | pd_weighted_lpl_energy$pimage )
    wtot <- (wsup / wmed/ winf) + plot_layout(heights = c(0.33,0.33,0.33))
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5)))
    plwidth <- 13
    plheight <- 13
  } else {
    wsup <- (pl_adj_spectrum | pl_lpl_spectrum)
    winf <- (pd_energy$pimage | pd_spect_rad$pimage | pd_lpl_energy$pimage)
    wtot <- (wsup / winf) + plot_layout(heights = c(0.5,0.5))
    plot_annotation(title = nname,
                    theme = theme(plot.title = element_text(size = 16,hjust=0.5))) 
    plwidth <- 13
    plheight <- 10
  }
  sweight <- "WEIGHTED"
  if (!weighted_network)
    sweight <- "BINARY"
  nfile <- paste0(odir,paste0(nname,"_ALLPLOTS_NESTEDMAGS_",sweight))
  png(paste0(nfile,".png"),width=plwidth*ppi,height=plheight*ppi,res=ppi)
  print(wtot)
  dev.off()
}  