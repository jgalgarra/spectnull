get_colors <- function(models,colorlist){
  lnames <- names(colorlist)
  lcolors <- colorlist
  for (i in lnames){
    if (!(i %in% models)){
      lnames <- lnames[-which(lnames==i)]
      lcolors <- lcolors[-which(names(lcolors)==i)]
    }
  }
  calc_values <- list("lnames"= lnames, "lcolors" = lcolors)
  return(calc_values)
}


colormodels <- c("NETWORK"="blue","RND"="orange", "WRND"="coral3","B/VAZ"="palegreen4",
                 "B/SHUFFLE"="khaki4",
                 "VAZ"="palegreen3","BVAZ"="palegreen","SHUFFLE"="khaki2",
                 "BSHUFFLE"="khaki3","SWAP"="steelblue4","MGEN"="lightblue",
                 "PATEFIELD"="grey", "SYTR"="pink","WSYTR"="darkorchid1","HNESTED"="red","NESTED"="red4","WNESTED"="darkorchid3")
