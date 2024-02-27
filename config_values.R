debugpref <- "DEB"
ppi <- 100
options( warn = -1 )
datadir <- "data/"
lweightrf <- c("sqrt","ln","none")
dbaseb <- "smodels/"
odirb <- "plots/"
rdirb <- "results/"
dirnullsb <- "nullmatrix/"
# dir.create(odirb, showWarnings = FALSE)
# dir.create(rdirb, showWarnings = FALSE)
# dir.create(dirnullsb, showWarnings = FALSE)

create_dirs <- function(weightrf){
  dbase <<- paste0(debugpref,dbaseb)
  dir.create(dbase, showWarnings = FALSE)
  dbase <<- paste0(dbase,weightrf,"/")
  dir.create(dbase, showWarnings = FALSE)
  odir <<- paste0(dbase,odirb)
  rdir <<- paste0(dbase,rdirb)
  dirnulls <<- paste0(dbase,dirnullsb)
  dir.create(odir, showWarnings = FALSE)
  dir.create(rdir, showWarnings = FALSE)
  dir.create(dirnulls, showWarnings = FALSE)
}
