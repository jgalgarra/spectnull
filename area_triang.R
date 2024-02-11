area_triangle <- function(base)
{
  area <- base*(base+1)/2
}

basetriangle <- function(area){
  # area <- base * (base+1) / 2 -> b = (-1 + sqrt(1+8*area))/2
  return((-1 + sqrt(1+8*area))/2)
}