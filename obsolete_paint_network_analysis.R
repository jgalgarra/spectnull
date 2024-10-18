library(ggplot2)
library(forcats)
library(patchwork)
library(dplyr)
library(ggrepel)
library(ggbeeswarm)
library("data.table")
set.seed(122)
normresults <- "results/normalized/"
results <- "results/"
networkmags <- fread(paste0(results,"networkmagnitudes.csv"))
normdistavgs <- fread(paste0(normresults,"meannormdistances.csv"))
networktypes <- fread(paste0(normresults,"networktypes.csv"))


netw <- "M_PL_004"
network_norm_dist <- fread(paste0(normresults,"MINMAX_MODS_",netw,".csv"))
network_nested_indexes <- fread(paste0(results,"NESTED_",netw,".csv"))