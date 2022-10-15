.libPaths( c( "~/Rlibs2", .libPaths()) )
setwd("~/LRRsims")
print(getwd())
nsims = 2500
library(npcausalML)
#library(future)
#plan(multisession, workers = 16)

