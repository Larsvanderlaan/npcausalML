.libPaths( c( .libPaths(), "~/Rlibs2") )
setwd("~/LRRsims")
print(getwd())
nsims = 2500
library(npcausalML)
#library(future)
#plan(multisession, workers = 16)

