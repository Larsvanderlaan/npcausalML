.libPaths( c( .libPaths(), "~/Rlibs2") )
setwd("~/LRRsims")
print(getwd())
nsims = 1000
library(npcausalML)
#library(future)
#plan(multisession, workers = 16)

