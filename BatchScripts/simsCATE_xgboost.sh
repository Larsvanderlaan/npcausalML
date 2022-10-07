
filename="simsCATE_xgboost"
nsims=2500
export R_LIBS=~/Rlibs
export R_LIBS_USER=~/Rlibs
sbatch  --export=filename=$filename,nsims=$nsims simsCATE_xgboost.sbatch
