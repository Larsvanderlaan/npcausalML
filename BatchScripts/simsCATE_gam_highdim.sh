
filename="simsCATE_gam_highdim"
nsims=2500
export R_LIBS=~/Rlibs
export R_LIBS_USER=~/Rlibs
sbatch  --export=filename=$filename,nsims=$nsims simsCATE_gam_highdim.sbatch
