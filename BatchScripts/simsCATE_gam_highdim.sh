#!/bin/usr/env bash

nsims=2500
export R_LIBS=~/Rlibs2
export R_LIBS_USER=~/Rlibs2
for n in 500 1000 2500 5000
do
  for pos in "TRUE" "FALSE"
  do
    for hard in "TRUE" "FALSE"
    do
      sbatch  --export=n=$n,pos=$pos,hard=$hard ~/LRRsims/BatchScripts/simsCATE_gam_highdim.sbatch
    done
  done
done
