#!/bin/csh
#SBATCH --time=1:00:00 # walltime, abbreviated by -t
#SBATCH --nodes=1      # number of cluster nodes, abbreviated by -N
#SBATCH -o slurm-%j.out-%N # name of the stdout, using the job number (%j) and the first node (%N)
#SBATCH -e slurm-%j.err-%N # name of the stderr, using job and first node values
#SBATCH --ntasks=1    # number of MPI tasks, abbreviated by -n
# additional information for allocated clusters
#SBATCH --account=lvanderl     # account - abbreviated by -A
#SBATCH --partition=lonepeak  # partition, abbreviated by -p

module load R
Rscript -e 'source("~/LRRsims/BatchScripts/test.R")'
