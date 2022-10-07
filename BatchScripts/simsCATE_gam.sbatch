#!/bin/bash
#SBATCH --job-name Lars_simsCATE_gam%j        # Set a name for your job. This is especially useful if you have multiple jobs queued.
#SBATCH --partition medium       # Slurm partition to use
#SBATCH --ntasks 10            # Number of tasks to run. By default, one CPU core will be allocated per task
#SBATCH --time 10-00:00          # Wall time limit in D-HH:MM
#SBATCH --mem-per-cpu=10000       # Memory limit for each tasks (in MB)
#SBATCH -o outfiles/simsCATE_gam.out      # File to which STDOUT will be written
#SBATCH -e outfiles/simsCATE_gam.err      # File to which STDERR will be written
#SBATCH --mail-type=NONE        # Type of email notification- NONE,BEGIN,END,FAIL,ALL
#SBATCH --mail-user=lvdlaan@uw.edu # Email to which notifications will be sent
echo "LOADING R"
module load R
echo "R LOADED"
Rscript -e 'source("R_setup.R"); source("~/LRRsims/FinalSimulationCode/simsCATE_gam.R");
source("~/LRRsims/FinalSimulationCode/performancePlot_gam.R")'
