#!/bin/bash
#SBATCH --job-name=sphere_sim_ROC
#SBATCH -t 48:00:00# Runtime in HH:MM:SS
#SBATCH -n 20 # number of cores used 
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=timothy_sudijono@brown.edu # Email to which notifications will be sent

nproc=$(($SLURM_JOB_CPUS_PER_NODE*$SLURM_NNODES))
echo $nproc threads
export OMP_NUM_THREADS=$nproc

module load gcc
module load lapack
module load openblas
module load R/3.5.1

Rscript R_scripts/Package_Setup.R

CAUSAL=(1 2 4 2 3 6)
SHARED=(2 3 6 1 2 4)

Rscript --vanilla R_scripts/Sphere_Simulations/Generate_ROC.R 4 6 10




