#!/bin/bash
#SBATCH --job-name=coneangle
#SBATCH -t 360:00:00# Runtime in HH:MM:SS
#SBATCH -N 1 # number of nodes used
#SBATCH -n 20 # number of cores used 
#SBATCH --mem=300g
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=timothy_sudijono@brown.edu # Email to which notifications will be sent

nproc=$(($SLURM_JOB_CPUS_PER_NODE*$SLURM_NNODES))
#echo $nproc threads
export OMP_NUM_THREADS=$nproc

module load gcc
module load lapack
module load openblas
module load R/3.5.2

Rscript R_scripts/Package_Setup.R

Rscript --vanilla R_scripts/Sphere_Simulations/Generate_ROC_cone_angle.R 
