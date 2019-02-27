#!/bin/bash
#SBATCH --job-name=sphere_sim_ROC
#SBATCH -t 48:00:00# Runtime in HH:MM:SS
#SBATCH -n 20 # number of cores used 
#SBATCH --mem 32G # amount of memomry allocated per node
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=timothy_sudijono@brown.edu # Email to which notifications will be sent
#SBATCH -a 0-12:4



nproc=$(($SLURM_JOB_CPUS_PER_NODE*$SLURM_NNODES))
echo $nproc threads
export OMP_NUM_THREADS=$nproc

module load gcc
module load lapack
module load openblas
module load R/3.5.2

Rscript R_scripts/Package_Setup.R

Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 2 1 10
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 3 2 10
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 6 4 10

Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 2 1 5
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 3 2 5
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 6 4 5

Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 1 2 10
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 2 3 10
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 4 6 10

Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 1 2 5
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 2 3 5
Rscript --vanilla R_scripts/ROC_sim_causal15_shared5_vertex.R 4 6 5