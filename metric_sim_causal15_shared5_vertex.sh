#!/bin/bash
#
#SBATCH -t 60-00:00 # Runtime in D-HH:MM
#SBATCH -n 20 # number of cores used 
#SBATCH -o sim.out # File to which STDOUT will be written
#SBATCH -e sim.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=timothy_sudijono@brown.edu # Email to which notifications will be sent
#SBATCH --J=sinatra_metric_sim # name of the job

nproc=$(($SLURM_JOB_CPUS_PER_NODE*$SLURM_NNODES))
echo $nproc threads
export OMP_NUM_THREADS=$nproc

module load gcc
module load lapack
module load openblas
module load R/3.5

Rscript metric_sim_causal15_shared5.R