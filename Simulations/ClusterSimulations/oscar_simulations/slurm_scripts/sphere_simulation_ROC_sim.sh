#!/bin/bash
#SBATCH --job-name=sphere_sim_ROC
#SBATCH -t 48:00:00# Runtime in HH:MM:SS
#SBATCH -n 20 # number of cores used 
#SBATCH --mem 32G # amount of memomry allocated per node
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=timothy_sudijono@brown.edu # Email to which notifications will be sent
#SBATCH --array=0-11%20

nproc=$(($SLURM_JOB_CPUS_PER_NODE*$SLURM_NNODES))
echo $nproc threads
export OMP_NUM_THREADS=$nproc

module load gcc
module load lapack
module load openblas
module load R/3.5.2

Rscript R_scripts/Package_Setup.R

CAUSAL=(1 2 4 2 3 6)
SHARED=(2 3 6 1 2 4)

for j in {5,10}
do
	for i in {0..5}
	do
		task_id=$(( $(( 6 * $(( $(( $j/5 - 1 )) )) )) + $i ))
		causalidx=${CAUSAL[$i]}
		sharedidx=${SHARED[$i]}

		if [ "$task_id" = "$SLURM_ARRAY_TASK_ID" ] 
		then 
			Rscript --vanilla R_scripts/Sphere_Simulations/Generate_ROC.R $causalidx $sharedidx $j
		fi
	done
done



