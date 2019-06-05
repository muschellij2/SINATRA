#!/bin/bash
#SBATCH --job-name=timing
#SBATCH -t 96:00:00# Runtime in HH:MM:SS
#SBATCH -N 1 # number of nodes used
#SBATCH -n 1 # number of cores used 
#SBATCH --mem=350g
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=timothy_sudijono@brown.edu # Email to which notifications will be sent
#SBATCH --array=0-35%40

nproc=$(($SLURM_JOB_CPUS_PER_NODE*$SLURM_NNODES))
#echo $nproc threads
export OMP_NUM_THREADS=$nproc

module load gcc
module load lapack
module load openblas
module load R/3.5.2

Rscript R_scripts/Package_Setup.R

SHAPES_PER_CLASS=(25 50 75)
NUM_CONES=(25 50 75)
DIR_PER_CONE=(5 10)
CURVE_LENGTH=(25 50)

for i in {0..2}
 do
 	for j in {0..2}
 	do
 		for k in {0,1}
 		do
 			for l in {0,1}
 			do
				task_id=$(( $(( 2 * $(( $(( 2 * $(( $(( 3 * $(( $i )) )) + $j )) )) + $k )) )) + $l ))
				shapes=${SHAPES_PER_CLASS[$i]}
				numcones=${NUM_CONES[$j]}
				dirpercone=${DIR_PER_CONE[$k]}
				curvelength=${CURVE_LENGTH[$l]}


				if [ "$task_id" = "$SLURM_ARRAY_TASK_ID" ] 
				then 
					Rscript --vanilla R_scripts/Sphere_Simulations/timing_simulation.R $shapes $numcones $dirpercone $curvelength
				fi
			done
		done
	done
done