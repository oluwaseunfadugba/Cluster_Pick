#!/bin/bash
#SBATCH --partition computeq
#SBATCH --cpus-per-task 1
#SBATCH --array 1-500%50
#SBATCH --time=25-00:00:00
#SBATCH --account=account_name
#SBATCH --job-name=88_20_1

module load matlab/R2018a
matlab -nodisplay -r 'i=${SLURM_ARRAY_TASK_ID};cluster_pick_Driver_HPC(i)'

