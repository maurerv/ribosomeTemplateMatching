#!/bin/bash

#SBATCH --job-name=process_coordinates
#SBATCH --output=process_coordinates_%A_%a.out
#SBATCH --error=process_coordinates_%A_%a.err
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=20G
#SBATCH --partition=htc-el8

source $HOME/.bashrc

dirs=( $BASE_DIR/*/ )
dir=${dirs[$SLURM_ARRAY_TASK_ID-1]}
dir=${dir%*/}

echo "Processing: $dir"

module load R/4.2.0-foss-2021b
Rscript process_coordinates.R $dir
