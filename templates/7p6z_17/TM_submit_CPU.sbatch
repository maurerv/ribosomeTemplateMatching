#!/bin/bash
#SBATCH --mem 200G
#SBATCH -p htc-el8
#SBATCH -N 1
#SBATCH -n 6
#SBATCH -t 48:00:00
#SBATCH -o /Users/vmaurer/src/ribosomeSpheres/templates/7p6z_17/slurm.%N.%j.out
#SBATCH -e /Users/vmaurer/src/ribosomeSpheres/templates/7p6z_17/slurm.%N.%j.err
module purge
module load PyTom/1.0b-foss-2021a-CUDA-11.3.1
time mpirun -n 6 ` which localization.py ` -j /Users/vmaurer/src/ribosomeSpheres/templates/7p6z_17/TM_job.xml -x 3 -y 3 -z 3
