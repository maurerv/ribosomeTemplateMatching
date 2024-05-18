#!/bin/bash
export BASE_DIR="/g/kosinski/vmaurer/ribosomePaper/templates_pytme2k"
num_dirs=$(find $BASE_DIR -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=1-$num_dirs callPeaks_pytme.sh