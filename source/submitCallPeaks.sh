#!/bin/bash
BASE_DIR="/g/kosinski/vmaurer/ribosomePaper/templates"
num_dirs=$(find $BASE_DIR -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=1-$num_dirs callPeaks.sh