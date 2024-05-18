#!/bin/bash

export BASE_DIR="/g/kosinski/vmaurer/ribosomePaper/templates_inverted_15k"
num_dirs=$(find $BASE_DIR -mindepth 1 -maxdepth 1 -type d | wc -l)
sbatch --array=1-$num_dirs run_process_coordinates.sh
