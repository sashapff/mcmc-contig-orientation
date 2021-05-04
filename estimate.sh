#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug
#SBATCH -o estimate_output.txt

python distance_estimate.py
