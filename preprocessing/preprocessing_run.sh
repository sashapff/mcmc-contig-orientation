#!/bin/sh
# one hour timelimit:
#SBATCH --time 4:00:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt


module load python

python main.py
