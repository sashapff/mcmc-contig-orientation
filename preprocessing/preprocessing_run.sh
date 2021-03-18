#!/bin/sh
# one hour timelimit:
#SBATCH --time 1:00:00
#SBATCH -p defq
#SBATCH -o output.txt


module load python

python main.py