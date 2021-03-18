#!/bin/sh
# one hour timelimit:
#SBATCH --time 6:00:00
#SBATCH -p tiny
#SBATCH -o output.txt


module load python

python main.py
