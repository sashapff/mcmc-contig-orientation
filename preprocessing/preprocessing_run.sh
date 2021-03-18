#!/bin/sh
# one hour timelimit:
#SBATCH --time 1:00:00
# default queue, 32 processors (two nodes worth)
#SBATCH -p defq -n 32
#SBATCH -o output.txt


module load python

python main.py