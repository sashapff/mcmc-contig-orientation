#!/bin/sh
#SBATCH --time 1:00:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt

module load Anaconda2
module load python

python main.py
