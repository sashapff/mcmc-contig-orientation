#!/bin/sh
#SBATCH --time 00:05:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt

module load Anaconda2

python main.py
