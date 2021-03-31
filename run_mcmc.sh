#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug-cpu
#SBATCH -o output.txt

python main.py

