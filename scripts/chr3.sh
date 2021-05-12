#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug-cpu
#SBATCH -o output/chr3.txt

python chr3.py

