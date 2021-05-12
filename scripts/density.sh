#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug
#SBATCH -o output/density.txt

python density.py
