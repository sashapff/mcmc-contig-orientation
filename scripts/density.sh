#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug
#SBATCH -o density_output.txt

python density.py
