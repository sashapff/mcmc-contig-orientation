#!/bin/sh
#SBATCH --time 01:00:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt

python main.py
