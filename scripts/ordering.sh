#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug
#SBATCH -o output/ordering_exp_swap.txt

python ordering.py
