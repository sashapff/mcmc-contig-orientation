#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug
#SBATCH -o output/ordering_exp_change_pos.txt

python main_ordering.py
