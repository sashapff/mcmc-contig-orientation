#!/bin/sh
#SBATCH --time 8:00:00
#SBATCH -p debug
#SBATCH -o output.txt

python main_orientation.py
#python main_ordering.py
