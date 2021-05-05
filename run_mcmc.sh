#!/bin/sh
#SBATCH --time 8:00:00
#SBATCH -p tiny
#SBATCH -o output.txt

#python main_orientation.py
python main_ordering.py
