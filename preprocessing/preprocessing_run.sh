#!/bin/sh
#SBATCH --time 4:00:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt


module load python

python get_contig_length.py
