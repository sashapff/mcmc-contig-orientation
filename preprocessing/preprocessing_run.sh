#!/bin/sh
#SBATCH --time 1:00:00
#SBATCH -p debug
#SBATCH -o output.txt


module load python

python get_contig_length.py
