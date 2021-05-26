#!/bin/sh
#SBATCH --time 1:00:00
#SBATCH -p debug
#SBATCH -o output/lengths.txt

module load python

python preprocess/filter_reads.py
python preprocess/get_lengths.py
