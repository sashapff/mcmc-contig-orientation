#!/bin/sh
#SBATCH --time 01:00:00
#SBATCH -p debug-cpu
#SBATCH -o output.txt

echo "CHROMOSOME 2"
python main.py 2

echo "CHROMOSOME 3"
python main.py 3

