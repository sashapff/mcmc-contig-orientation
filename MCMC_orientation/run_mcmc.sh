#!/bin/sh
#SBATCH --time 1:00:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt


module load python
module load pandas

python --version
