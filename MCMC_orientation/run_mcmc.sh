#!/bin/sh
#SBATCH --time 00:05:00
#SBATCH -p debug-gpu
#SBATCH -o output.txt

conda create --name pandas
conda activate pandas

python main.py
