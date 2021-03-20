#!/bin/sh
#SBATCH --time 1-00:00:00
#SBATCH -p short
#SBATCH -o output.txt

echo "CHROMOSOME 2"
python main.py 2

echo "CHROMOSOME 3"
python main.py 3

echo "CHROMOSOME 4"
python main.py 4

echo "CHROMOSOME 5"
python main.py 5

echo "CHROMOSOME 6"
python main.py 6

echo "CHROMOSOME 7"
python main.py 7

echo "CHROMOSOME 8"
python main.py 8

echo "CHROMOSOME 9"
python main.py 9

echo "CHROMOSOME 10"
python main.py 10

echo "CHROMOSOME 11"
python main.py 11

echo "CHROMOSOME 12"
python main.py 12

echo "CHROMOSOME 13"
python main.py 13

echo "CHROMOSOME 14"
python main.py 14

echo "CHROMOSOME 15"
python main.py 15

echo "CHROMOSOME 16"
python main.py 16

echo "CHROMOSOME 17"
python main.py 17

echo "CHROMOSOME 18"
python main.py 18

echo "CHROMOSOME 19"
python main.py 19

echo "CHROMOSOME 20"
python main.py 20

echo "CHROMOSOME 21"
python main.py 21

echo "CHROMOSOME 22"
python main.py 22

echo "CHROMOSOME MT"
python main.py MT

echo "CHROMOSOME X"
python main.py X

