# HiC_MCMC
The orientation of contigs using Hi-C data

## Purpose of the project:
Existing scaffolding algorithms using Hi-C data is good in finding orders of contigs. However, this builds contain sufficient numbers of errors at orientations. We apply a probability model for resolving this problem.

## Usage:
I. load data with function get_contigs_and_pairs

II. estimate density of distance between pieces of reads

III. run the MSMS algorithm

(see HiC_MCMC/main)

## Dependencies:
* Python 3.7

Python dependencies can be installed with pip:
 
 `
 pip install -r requirements.txt
 `

## Input data:
Data has to contain information about HiC reads in to follow a format
1) pairs.txt

| * | name of contig which contains the first piece of read | position first piece of read in contigs | name of contig which contains the second piece of read | position second  piece of read in contigs | * | * |
2) len.tsv

| name of contig | his length |
3) layout.txt has to contain order of contigs with the current orientation

(see data)

## Links 
1. MCMC algorithms
http://statweb.stanford.edu/~cgates/PERSI/papers/MCMCRev.pdf
