# Contigs ordering using Hi-C data

## Purpose of the project:
Based on collected contigs and a Hi-C signal, it is necessary to assemble the genome, that is, to orient and order the contigs and estimate the distances between them.

## Methods
* Randomly choose the initial state (orientation or order)
* Repeat iteratively
** Making a random change (changing orientation or order)
** If the likelihood function (P) has increased, then we accept the change, otherwise we accept with probability: P(new state)/P(old state)

## Usage:


## Dependencies:
* Python 3.7

Python dependencies can be installed with pip:
 
 `
 pip install -r requirements.txt
 `

## Input data:
Data has to contain information about Hi-C reads in to follow a format

1) pairs.txt

| * | name of contig which contains the first piece of read | position first piece of read in contigs | name of contig which contains the second piece of read | position second  piece of read in contigs | * | * |

2) len.tsv

| name of contig | his length |

3) layout.txt has to contain order of contigs with the current orientation

(see data)

## Links 
1. MCMC algorithms
http://statweb.stanford.edu/~cgates/PERSI/papers/MCMCRev.pdf
