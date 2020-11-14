import numpy as np
import pandas as pd
from load import get_contigs_and_pairs
from prob import density, toy_density
from model import MCMC
from tools import get_orientation, simulation


if __name__ == "__main__":
    simulation(100, 5000, 5)
    path_layout = "../data/simulation.layout.txt"
    path_lens = "../data/simulation.lens.tsv"
    path_pairs = "../data/simulation.pairs.txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig, in_contigs = get_contigs_and_pairs(path_layout, path_lens, path_pairs,
                                                                                  long_contig=True,
                                                                                  all_contigs=True,
                                                                                  min_len=0)

    print("Estimation of density...")
    P, f = toy_density(longest_contig)
    print("Estimation of density is done")

    get_orientation([0 for i in range(len(contigs))], pairs, contigs)
    print("MCMC is running...")
    MCMC(pairs, contigs, P, 500, True, n_chains=1)
    print("Have found follow orientations:", [contigs[i].o for i in range(len(contigs))])
