from load import get_contigs_and_pairs
from prob import toy_density
from prob import density
from model import MCMC
from tools import get_orientation
import numpy as np


if __name__ == "__main__":
    # Example
    path = "/Users/alexandra/bioinf/mcmc/data/"

    path_layout = path + "chr1.layout.txt"
    path_lens = path + "comp18_lens.tsv"
    path_pairs = path + "pairs18.txt"

    # path_layout = path + "simulation.layout.txt"
    # path_lens = path + "simulation.lens.tsv"
    # path_pairs = path + "simulation.pairs.txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig = get_contigs_and_pairs(path_layout, path_lens, path_pairs, long_contig=True)
    # pairs, contigs, id_contig, longest_contig, in_contigs = get_contigs_and_pairs(path_layout, path_lens, path_pairs,
    #                                                                               long_contig=True,
    #                                                                               all_contigs=True,
    #                                                                               min_len=0)

    correct_contigs = contigs.copy()

    print("Estimation of density...")
    P, f = density(longest_contig)
    # P, f = toy_density(longest_contig)
    print("Estimation of density is done")

    print("MCMC is running...")
    get_orientation([0 for i in range(len(contigs))], pairs, contigs)
    MCMC(pairs, contigs, P, 100, n_chains=1)
    print("Have found follow orientation:", [contigs[i].o for i in range(len(contigs))])


    with open(path + "final.layout.txt", "w") as file:
        sign = lambda x: "+" if x == 1 else "-"
        file.write(",".join([contig.name + sign(contig.o) for contig in contigs]))
    print("Result has been saved")

    print("Count correctness...")
    correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]].o for contig in contigs]).sum()
    print(f"{correct_number / len(contigs) * 100}% contigs were oriented correctly")