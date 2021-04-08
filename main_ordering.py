from tools.load import get_contigs_and_pairs
from ordering.model import MCMC
from ordering.tools import get_ordering
from tools.prob import density, toy_density
import numpy as np

if __name__ == "__main__":
    print("Start!")
    # chr_ind = sys.argv[1]
    # Example

    path_layout = "/Users/alexandra/bioinf/mcmc/data/chr1.layout.txt"
    # path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"

    path_lens = "/Users/alexandra/bioinf/mcmc/data/comp18_lens.tsv"
    # path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"

    path_pairs = "/Users/alexandra/bioinf/mcmc/data/pairs18.txt"
    # path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig, longest_contig_name = get_contigs_and_pairs(path_layout, path_lens, path_pairs, long_contig=True)
    # pairs, contigs, id_contig, longest_contig, longest_contig_name, in_contigs = get_contigs_and_pairs(path_layout, path_lens, path_pairs,
    #                                                                               long_contig=True,
    #                                                                               all_contigs=True,
    #                                                                               min_len=0)

    correct_contigs = [contig.pos for contig in contigs]

    print("Estimation of density...")
    P, f = density(longest_contig)
    # P, f = toy_density(longest_contig)
    print("Estimation of density is done")

    print("MCMC is running...")
    get_ordering(np.random.choice(len(contigs), len(contigs), replace=False), pairs, contigs)
    MCMC(pairs, contigs, P, 5000, n_chains=1)
    print("Have found follow ordering:", [contigs[i].pos for i in range(len(contigs))])

    # with open("/Users/alexandra/bioinf/mcmc/data/final1.layout.txt", "w") as file:
    #     sign = lambda x: "+" if x == 1 else "-"
    #     file.write(",".join([contig.name + sign(contig.o) for contig in sorted(contigs, key=lambda contig: contig.pos)]))
    # print("Result has been saved")

    # print("Count correctness...")
    # correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]] for contig in contigs]).sum()
    # print(f"{correct_number}/{len(contigs)} contigs were oriented correctly ({correct_number / len(contigs) * 100}%)")
    #
    # with open("/lustre/groups/cbi/Users/aeliseev/aivanova/data/stat." + chr_ind + ".txt", "w") as file:
    #     file.write(f"CHROMOSOME {chr_ind}\n")
    #     file.write(f"{correct_number}/{len(contigs)} contigs were oriented correctly ({correct_number / len(contigs) * 100}%)\n")
