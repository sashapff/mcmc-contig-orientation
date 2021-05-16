from tools.load import get_contigs_and_pairs
from tools.prob import density, toy_density
from tools.tools import get_longest_contig, filter_pairs, log_likelihood, get_distance_one_contig
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

if __name__ == "__main__":
    print("Start!")

    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/distance'
    # path_to_output = f'/Users/alexandra/bioinf/mcmc/output/{min_contig_length_name}'

    chr_ind = '4'
    path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"
    path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"
    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"

    # # path_layout = "/Users/alexandra/bioinf/mcmc/data/chr1.layout.txt"
    # path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"
    # #
    # # path_lens = "/Users/alexandra/bioinf/mcmc/data/comp18_lens.tsv"
    # path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"
    # #
    # # path_pairs = "/Users/alexandra/bioinf/mcmc/data/pairs18.txt"
    # path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig_pairs, longest_contig = get_contigs_and_pairs(path_layout, path_lens,
                                                                                            path_pairs,
                                                                                            long_contig=True,
                                                                                            # min_len=0,
                                                                                            from_one_contig=True)

    # P = toy_density(longest_contig_pairs)
    P, f = density(longest_contig_pairs)

    ind = (pairs[:, 1] >= pairs[:, 3])
    pairs[ind, 1], pairs[ind, 3] = pairs[ind, 3], pairs[ind, 1]

    D = 1000
    # d = 5
    left = longest_contig.length // 2 - D // 2
    right = left + D

    filtered_pairs = filter_pairs(pairs, id_contig[longest_contig.name], left, right)
    filtered_pairs[:, 3] -= right

    log_likelihood_arr = []
    log_likelihood_range = [1, 10, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000]

    for d in log_likelihood_range:
        filtered_pairs[:, 6] = d
        coeff = np.log(quad(np.exp(P), d, 10**10))
        print(f'd={d}, coeff={coeff}')
        ll = P(get_distance_one_contig(filtered_pairs, left)).sum() / coeff
        log_likelihood_arr.append(ll)
        print(f'likelihood for d={d} is {ll}')

    plt.plot(log_likelihood_range, log_likelihood_arr)
    plt.xlabel('d, distance estimate')
    plt.ylabel('log likelihood')
    plt.legend()
    plt.xscale('log')
    plt.savefig(f'{path_to_output}/log_likelihood.png')




