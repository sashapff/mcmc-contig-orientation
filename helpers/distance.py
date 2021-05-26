import matplotlib.pyplot as plt
import numpy as np

from utils.load import get_contigs_and_pairs
from utils.prob import normalize, estimate_density
from utils.tools import filter_pairs, get_distance_one_contig

if __name__ == "__main__":
    longest_contig_len = 20_000

    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/distance'

    chr_ind = '4'
    path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"
    path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"
    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"
    # longest_contig
    pairs, contigs, id_contig, longest_contig_pairs, longest_contig = get_contigs_and_pairs(path_layout, path_lens,
                                                                                            path_pairs,
                                                                                            long_contig=True,
                                                                                            min_len=0,
                                                                                            from_one_contig=True)

    longest_contig_len = longest_contig.length

    ind = (pairs[:, 1] >= pairs[:, 3])
    pairs[ind, 1], pairs[ind, 3] = pairs[ind, 3], pairs[ind, 1]

    D = 1000
    left = longest_contig.length // 2 - D // 2
    right = left + D

    filtered_pairs = filter_pairs(pairs, id_contig[longest_contig.name], left, right)

    P, f = estimate_density(filtered_pairs)
    P = normalize(P, 0, np.inf)

    p_range = range(1000)
    plt.title('The approximate log density')
    plt.plot(p_range, [P(x) for x in p_range])
    plt.xlabel('x, distance')
    plt.ylabel('p(x), approximate log density')
    plt.show()
    plt.clf()

    filtered_pairs[:, 3] -= right

    log_likelihood_arr = []
    log_likelihood_range = range(0, longest_contig_len, 100)

    for d in log_likelihood_range:
        filtered_pairs[:, 6] = d
        P_norm = normalize(P, d, np.inf)
        ll = P_norm(get_distance_one_contig(filtered_pairs, left)).sum()
        log_likelihood_arr.append(ll)
    plt.plot(log_likelihood_range, log_likelihood_arr)
    plt.xlabel('d, distance estimate')
    plt.ylabel('log likelihood')
    plt.legend()
    plt.savefig(f'{path_to_output}/log_likelihood.png')
