from tools.load import get_contigs_and_pairs
from tools.prob import density, toy_density
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/density'

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

    P, f = density(longest_contig_pairs)

    p_range = range(200000)
    plt.title('The approximate log density')
    plt.plot(p_range, [P(x) for x in p_range])
    plt.xlabel('x, distance')
    plt.ylabel('p(x), approximate log density')
    plt.savefig(f'{path_to_output}/approximate_200k.png')

    plt.clf()

    p_range = range(5000)
    plt.title('The approximate log density')
    plt.plot(p_range, [P(x) for x in p_range])
    plt.xlabel('x, distance')
    plt.ylabel('p(x), approximate log density')
    plt.savefig(f'{path_to_output}/approximate_5k.png')

    plt.clf()

    p_range = range(200000)
    plt.title('The estimated log density')
    plt.plot(p_range, f(np.array(p_range)))
    plt.xlabel('x, distance')
    plt.ylabel('f(x), estimated log density')
    plt.savefig(f'{path_to_output}/estimate_200k.png')

    plt.clf()

    p_range = range(5000)
    plt.title('The estimated log density')
    plt.plot(p_range, f(np.array(p_range)))
    plt.xlabel('x, distance')
    plt.ylabel('f(x), estimated log density')
    plt.savefig(f'{path_to_output}/estimate_5k.png')
