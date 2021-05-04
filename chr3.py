import numpy as np
from orientation.tools import get_orientation
from tools.load import check_reads, get_contigs_and_pairs
from tools.prob import density
from tools.tools import get_distance, log_likelihood
import matplotlib.pyplot as plt

if __name__ == "__main__":
    chr_ind = '3'
    path_pairs = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs{chr_ind}.txt'
    output_path = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/'
    # path_pairs = '/Users/alexandra/bioinf/mcmc/chr3/chr_pairs3.txt'

    # check_reads(path_pairs, output_path, chr_ind)

    path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"
    # path_layout = "/Users/alexandra/bioinf/mcmc/data/chr1.layout.txt"
    # path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"

    # if chr_ind == '9':
    # path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr9.partial_layout.txt"

    path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"
    # path_lens = "/Users/alexandra/bioinf/mcmc/data/comp18_lens.tsv"
    # path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"

    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"
    # path_pairs = "/Users/alexandra/bioinf/mcmc/data/pairs18.txt"
    # path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"

    pairs, contigs, id_contig, longest_contig, longest_contig_name = get_contigs_and_pairs(path_layout, path_lens,
                                                                                           path_pairs,
                                                                                           long_contig=True)
    print(f"Analyse {len(contigs)} contigs..")

    correct_contigs = [contig.o for contig in contigs]

    correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]] for contig in contigs]).sum()

    P, f = density(longest_contig)

    for i in range(2):
        for j in range(2):
            get_orientation([i, j], pairs, contigs)
            correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]] for contig in contigs]).sum()
            print(f"{[contig.o for contig in contigs]}:")
            print(f'\tcorrectness: {correct_number}/{len(contigs)} ({correct_number / len(contigs) * 100}%)')
            print(f'\tlog_likelihood: {log_likelihood(pairs, contigs, P)}')

            distances = get_distance(pairs, contigs)
            lengths = []
            min_left = 10 ** 10
            max_left = 0
            min_right = 10 ** 10
            max_right = 0
            for (k, pair) in enumerate(pairs):
                lengths.append(distances[k])
                min_left = min(min_left, pair[1])
                max_left = max(max_left, pair[1])
                min_right = min(min_right, pair[3])
                max_right = max(max_right, pair[3])

            print(
                f'\tmin pos in left read: {min_left}, max pos in left read: {max_left}, contig total length: {contigs[0].length}')
            print(
                f'\tmin pos in right read: {min_right}, max pos in right read: {max_right}, contig total length: {contigs[1].length}')

            plt.hist(lengths)
            plt.title(f'Orientation [contig.o for contig in contigs], {correct_number / len(contigs) * 100}%')
            plt.xlabel('distance between reads')
            plt.ylabel(f'number of reads')
            plt.savefig(f'{output_path}/plots_check/distances{[contig.o for contig in contigs]}.png')
            plt.clf()
