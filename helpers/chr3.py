import matplotlib.pyplot as plt
import numpy as np

from orientation.tools import get_orientation
from utils.load import check_reads, get_contigs_and_pairs
from utils.prob import density
from utils.tools import get_distance, log_likelihood

if __name__ == "__main__":
    chr_ind = '3'
    output_path = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/'

    path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"

    path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"

    path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"

    pairs, contigs, id_contig, longest_contig, longest_contig_name = get_contigs_and_pairs(path_layout, path_lens,
                                                                                           path_pairs,
                                                                                           min_len=0,
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
            for (k, pair) in enumerate(pairs):
                if distances[k] < 100_000:
                    lengths.append(distances[k])

            plt.hist(lengths)
            print(f'\tdistance mean: {np.array(lengths).mean()}')
            plt.title(f'Orientation {[contig.o for contig in contigs]}, {correct_number / len(contigs) * 100}%')
            plt.xlabel('distance between reads')
            plt.ylabel(f'number of reads')
            plt.savefig(f'{output_path}/plots_check/distances{[contig.o for contig in contigs]}.png')
            plt.clf()

    min_left = 10 ** 10
    max_left = 0
    min_right = 10 ** 10
    max_right = 0
    pos_left = []
    pos_right = []
    for (k, pair) in enumerate(pairs):
        min_left = min(min_left, pair[1])
        max_left = max(max_left, pair[1])
        min_right = min(min_right, pair[3])
        max_right = max(max_right, pair[3])
        pos_left.append(pair[1])
        pos_right.append(pair[3])

    print(
        f'\tmin read pos in left contig: {min_left}, max read pos in left contig: {max_left} ({contigs[0].length - max_left}), contig total length: {contigs[0].length}')
    print(
        f'\tmin read pos in right contig: {min_right}, max read pos in right contig: {max_right} ({contigs[1].length - max_right}), contig total length: {contigs[1].length}')

    plt.hist(pos_left)
    plt.title(f'Orientation {[contig.o for contig in contigs]}, {correct_number / len(contigs) * 100}%')
    plt.xlabel('read position in the first contig')
    plt.ylabel(f'number of reads')
    plt.savefig(f'{output_path}/positions_0.png')
    plt.clf()

    plt.hist(pos_right)
    plt.title(f'Orientation {[contig.o for contig in contigs]}, {correct_number / len(contigs) * 100}%')
    plt.xlabel('read position in the second contig')
    plt.ylabel(f'number of reads')
    plt.savefig(f'{output_path}/positions_1.png')
    plt.clf()

    check_reads(path_pairs, output_path)
