import numpy as np
from orientation.tools import get_orientation
from tools.load import check_reads, get_contigs_and_pairs
from tools.tools import get_distance
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
    print(
        f"{[contig.o for contig in contigs]}: {correct_number}/{len(contigs)} ({correct_number / len(contigs) * 100}%)")

    for i in range(2):
        for j in range(2):
            get_orientation([i, j], pairs, contigs)
            correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]] for contig in contigs]).sum()
            print(
                f"{[contig.o for contig in contigs]}: {correct_number}/{len(contigs)} ({correct_number / len(contigs) * 100}%)")

    distances = get_distance(pairs, contigs)
    lengths = []
    for (i, pair) in enumerate(pairs):
        if pair[0] == chr_ind or pair[2] == chr_ind:
            lengths.append(distances[i])

    plt.hist(lengths)
    plt.xlabel('distance between reads')
    plt.ylabel(f'number of reads')
    plt.savefig(f'{output_path}/plots_check/distances.png')
