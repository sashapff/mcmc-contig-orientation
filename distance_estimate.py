from tools.load import get_contigs_and_pairs
from orientation.model import MCMC
from orientation.tools import get_orientation
from tools.prob import density, toy_density
from tools.tools import get_longest_contig, filter_pairs
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("Start!")

    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/distance'
    # path_to_output = f'/Users/alexandra/bioinf/mcmc/output/{min_contig_length_name}'

    chr_ind = '1'
    # path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"
    # path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"
    # path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"

    # path_layout = "/Users/alexandra/bioinf/mcmc/data/chr1.layout.txt"
    path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"

    # path_lens = "/Users/alexandra/bioinf/mcmc/data/comp18_lens.tsv"
    path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"

    # path_pairs = "/Users/alexandra/bioinf/mcmc/data/pairs18.txt"
    path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"

    # longest_contig
    pairs, contigs, id_contig, longest_contig_pairs, longest_contig_pairs = get_contigs_and_pairs(path_layout, path_lens,
                                                                                            path_pairs,
                                                                                            long_contig=True,
                                                                                            min_len=0)

    P = toy_density(longest_contig_pairs)

    d = longest_contig_pairs.length / 10
    left = longest_contig_pairs.length // 2 - d // 2
    right = left + d

    filtered_pairs = filter_pairs(longest_contig_pairs, left, right)
    print(filtered_pairs)




