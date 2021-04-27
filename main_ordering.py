from tools.load import get_contigs_and_pairs
from ordering.model import MCMC
from ordering.tools import shuffle_ordering, change_position, change_position_log_likelihood, swap, swap_log_likelihood
from tools.prob import density, toy_density
import numpy as np
import matplotlib.pyplot as plt

from tools.tools import log_likelihood

if __name__ == "__main__":
    print("Start!")
    # chr_ind = sys.argv[1]
    # Example
    #
    # path_layout = "/Users/alexandra/bioinf/mcmc/data/chr1.layout.txt"
    # # path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"
    # #
    # path_lens = "/Users/alexandra/bioinf/mcmc/data/comp18_lens.tsv"
    # # path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"
    # #
    # path_pairs = "/Users/alexandra/bioinf/mcmc/data/pairs18.txt"
    # # path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"
    #
    # path_to_output = "/Users/alexandra/bioinf/mcmc/output"

    # longest_contig
    # pairs, contigs, id_contig, longest_contig, longest_contig_name = get_contigs_and_pairs(path_layout, path_lens, path_pairs, long_contig=True)
    # pairs, contigs, id_contig, longest_contig, longest_contig_name, in_contigs = get_contigs_and_pairs(path_layout,
    #                                                                                                    path_lens,
    #                                                                                                    path_pairs,
    #                                                                                                    long_contig=True,
    #                                                                                                    all_contigs=True,
    #                                                                                                    min_len=0)

    min_contig_length = 100_000
    min_contig_length_name = '100k'

    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/{min_contig_length_name}'
    # path_to_output = f'/Users/alexandra/bioinf/mcmc/output/{min_contig_length_name}'

    # chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13',
    #                '14', '15', '16', '17', '18', '19', '20', '21', '22']
    chromosomes = ['1']

    for chr_ind in chromosomes:
        print(f'Chromosome {chr_ind}')

        path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"

        # if chr_ind == '9':
        # path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr9.partial_layout.txt"

        path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"
        path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"

        # longest_contig
        pairs, contigs, id_contig, longest_contig_pairs, longest_contig = get_contigs_and_pairs(path_layout, path_lens,
                                                                                                path_pairs,
                                                                                                long_contig=True,
                                                                                                min_len=min_contig_length)

        correct_contigs = [contig.pos for contig in contigs]

        print("Estimation of density...")
        P, f = density(longest_contig)
        # P, f = toy_density(longest_contig)
        print("Estimation of density is done")

        print("MCMC is running...")

        shuffle_ordering(pairs, contigs, P, n_iterations=len(contigs) * 5)

        log_likelihood_arr = MCMC(pairs, contigs, P, number_it=100)
        print("Have found follow ordering:", [contigs[i].pos for i in range(len(contigs))])

        plt.clf()
        plt.plot(log_likelihood_arr, label=f'{len(contigs)} contigs')
        plt.xlabel('iteration number')
        plt.ylabel('log_likelihood')
        plt.legend()
        plt.title(f'Ordering log likelihood')
        plt.savefig(f'{path_to_output}/plots/log_likelihood.png')

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
