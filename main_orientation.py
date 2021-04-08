from tools.load import get_contigs_and_pairs
from orientation.model import MCMC
from orientation.tools import get_orientation
from tools.prob import density, toy_density
from tools.tools import get_longest_contig
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":
    print("Start!")

    correct_total = 0
    contigs_total = 0
    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13',
                   '14', '15', '16', '17', '18', '19', '20', '21', '22']

    # chromosomes = ['1']

    pairs_arr, contigs_arr, id_contig_arr, longest_contig_pairs_arr, longest_contig_arr, correct_contigs_arr = [], [], [], [], [], []

    min_contig_length = 50_000
    min_contig_length_name = '50k'

    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/{min_contig_length_name}'
    # path_to_output = f'/Users/alexandra/bioinf/mcmc/output/{min_contig_length_name}'

    with open(f"{path_to_output}/stat.txt", "w") as file:
        pass

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

        # path_layout = "/Users/alexandra/bioinf/mcmc/data/chr1.layout.txt"
        # # path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"
        #
        # path_lens = "/Users/alexandra/bioinf/mcmc/data/comp18_lens.tsv"
        # # path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"
        #
        # path_pairs = "/Users/alexandra/bioinf/mcmc/data/pairs18.txt"
        # # path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"
        #
        # # longest_contig
        # pairs, contigs, id_contig, longest_contig_pairs, longest_contig = get_contigs_and_pairs(path_layout, path_lens,
        #                                                                                         path_pairs,
        #                                                                                         long_contig=True)
        # # pairs, contigs, id_contig, longest_contig_pairs, longest_contig, in_contigs \
        # #     = get_contigs_and_pairs(path_layout,
        # #                             path_lens,
        # #                             path_pairs,
        # #                             long_contig=True,
        # #                             all_contigs=True,
        # #                             min_len=0)

        correct_contigs = [contig.o for contig in contigs]

        pairs_arr.append(pairs)
        contigs_arr.append(contigs)
        id_contig_arr.append(id_contig)
        longest_contig_pairs_arr.append(longest_contig_pairs)
        longest_contig_arr.append(longest_contig)
        correct_contigs_arr.append(correct_contigs)

    print("Estimation of density...")

    P, f = density(get_longest_contig(longest_contig_pairs_arr, longest_contig_arr))
    # P, f = toy_density(get_longest_contig(longest_contig_pairs_arr, longest_contig_arr))
    print("Estimation of density is done")

    for (j, chr_ind) in enumerate(chromosomes):
        pairs, contigs, id_contig, correct_contigs \
            = pairs_arr[j], contigs_arr[j], id_contig_arr[j], correct_contigs_arr[j]

        print("MCMC is running...")
        get_orientation([0 for i in range(len(contigs))], pairs, contigs)
        accuracy_arr, log_likelihood_arr = MCMC(pairs, contigs, P, number_it=100, n_chains=1,
                                                correct_contigs=correct_contigs)
        print("Have found follow orientation:", [contigs[i].o for i in range(len(contigs))])

        with open(f"{path_to_output}/final_chr{chr_ind}.layout.txt", "w") as file:
            sign = lambda x: "+" if x == 1 else "-"
            file.write(",".join([contig.name + sign(contig.o) for contig in contigs]))
        print("Result has been saved")

        print("Count correctness...")
        correct_number = np.array([contig.o == correct_contigs[id_contig[contig.name]] for contig in contigs]).sum()
        print(
            f"{correct_number}/{len(contigs)} contigs were oriented correctly ({correct_number / len(contigs) * 100}%)")

        with open(f"{path_to_output}/stat.txt", "a") as file:
            file.write(f"CHROMOSOME {chr_ind}\n")
            file.write(
                f"{correct_number}/{len(contigs)} contigs were oriented correctly ({correct_number / len(contigs) * 100}%)\n")

        correct_total += correct_number
        contigs_total += len(contigs)

        plt.clf()
        plt.plot(accuracy_arr, label=f'{len(contigs)} contigs')
        plt.xlabel('iteration number')
        plt.ylabel('accuracy')
        plt.savefig(f'{path_to_output}/plots/accuracy_chr{chr_ind}.png')

        plt.clf()
        plt.plot(log_likelihood_arr, label=f'{len(contigs)} contigs')
        plt.xlabel('iteration number')
        plt.ylabel('log_likelihood')
        plt.savefig(f'{path_to_output}/plots/log_likelihood_chr{chr_ind}.png')

    with open(f"{path_to_output}/stat.txt", "a") as file:
        file.write(f"TOTAL ACCURACY\n")
        file.write(
            f"{correct_total}/{contigs_total} contigs were oriented correctly ({correct_total / contigs_total * 100}%)\n")
