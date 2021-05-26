import matplotlib.pyplot as plt

from ordering.model import MCMC
from ordering.tools import shuffle_ordering
from utils.load import get_contigs_and_pairs
from utils.prob import density
from utils.tools import log_likelihood

if __name__ == "__main__":
    print("Start!")
    min_contig_length = 100_000
    min_contig_length_name = '100k'

    path_to_output = f'/lustre/groups/cbi/Users/aeliseev/aivanova/data/{min_contig_length_name}'

    chromosomes = ['1', '2', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13',
                   '14', '15', '16', '17', '18', '19', '20', '21', '22']

    for chr_ind in chromosomes:
        print(f'Chromosome {chr_ind}')

        path_layout = "/GWSPH/groups/cbi/Users/pavdeyev/HiCProject/layouts/chr" + chr_ind + ".layout.txt"

        path_lens = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/contig_length/contig.length." + chr_ind + ".txt"

        path_pairs = "/lustre/groups/cbi/Users/aeliseev/aivanova/data/pairs/chr_pairs" + chr_ind + ".txt"

        # longest_contig
        pairs, contigs, id_contig, longest_contig, longest_contig_name = get_contigs_and_pairs(path_layout, path_lens,
                                                                                               path_pairs,
                                                                                               min_len=min_contig_length,
                                                                                               long_contig=True)
        correct_contigs = [contig.pos for contig in contigs]

        print("Estimation of density...")
        P, f = density(longest_contig)
        print("Estimation of density is done")

        print("MCMC is running...")

        correct_log_likelihood = log_likelihood(pairs, contigs, P)

        accuracy = 0
        n_it = 10
        for it in range(n_it):
            it_accuracy = 0

            shuffle_ordering(pairs, contigs, P, n_iterations=len(contigs) * 5)

            log_likelihood_arr = MCMC(pairs, contigs, P, number_it=500)
            print("Have found follow ordering:", [contigs[i].pos for i in range(len(contigs))])

            print(log_likelihood_arr)

            print('Plots...')
            plt.plot(log_likelihood_arr, label=f'{len(contigs)} contigs')
            plt.plot([correct_log_likelihood for _ in range(len(log_likelihood_arr))], label='correct order')
            plt.xlabel('iteration number')
            plt.ylabel('log_likelihood')
            plt.legend()
            plt.title(f'Ordering log likelihood')
            plt.savefig(f'{path_to_output}/plots_ordering/log_likelihood_chr{chr_ind}.png')
            plt.clf()

            print("Count correctness...")
            with open(f"{path_to_output}/final_ordering/final_chr{chr_ind}.layout.txt", "w") as file:
                file.write(",".join([str(contigs[i].pos) for i in range(len(contigs))]))

            for i in range(len(contigs) - 1):
                if contigs[i].pos + 1 == contigs[i + 1].pos:
                    it_accuracy += 1
            it_accuracy /= len(contigs) - 1
            accuracy += it_accuracy
            print(f'it: {it}\taccuracy: {it_accuracy}')

        accuracy /= n_it
        print(f'n_it: {n_it}\taccuracy: {accuracy}')
