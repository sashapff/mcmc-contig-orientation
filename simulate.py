from scipy.stats import expon
from ordering.model import MCMC
from ordering.tools import shuffle_ordering
from orientation.tools import simulation
from tools.load import get_contigs_and_pairs
from tools.prob import toy_density
import matplotlib.pyplot as plt

if __name__ == "__main__":
    simulation(max_len=1000, n_reads=10_000, n_contigs=5, p_distr=expon, output_path='../data_sim/')

    # print("Start!")
    # path_to_output = f'/Users/alexandra/bioinf/mcmc/ordering_sim/'
    # chromosomes = ['Sim']
    #
    # for chr_ind in chromosomes:
    #     print(f'Chromosome {chr_ind}')
    #
    #     path_layout = "/Users/alexandra/bioinf/mcmc/data/simulation.layout.txt"
    #
    #     path_lens = "/Users/alexandra/bioinf/mcmc/data/simulation.lens.tsv"
    #
    #     path_pairs = "/Users/alexandra/bioinf/mcmc/data/simulation.pairs.txt"
    #
    #     pairs, contigs, id_contig, longest_contig, longest_contig_name, in_contigs = get_contigs_and_pairs(path_layout,
    #                                                                                                        path_lens,
    #                                                                                                        path_pairs,
    #                                                                                                        long_contig=True,
    #                                                                                                        all_contigs=True,
    #                                                                                                        min_len=0)
    #
    #     correct_contigs = [contig.pos for contig in contigs]
    #
    #     print("Estimation of density...")
    #     P, f = toy_density(longest_contig)
    #     print("Estimation of density is done")
    #
    #     print("MCMC is running...")
    #
    #     shuffle_ordering(pairs, contigs, P, n_iterations=len(contigs) * 5)
    #
    #     log_likelihood_arr = MCMC(pairs, contigs, P, number_it=1000)
    #     print("Have found follow ordering:", [contigs[i].pos for i in range(len(contigs))])
    #
    #     plt.clf()
    #     plt.plot(log_likelihood_arr, label=f'{len(contigs)} contigs')
    #     plt.xlabel('iteration number')
    #     plt.ylabel('log_likelihood')
    #     plt.legend()
    #     plt.title(f'Ordering log likelihood')
    #     plt.savefig(f'{path_to_output}/plots_ordering/log_likelihood_chr{chr_ind}.png')
    #
    #     with open(f"{path_to_output}/final_ordering/final_chr{chr_ind}.layout.txt", "w") as file:
    #         file.write(",".join([str(contigs[i].pos) for i in range(len(contigs))]))