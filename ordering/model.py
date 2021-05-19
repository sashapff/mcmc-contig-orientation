import numpy as np
import random
from tqdm import tqdm
from ordering.tools import shuffle_ordering, change_position, swap, swap_log_likelihood, \
    swap_neigh, change_position_log_likelihood
from tools.tools import log_likelihood


def MCMC(pairs, contigs, P, number_it=500):
    """
    Changing orientation for better likelihood
    :param pairs: array of reads
    :param contigs: list of contigs
    :param P: density
    :param number_it: number of iterations
    :return None, change on place
    """
    lk_old = log_likelihood(pairs, contigs, P)

    print('BEFORE', lk_old, [contig.pos for contig in contigs])

    log_likelihood_arr = []

    for _ in tqdm(range(number_it)):
        number_contig_1 = np.random.randint(0, len(contigs))
        number_contig_2 = np.random.randint(0, len(contigs))

        last_contig_pos = contigs[number_contig_1].pos

        # lk_new = swap_log_likelihood(lk_old, number_contig_1, number_contig_2, pairs, contigs, P)
        # lk_new = change_position_log_likelihood(lk_old, number_contig_1, number_contig_2, pairs, contigs, P)
        change_position(number_contig_1, number_contig_2, pairs, contigs)
        lk_new = log_likelihood(pairs, contigs, P)

        if random.random() > np.exp(lk_new - lk_old):
            # Decline
            # print('decline')
            # swap(number_contig_1, number_contig_2, pairs, contigs)
            change_position(number_contig_1, last_contig_pos, pairs, contigs)
        else:
            # Accept
            # print('accept')
            lk_old = lk_new

        log_likelihood_arr.append(lk_old)
        # print(new_contigs, [contig.pos for contig in contigs], lk_old, lk_new)

    # assert new_contigs == [contig.pos for contig in contigs]
    # get_ordering(new_contigs, pairs, contigs)

    return log_likelihood_arr
